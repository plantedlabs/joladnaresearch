#!/usr/bin/env bash
# =============================================================================
# 00_setup_reference.sh: fetch ONLY chrY (and chrX, for X/Y homology checks)
# from the exact reference FASTA the GGVP CRAMs were encoded against, using
# HTTP byte-range requests. Verify MD5 against the CRAM @SQ M5 tag. Populate a
# local htslib REF_CACHE so samtools never falls back to the slow ENA server.
#
# Why: the ENA refget endpoint (www.ebi.ac.uk/ena/cram/md5/<md5>) serves at
# ~85 KB/s, so htslib needs ~11 min to fetch the 57 MB chrY before it can
# decode a single read. That caused the 45-60 s timeouts (exit 124).
#
# Downloads: ~58 MB (chrY) + ~158 MB (chrX) + ~13 MB UCSC annotation.
# Usage:  bash 00_setup_reference.sh [REF_DIR]   (default: ../ref)
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REF_DIR="${1:-$HERE/../ref}"
mkdir -p "$REF_DIR"; cd "$REF_DIR"

FA_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CURL="curl -sS --fail --retry 4 --retry-delay 2 --connect-timeout 20"

# Expected MD5s, from the CRAM header (@SQ M5) of SC_GMJOL5309804
declare -A EXPECT_MD5=( [chrY]=ce3e31103314a704255f3cd90369ecce [chrX]=2b3a55ff7f58eb308420c8a9b11cac50 )

[[ -s GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ]] || \
  $CURL --max-time 60 -o GRCh38_full_analysis_set_plus_decoy_hla.fa.fai "$FA_URL.fai"

fetch_contig () {
  local ctg=$1 out="${1}_analysisset.fa"
  if [[ -s "$out" && -s "$out.md5ok" ]]; then echo "[ok] $out already verified"; return; fi
  read -r _ len off lb lw < <(awk -v c="$ctg" '$1==c' GRCh38_full_analysis_set_plus_decoy_hla.fa.fai)
  local nb=$(( len/lb*lw + (len%lb>0 ? len%lb+1 : 0) ))
  echo "[fetch] $ctg len=$len bytes=$nb offset=$off"
  { echo ">$ctg"; $CURL --max-time 900 -r "$off-$((off+nb-1))" "$FA_URL"; } > "$out.tmp"
  # MD5 per SAM spec: uppercase sequence, no newlines
  local md5; md5=$(grep -v '^>' "$out.tmp" | tr -d '\n' | tr a-z A-Z | md5sum | cut -d' ' -f1)
  if [[ "$md5" != "${EXPECT_MD5[$ctg]}" ]]; then
     echo "[FAIL] $ctg md5 $md5 != ${EXPECT_MD5[$ctg]}"; exit 1; fi
  mv "$out.tmp" "$out"; samtools faidx "$out"; echo "$md5" > "$out.md5ok"
  # htslib REF_CACHE layout: %2s/%2s/%s = md5[0:2]/md5[2:4]/md5[4:]
  mkdir -p "cache/${md5:0:2}/${md5:2:2}"
  grep -v '^>' "$out" | tr -d '\n' | tr a-z A-Z > "cache/${md5:0:2}/${md5:2:2}/${md5:4}"
  echo "[ok] $ctg md5 verified: $md5"
}
fetch_contig chrY
fetch_contig chrX

# UCSC hg38 chrY (soft-masked: lowercase = RepeatMasker/TRF) and segmental dups
[[ -s ucsc_hg38_chrY.fa ]] || { $CURL --max-time 300 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz | gunzip -c > ucsc_hg38_chrY.fa; samtools faidx ucsc_hg38_chrY.fa; }
[[ -s genomicSuperDups.chrY.tsv ]] || $CURL --max-time 300 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
   | gunzip -c | awk -F'\t' '$2=="chrY"' > genomicSuperDups.chrY.tsv

echo
echo "Reference ready. Use:"
echo "  export REF_PATH=$REF_DIR/cache/%2s/%2s/%s REF_CACHE=$REF_DIR/cache/%2s/%2s/%s"
echo "  samtools view -T $REF_DIR/chrY_analysisset.fa <CRAM_URL> chrY:POS-POS"
