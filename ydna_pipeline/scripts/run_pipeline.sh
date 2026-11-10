#!/usr/bin/env bash
# One-command reproduction. Dependencies: samtools>=1.10, bwa (0.7.x), python3 (stdlib only), curl.
#   sudo apt-get install -y samtools bwa curl
# Network: ftp.1000genomes.ebi.ac.uk, hgdownload.soe.ucsc.edu (HTTPS).
# Disk: ~600 MB in ydna_pipeline/ref (chrX+chrY FASTA + bwa index). No full-genome download.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results logs

bash scripts/00_setup_reference.sh ref 2>&1 | tee logs/00_setup_reference.log
[[ -s ref/chrXY_analysisset.fa.bwt ]] || { cat ref/chrY_analysisset.fa ref/chrX_analysisset.fa > ref/chrXY_analysisset.fa; bwa index ref/chrXY_analysisset.fa; }

# Local cache only: never fall back to the slow ENA refget server
export REF_PATH="$PWD/ref/cache/%2s/%2s/%s" REF_CACHE="$PWD/ref/cache/%2s/%2s/%s"

python3 scripts/01_annotate_sites.py --sites config/private_variants.tsv config/tree_markers.tsv \
  --refdir ref --out results/site_annotation.csv

python3 scripts/02_cram_site_evidence.py --samples config/samples.tsv \
  --roles "${ROLES:-target_P1,target_P2,target_P3}" \
  --sites config/private_variants.tsv config/tree_markers.tsv \
  --ref ref/chrY_analysisset.fa --mappability results/site_annotation.csv \
  --out results/jola_site_evidence.csv --reads-out results/jola_per_read_evidence.tsv \
  2>&1 | tee logs/02_cram_site_evidence.log

python3 scripts/03_summarize.py results/jola_site_evidence.csv results/site_annotation.csv > results/evidence_matrix.md

# Reliability strengthening (published genome-wide tracks, strict filter, YBrowse-named tree markers)
python3 scripts/01_annotate_sites.py --sites config/private_variants.tsv config/tree_markers.tsv --refdir ref --out results/site_annotation.csv
python3 scripts/04_genomewide_mappability.py --sites config/private_variants.tsv config/tree_markers.tsv --out results/genomewide_mappability.csv
python3 scripts/02_cram_site_evidence.py --strict --samples config/samples.tsv --roles "${ROLES:-target_P1,target_P2,target_P3}" \
  --sites config/private_variants.tsv config/tree_markers.tsv --ref ref/chrY_analysisset.fa \
  --mappability results/site_annotation.csv --out results/jola_site_evidence_strict.csv \
  --reads-out results/jola_per_read_evidence_strict.tsv 2>&1 | tee logs/02_cram_site_evidence_strict.log
python3 scripts/06_tree_concordance.py results/jola_site_evidence_strict.csv config/tree_markers.tsv > results/tree_concordance.md
python3 scripts/05_rank_private_snps.py
