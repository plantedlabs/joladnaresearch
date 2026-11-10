#!/usr/bin/env bash
# Confirm FTF75935-derived candidates: strict evidence at 13 YBrowse-named tree markers
# (M132 ... Z15174 block ... Z5991/FT402211/Z15196/FTF75935) + 11 Tier-A private SNPs of YF133256.
# Usage: bash scripts/09_confirm_candidates.sh SAMPLE1,SAMPLE2,...
set -euo pipefail
cd "$(dirname "$0")/.."
export REF_PATH="$PWD/ref/cache/%2s/%2s/%s" REF_CACHE="$PWD/ref/cache/%2s/%2s/%s"
IFS=, read -ra S <<< "$1"
out=config/candidates.tsv; echo -e "#sample\tpopulation\trole\tyfull_branch\tcram_url" > $out
for s in "${S[@]}"; do awk -v s=$s -F'\t' '$1==s{print $1"\tGGVP "$2"\tcandidate\tunknown\t"$3}' config/ggvp_manifest.tsv >> $out; done
python3 scripts/02_cram_site_evidence.py --strict --samples $out --sites config/tree_markers.tsv config/tierA_private_variants.tsv \
  --ref ref/chrY_analysisset.fa --mappability results/site_annotation.csv \
  --out results/ftf75935_screen/candidate_confirmation.csv --reads-out results/ftf75935_screen/candidate_confirmation_reads.tsv
