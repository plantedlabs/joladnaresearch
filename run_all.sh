#!/bin/bash
# =============================================================================
# Jola Y-DNA Haplogroup Scanner — Full Dataset
# =============================================================================
# Screens all 40 unchecked male samples from the Gambian Genome Variation
# Project (PRJEB3252) for haplogroup E-Z15174 — the Jola-Fonyi lineage
# confirmed in Kelland Drumgoole's paternal ancestry.
#
# Requirements:
#   - samtools installed (sudo apt-get install samtools)
#   - Internet access to ftp.1000genomes.ebi.ac.uk
#   - chrY reference genome (downloaded by setup below)
#
# Usage:
#   chmod +x run_all.sh
#   ./run_all.sh
#
# Results saved to: all_results.txt
#
# Researcher: Kelland Drumgoole (YF133256)
# Haplogroup: E-FTF75935 / E-Z15174*
# Date: September 5, 2026
# =============================================================================

set -e

# Configuration
CRAM_BASE="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/data/JOLA"
CRAM_PATTERN=".alt_bwamem_GRCh38DH.20151208.JOLA.gambian_lowcov.cram"
REF="chrY.fa"
OUT="all_results.txt"

# SNP positions (GRCh38 / hg38)
# Format: "NAME:REGION:REF:DERIVED"
SNPS=(
    "Z15174:chrY:4202358-4202358:G:T"
    "Z15234:chrY:8414659-8414659:C:T"
    "Z15231:chrY:14719978-14719978:G:A"
    "Z15232:chrY:14765245-14765245:T:C"
)

# All 40 unchecked males from GGVP GWJ dataset
SAMPLES=(
    SC_GMJOL5309805 SC_GMJOL5309827 SC_GMJOL5309828
    SC_GMJOL5309844 SC_GMJOL5309852 SC_GMJOL5309868
    SC_GMJOL5309875 SC_GMJOL5309876 SC_GMJOL5309897
    SC_GMJOL5309903 SC_GMJOL5309904 SC_GMJOL5309905
    SC_GMJOL5309919 SC_GMJOL5309921 SC_GMJOL5309925
    SC_GMJOL5309926 SC_GMJOL5309927 SC_GMJOL5309928
    SC_GMJOL5309937 SC_GMJOL5309940 SC_GMJOL5309941
    SC_GMJOL5309942 SC_GMJOL5309944 SC_GMJOL5309945
    SC_GMJOL5309948 SC_GMJOL5309949 SC_GMJOL5309950
    SC_GMJOL5309951 SC_GMJOL5309952 SC_GMJOL5309953
    SC_GMJOL5309963 SC_GMJOL5309967 SC_GMJOL5309968
    SC_GMJOL5309971 SC_GMJOL5309972 SC_GMJOL5309973
    SC_GMJOL5309974 SC_GMJOL5309975 SC_GMJOL5309976
    SC_GMJOL5309977
)

# Already confirmed E-Z15174 carriers (not re-checked here)
# SC_GMJOL5309804 — E-Z15196, Jola-Fonyi (dyo)
# SC_GMJOL5309829 — E-Z15196, Jola-Fonyi (dyo)
# SC_GMJOL5309851 — E-FT402211

# =============================================================================
# SETUP
# =============================================================================

echo "=============================================="
echo "Jola Y-DNA Haplogroup Scanner"
echo "Target: E-Z15174 (Jola-Fonyi lineage)"
echo "=============================================="
echo ""

# Check samtools
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found. Install with:"
    echo "  sudo apt-get install -y samtools"
    exit 1
fi
echo "samtools version: $(samtools --version | head -1)"

# Download chrY reference if needed
if [ ! -f "$REF" ]; then
    echo ""
    echo "Downloading chrY reference genome (GRCh38)..."
    wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz" -O chrY.fa.gz
    gunzip chrY.fa.gz
    samtools faidx chrY.fa
    echo "chrY reference ready"
else
    echo "chrY reference already present"
fi

# Initialize output file
> "$OUT"

# =============================================================================
# SCAN
# =============================================================================

echo ""
echo "Starting scan of ${#SAMPLES[@]} Jola males..." | tee -a "$OUT"
echo "Target SNPs: Z15174, Z15234, Z15231, Z15232" | tee -a "$OUT"
echo "Started: $(date)" | tee -a "$OUT"
echo "==============================================" | tee -a "$OUT"

POSSIBLE_MATCHES=()

for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    num=$((i+1))
    total=${#SAMPLES[@]}

    echo "" | tee -a "$OUT"
    echo "[$num/$total] $sample" | tee -a "$OUT"

    CRAM="${CRAM_BASE}/${sample}/alignment/${sample}${CRAM_PATTERN}"
    derived_count=0

    for snp_def in "${SNPS[@]}"; do
        name=$(echo "$snp_def" | cut -d: -f1)
        region=$(echo "$snp_def" | cut -d: -f2-3)
        ref_allele=$(echo "$snp_def" | cut -d: -f4)
        derived_allele=$(echo "$snp_def" | cut -d: -f5)

        result=$(samtools mpileup -f "$REF" \
            -r "$region" \
            --no-BAQ -q 0 -Q 0 \
            "$CRAM" 2>/dev/null)

        if [ -z "$result" ]; then
            echo "  $name: NO_COVERAGE" | tee -a "$OUT"
        else
            # Check if derived allele is present
            bases=$(echo "$result" | awk '{print $5}' | tr '[:lower:]' '[:upper:]')
            derived_reads=$(echo "$bases" | grep -o "$derived_allele" | wc -l || true)
            total_reads=$(echo "$result" | awk '{print $4}')

            echo "  $name: $result" | tee -a "$OUT"

            if [ "$derived_reads" -gt 0 ] && [ "$total_reads" -gt 0 ]; then
                pct=$((derived_reads * 100 / total_reads))
                if [ "$pct" -ge 50 ]; then
                    echo "    *** DERIVED ALLELE FOUND ($derived_reads/$total_reads reads = ${pct}%) ***" | tee -a "$OUT"
                    ((derived_count++)) || true
                fi
            fi
        fi
    done

    if [ "$derived_count" -ge 2 ]; then
        echo "  *** POTENTIAL E-Z15174 MATCH — $sample ***" | tee -a "$OUT"
        POSSIBLE_MATCHES+=("$sample")
    fi
done

# =============================================================================
# SUMMARY
# =============================================================================

echo "" | tee -a "$OUT"
echo "==============================================" | tee -a "$OUT"
echo "SCAN COMPLETE: $(date)" | tee -a "$OUT"
echo "" | tee -a "$OUT"
echo "Samples checked: ${#SAMPLES[@]}" | tee -a "$OUT"
echo "Potential new E-Z15174 carriers: ${#POSSIBLE_MATCHES[@]}" | tee -a "$OUT"

if [ ${#POSSIBLE_MATCHES[@]} -gt 0 ]; then
    echo "" | tee -a "$OUT"
    echo "POTENTIAL NEW JOLA RELATIVES:" | tee -a "$OUT"
    for match in "${POSSIBLE_MATCHES[@]}"; do
        echo "  *** $match ***" | tee -a "$OUT"
    done
else
    echo "" | tee -a "$OUT"
    echo "No additional E-Z15174 carriers found." | tee -a "$OUT"
    echo "Confirmed carriers remain: SC_GMJOL5309804, SC_GMJOL5309829, SC_GMJOL5309851" | tee -a "$OUT"
    echo "E-Z15174 frequency in GGVP Jola dataset: 3/47 = 6.4%" | tee -a "$OUT"
fi

echo "" | tee -a "$OUT"
echo "Full results saved to: $OUT"
