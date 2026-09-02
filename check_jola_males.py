#!/usr/bin/env python3
"""
Jola Y-DNA Haplogroup Scanner
==============================
Checks 40 unchecked Jola males from the Gambian Genome Variation Project
for presence of E-Z15174 haplogroup SNPs that match Kelland Drumgoole's
confirmed paternal lineage.

Requires: samtools installed and accessible in PATH
Run in WSL2 (Ubuntu) or Linux terminal.

Author: Kelland Drumgoole Research Project
Date: 2026
"""

import subprocess
import sys
import os
from datetime import datetime

# ─────────────────────────────────────────────────────────────────────────────
# KEY SNP POSITIONS (GRCh38 / hg38)
# These are the positions that define Kelland's Jola haplogroup E-Z15174
# ─────────────────────────────────────────────────────────────────────────────
SNPS = {
    "Z15174": "chrY:4202358-4202358",
    "Z15234": "chrY:8414659-8414659",
    "Z15231": "chrY:14719978-14719978",
    "Z15232": "chrY:14765245-14765245",
}

# Reference genome (GRCh38) - required to decode CRAM files
REFERENCE = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
)

# FTP base path for JOLA samples
FTP_BASE = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "gambian_genome_variation_project/data/JOLA"
)

# CRAM filename pattern
CRAM_PATTERN = "{sample}.alt_bwamem_GRCh38DH.20151208.JOLA.gambian_lowcov.cram"

# ─────────────────────────────────────────────────────────────────────────────
# 40 UNCHECKED MALES (identified from IGSR GWJ manifest)
# Confirmed matches already on YFull tree are EXCLUDED
# ─────────────────────────────────────────────────────────────────────────────
UNCHECKED_MALES = [
    "SC_GMJOL5309805",
    "SC_GMJOL5309827",
    "SC_GMJOL5309828",
    "SC_GMJOL5309844",
    "SC_GMJOL5309852",
    "SC_GMJOL5309868",
    "SC_GMJOL5309875",
    "SC_GMJOL5309876",
    "SC_GMJOL5309897",
    "SC_GMJOL5309903",
    "SC_GMJOL5309904",
    "SC_GMJOL5309905",
    "SC_GMJOL5309919",
    "SC_GMJOL5309921",
    "SC_GMJOL5309925",
    "SC_GMJOL5309926",
    "SC_GMJOL5309927",
    "SC_GMJOL5309928",
    "SC_GMJOL5309937",
    "SC_GMJOL5309940",
    "SC_GMJOL5309941",
    "SC_GMJOL5309942",
    "SC_GMJOL5309944",
    "SC_GMJOL5309945",
    "SC_GMJOL5309948",
    "SC_GMJOL5309949",
    "SC_GMJOL5309950",
    "SC_GMJOL5309951",
    "SC_GMJOL5309952",
    "SC_GMJOL5309953",
    "SC_GMJOL5309963",
    "SC_GMJOL5309967",
    "SC_GMJOL5309968",
    "SC_GMJOL5309971",
    "SC_GMJOL5309972",
    "SC_GMJOL5309973",
    "SC_GMJOL5309974",
    "SC_GMJOL5309975",
    "SC_GMJOL5309976",
    "SC_GMJOL5309977",
]

# Already confirmed on YFull tree (for reference)
CONFIRMED_MATCHES = [
    "SC_GMJOL5309804",  # E-Z15196/dyo
    "SC_GMJOL5309829",  # E-Z15196/dyo
    "SC_GMJOL5309851",  # E-FT402211
]

# Unchecked fathers (also run these)
UNCHECKED_FATHERS = [
    "SC_GMJOL5309898",
    "SC_GMJOL5309899",
    "SC_GMJOL5309922",
]


def check_samtools():
    """Verify samtools is installed."""
    try:
        result = subprocess.run(
            ["samtools", "--version"],
            capture_output=True, timeout=10
        )
        version = result.stdout.decode("utf-8", errors="replace").split('\n')[0]
        print(f"✓ samtools found: {version}")
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("✗ samtools not found. Install with:")
        print("  sudo apt-get install -y samtools")
        return False


def get_pileup(sample_id, snp_name, region):
    """
    Run samtools mpileup for a specific sample and SNP position.
    Returns the base calls at that position, or None if failed.
    """
    cram_filename = CRAM_PATTERN.format(sample=sample_id)
    cram_url = f"{FTP_BASE}/{sample_id}/{cram_filename}"

    cmd = [
        "samtools", "mpileup",
        "-r", region,
        "-f", REFERENCE,
        cram_url,
        "--no-BAQ",
        "-q", "0",
        "-Q", "0",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout per sample per SNP
        )

        if result.returncode != 0:
            return None, "error"

        lines = [l for l in result.stdout.strip().split('\n') if l]
        if not lines:
            return None, "no_data"

        # Parse pileup: chrom, pos, ref, depth, bases, quals
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 5:
                ref_base = parts[2].upper()
                depth = int(parts[3])
                bases = parts[4].upper()

                # Clean bases - remove markers
                clean = ""
                i = 0
                while i < len(bases):
                    c = bases[i]
                    if c in [',', '.']:
                        clean += ref_base  # reference allele
                    elif c in ['A', 'T', 'G', 'C']:
                        clean += c
                    elif c in ['+', '-']:
                        # Skip indel
                        i += 1
                        num = ""
                        while i < len(bases) and bases[i].isdigit():
                            num += bases[i]
                            i += 1
                        if num:
                            i += int(num) - 1
                    elif c in ['^']:
                        i += 1  # skip mapping quality
                    i += 1

                # Count base frequencies
                from collections import Counter
                counts = Counter(clean)
                return {
                    "ref": ref_base,
                    "depth": depth,
                    "bases": dict(counts),
                    "raw": bases[:50]
                }, "ok"

        return None, "no_coverage"

    except subprocess.TimeoutExpired:
        return None, "timeout"
    except Exception as e:
        return None, f"exception: {e}"


def interpret_result(result, status, snp_name):
    """
    Interpret the pileup result.
    Returns: DERIVED, ANCESTRAL, LOW_COVERAGE, or ERROR
    """
    if status != "ok" or result is None:
        return status.upper()

    depth = result["depth"]
    if depth < 2:
        return "LOW_COVERAGE"

    ref = result["ref"]
    bases = result["bases"]
    total = sum(bases.values())

    if total == 0:
        return "NO_DATA"

    # Check if any non-reference allele is present at >20% frequency
    for base, count in bases.items():
        if base != ref and count / total > 0.20:
            return f"POSSIBLE_DERIVED ({base}:{count}/{total})"

    # All reads match reference
    return f"ANCESTRAL (ref={ref}, depth={depth})"


def run_analysis(samples_to_check, label="SAMPLES"):
    """Run the full haplogroup check for a list of samples."""
    results = {}

    print(f"\n{'='*60}")
    print(f"Checking {len(samples_to_check)} {label}")
    print(f"SNPs: {', '.join(SNPS.keys())}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}\n")

    for i, sample in enumerate(samples_to_check, 1):
        print(f"[{i}/{len(samples_to_check)}] {sample}")
        sample_results = {}

        for snp_name, region in SNPS.items():
            print(f"  Checking {snp_name} ({region})...", end=" ", flush=True)
            result, status = get_pileup(sample, snp_name, region)
            interpretation = interpret_result(result, status, snp_name)
            sample_results[snp_name] = interpretation
            print(interpretation)

        results[sample] = sample_results

        # Flag if any possible derived alleles found
        if any("POSSIBLE_DERIVED" in v for v in sample_results.values()):
            print(f"  *** POTENTIAL MATCH — {sample} shows derived alleles ***")

        print()

    return results


def write_report(results, filename="jola_haplogroup_results.txt"):
    """Write a formatted report of all results."""
    with open(filename, 'w') as f:
        f.write("JOLA Y-DNA HAPLOGROUP SCAN RESULTS\n")
        f.write("="*60 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Target haplogroup: E-Z15174 (Kelland Drumgoole lineage)\n")
        f.write(f"SNPs checked: {', '.join(SNPS.keys())}\n")
        f.write(f"Samples checked: {len(results)}\n\n")

        # Summary — potential matches first
        f.write("POTENTIAL MATCHES (possible derived alleles found):\n")
        f.write("-"*40 + "\n")
        matches_found = False
        for sample, snp_results in results.items():
            if any("POSSIBLE_DERIVED" in v for v in snp_results.values()):
                f.write(f"\n*** {sample} ***\n")
                for snp, result in snp_results.items():
                    f.write(f"  {snp}: {result}\n")
                matches_found = True

        if not matches_found:
            f.write("No potential matches found in this scan.\n")

        # Full results
        f.write("\n\nFULL RESULTS:\n")
        f.write("-"*40 + "\n")
        for sample, snp_results in results.items():
            f.write(f"\n{sample}:\n")
            for snp, result in snp_results.items():
                f.write(f"  {snp}: {result}\n")

    print(f"\nReport saved to: {filename}")
    return filename


def main():
    print("\n" + "="*60)
    print("JOLA Y-DNA HAPLOGROUP SCANNER")
    print("Drumgoole Family Research Project")
    print("Target: E-Z15174* (Jola-Fonyi lineage)")
    print("="*60)

    # Check samtools
    if not check_samtools():
        sys.exit(1)

    # Determine which samples to run
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        # Test mode — just run 3 samples
        print("\nTEST MODE — running 3 samples only")
        samples = UNCHECKED_MALES[:3]
    elif len(sys.argv) > 1 and sys.argv[1] == "--fathers":
        # Just check the 3 unchecked fathers
        samples = UNCHECKED_FATHERS
    elif len(sys.argv) > 1 and sys.argv[1] == "--all":
        # Everything including unchecked fathers
        samples = UNCHECKED_MALES + UNCHECKED_FATHERS
    else:
        # Default — all 40 unchecked males
        samples = UNCHECKED_MALES

    # Run analysis
    results = run_analysis(samples, "UNCHECKED JOLA MALES")

    # Write report
    report_file = write_report(results)

    # Print summary
    print("\n" + "="*60)
    print("SCAN COMPLETE")
    print("="*60)
    potential_matches = [
        s for s, r in results.items()
        if any("POSSIBLE_DERIVED" in v for v in r.values())
    ]
    print(f"Samples checked: {len(results)}")
    print(f"Potential matches found: {len(potential_matches)}")
    if potential_matches:
        print("\nPOTENTIAL NEW JOLA RELATIVES:")
        for m in potential_matches:
            print(f"  *** {m} ***")
    print(f"\nFull report: {report_file}")


if __name__ == "__main__":
    main()
