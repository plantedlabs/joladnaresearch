# Technical Log: GGVP Jola chrY 19-variant comparison

Run date: 2026-11-10. Execution host: Claude Code cloud container (Ubuntu 24.04, Linux x86_64).
The same commands run unchanged on Ubuntu 24.04 under WSL2.

## Software
| Tool | Version | Source |
|---|---|---|
| samtools / htslib | 1.19.2 / 1.19 | Ubuntu apt `samtools` |
| bcftools | 1.19 | Ubuntu apt (installed, not used in final pipeline) |
| bwa | 0.7.17-r1188 | Ubuntu apt `bwa` |
| Python | 3.x, standard library only | system |
| curl | system | system |

The exact versions captured at run time are in `logs/versions.txt`.

## Data sources (all public, HTTPS)
| Item | URL / identifier |
|---|---|
| CRAMs | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/data/JOLA/<ID>/alignment/<ID>.alt_bwamem_GRCh38DH.20151208.JOLA.gambian_lowcov.cram` (paths for 5309851 / 5309829 confirmed in `gambian_genome_variation_project.GRCh38DH.alignment.index`, 395 lines; 95 FULA, 100 JOLA, 99 MANDINKA, 100 WOLOFF CRAMs) |
| Reference | `…/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa` (only the chrY and chrX byte ranges were fetched, via `.fai` offsets) |
| Repeat masking | UCSC hg38 `chromosomes/chrY.fa.gz` (soft-masked) |
| Segmental dups | UCSC hg38 `database/genomicSuperDups.txt.gz` (filtered to chrY) |
| Tree | YFull YTree v14.07.00 (02 Nov 2026), `https://www.yfull.com/tree/E-Z15174/` |

## Diagnosis of the CRAM decoding timeout (exit 124)

| Test | Result |
|---|---|
| `curl -I` on CRAM | 200, `Accept-Ranges: bytes`, Content-Length 12,724,790,752 |
| `samtools view -H` remote CRAM | 6.4 s; 3,366 @SQ lines |
| chrY @SQ | `LN:57227415 M5:ce3e31103314a704255f3cd90369ecce UR:ftp://…GRCh38_full_analysis_set_plus_decoy_hla.fa` |
| ENA refget `www.ebi.ac.uk/ena/cram/md5/ce3e3110…` | HTTP 200, but **1,711,618 bytes in 20 s (~85 KB/s)**; 57 MB would need ~11 min |
| Byte-range fetch of chrY slice from the 1000G FASTA (offset 3,074,345,836, 58,044,950 bytes) | **5.1 s**; MD5 of uppercase sequence = `ce3e31103314a704255f3cd90369ecce` ✔ |
| `samtools view -T chrY_analysisset.fa <CRAM> chrY:7013476-7013476` with `REF_PATH` disabled | **8.2 s**, 17 reads decoded |

**Conclusion.** Remote indexed access was never the problem. `samtools view -c` succeeds because counting
records needs no reference. Decoding bases needs the reference, and htslib was fetching all of chrY from
the ENA refget server at about 85 KB/s, so the 45–60 s timeouts always fired first. The fix is to
pass a local, MD5-verified chrY with `-T` and point `REF_PATH`/`REF_CACHE` at a local cache only.

Why not UCSC `chrY.fa`? The earlier repo scan used it with `-f`. It differs from the analysis set: UCSC
is soft-masked and its PARs are not N-masked, so its MD5 will not match `M5:` in the header. Decoding
mostly still works because htslib warns but continues, but it is not guaranteed for reads that
reference-diff across masked regions. The analysis-set slice is the exact reference the CRAMs use.

Environment note: the cloud sandbox initially blocked EBI and UCSC hosts. A personal environment
("joladnaresearch") was created that allows `ftp.1000genomes.ebi.ac.uk`, `www.ebi.ac.uk`,
`hgdownload.soe.ucsc.edu`, `www.yfull.com` and `discover.familytreedna.com`. On WSL2 this does not apply.

## Existing chrY VCFs
`release/20200217_biallelic_SNV/` contains chr1–22 and chrX only. **No chrY genotype calls are
published for GGVP**, so CRAM pileup is required.

## Pipeline steps
1. `00_setup_reference.sh` range-fetches chrY and chrX, verifies MD5 against the CRAM header
   (chrX `2b3a55ff7f58eb308420c8a9b11cac50`), builds `faidx`, populates `ref/cache/` in htslib
   md5 layout, and downloads the UCSC soft-masked chrY and chrY segdups.
2. `bwa index` on chrY+chrX (~213 Mb; a few minutes).
3. `01_annotate_sites.py` checks the reference base, repeat mask, segdups and PAR, and runs the
   **preliminary** chrX/chrY-only 100-mer mappability test (see the caveats in the script header).
4. `02_cram_site_evidence.py` runs per-read evidence at each site in each CRAM (see the QC rules below).
5. `03_summarize.py` builds the Markdown evidence matrix.

## QC decisions
* **Read exclusion:** flags 0x4, 0x100, 0x200, 0x400 (BioBamBam-marked duplicates) and 0x800, plus
  reads with a deletion or ref-skip at the site.
* **MAPQ < 20** and **base quality < 20** are counted but not used for calling. GGVP base qualities
  are GATK-BQSR-recalibrated and binned.
* Reads with the site within 5 bp of a read end, and soft-clipped reads, are counted and reported, not removed.
* **Haploid model:** Y is hemizygous, so a genuine site should be about 100% one allele. The ≥90%
  threshold allows for roughly 1 error read in 10.
* **Minimum 2 filtered reads** for any call. A single read is reported as INSUFFICIENT_COVERAGE
  (the base it shows is noted) and **never** called reference. Zero coverage is a no-call.
* **UNRELIABLE_MAPPING** applies when more than 50% of primary non-duplicate reads have MAPQ < 20,
  or when the preliminary mappability flag is `fail`.
* **Confidence:** high = ≥4 supporting reads on both strands and mappability `pass`;
  moderate = ≥3 supporting reads and not `caution`; otherwise low.
* **Error handling:** each query has a 180 s timeout with 2 retries and exponential backoff, and
  failures are recorded as `NO_DATA_ERROR`.

## Errors encountered during this session
* `curl` → 403 from the sandbox egress proxy for `ftp.1000genomes.ebi.ac.uk`; resolved by the network allowlist above.
* YFull `snp-list-search` only resolves Y/YP-series names, so Z-series coordinates could not be confirmed there.
  The four branch-marker coordinates come from the repo's `run_all.sh` and are labelled unverified.
  The CRAM results serve as an empirical check: if they are real E-Z15174 markers, all three Jola men must be derived.
* First run of `01_annotate_sites.py` used the wrong segdup column (alignment length instead of
  fracMatch) and counted short partial `bwa` alt alignments. Both were fixed: column 27, and ≥90 aligned bp.
