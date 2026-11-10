# ydna_pipeline: chrY private-variant comparison against public African CRAMs

Reproducible pipeline comparing the 19 Big Y private variants of YF133256 (E-FTF75935 / YFull E-Z15174*)
with GGVP Jola men, using indexed remote access to the CRAMs. No full-genome or full-CRAM download is needed.

**Start here:** [`docs/REPORT.md`](docs/REPORT.md) (findings), [`results/evidence_matrix.md`](results/evidence_matrix.md),
[`docs/TECHNICAL_LOG.md`](docs/TECHNICAL_LOG.md), [`docs/EXPANSION_PLAN.md`](docs/EXPANSION_PLAN.md).

## Run
```bash
sudo apt-get install -y samtools bwa curl python3
bash scripts/run_pipeline.sh            # 3 target Jola samples (~10 min, ~600 MB ref/ on disk)
ROLES=control bash scripts/run_pipeline.sh   # optional negative-control Jola men
```

## The CRAM timeout fix (TL;DR)
htslib was fetching chrY from the ENA refget server at ~85 KB/s. `00_setup_reference.sh` instead byte-range-fetches
chrY from the exact 1000G analysis-set FASTA (5 s), verifies the MD5 against the CRAM header
(`ce3e31103314a704255f3cd90369ecce`), and builds a local `REF_CACHE`. Then:
```bash
export REF_PATH=$PWD/ref/cache/%2s/%2s/%s REF_CACHE=$PWD/ref/cache/%2s/%2s/%s
samtools view -T ref/chrY_analysisset.fa "$CRAM_URL" chrY:7013476-7013476   # ~2-8 s
```

## Layout
| Path | Content |
|---|---|
| `config/` | sites (19 private, 4 branch markers) and samples (CRAM URLs from the GGVP alignment index) |
| `scripts/00_setup_reference.sh` | chrY/chrX fetch, MD5 verification, htslib cache, UCSC repeat/segdup annotation |
| `scripts/01_annotate_sites.py` | reference check, repeats, segdups, **preliminary** chrX/chrY 100-mer mappability |
| `scripts/02_cram_site_evidence.py` | per-read evidence, filters, haploid classification |
| `scripts/03_summarize.py` | Markdown evidence matrix |
| `results/jola_site_evidence.csv` | main machine-readable table (one row per sample × site) |
| `results/jola_per_read_evidence.tsv` | every read examined, with its filter status |
| `results/site_annotation.csv` | site-level QC |
| `ref/` | (git-ignored) reference slices and bwa index |
