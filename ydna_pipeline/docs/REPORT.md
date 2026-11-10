# Comparing 19 private Y-chromosome variants of YF133256 with three GGVP Jola men

**Date:** 2026-11-10
**Data:** Gambian Genome Variation Project low-coverage GRCh38 CRAMs (ENA ERP001781)
**Tree reference:** YFull YTree v14.07.00 (02 Nov 2026)
**Files:** `results/jola_site_evidence.csv` (machine-readable, 69 rows: 57 private + 12 marker cells),
`results/evidence_matrix.md` (matrix), `results/site_annotation.csv` (site QC),
`results/jola_per_read_evidence.tsv` (every read examined).

---

## 1. Bottom line

* **None of the three Jola men carries any of your 19 private variants.** Across the 57
  sample × site cells there are **289 high-quality reads, and every one shows the reference
  (ancestral) allele.** There are **0** derived and **0** third-allele reads.
* **46/57 cells are positive ancestral calls** (23 high, 10 moderate, 13 low confidence). 5 are
  no-calls for insufficient coverage, and 6 are no-calls for unreliable mapping. No cell was called
  ancestral merely for lacking alternate reads: every ancestral call needs ≥2 concordant high-quality reads.
* **All three men are confirmed derived for the E-Z15174 block** at every marker with coverage
  (45 derived reads, 0 ancestral; 9/12 marker cells callable). This also empirically supports the
  marker coordinates taken from the earlier repo script.
* This is **exactly what the published tree predicts.** YFull places you on E-Z15174* and the
  three Jola men on the sister clade E-Z5991 (E-FT402211: 5309804 and 5309851; E-Z15196: 5309829).
  Your private variants arose on your own line **after** that split.
  **There is no evidence of a closer paternal relationship than YFull already shows**, a shared
  ancestor about 7,000 years ago (YFull TMRCA 95% CI 8,400–5,700 ybp).

## 2. What was tested

| Set | Sites | Purpose |
|---|---|---|
| Private variants (PV01–PV19) | 19 | Does any Jola man share a mutation from your terminal line? |
| E-Z15174 block markers (Z15174, Z15234, Z15231, Z15232) | 4 | Positive control: all three men should be derived. Also confirms coordinates. |

Samples: SC_GMJOL5309804 (E-FT402211), SC_GMJOL5309851 (E-FT402211), SC_GMJOL5309829 (E-Z15196).
The CRAM paths come from the official GGVP alignment index. Note that the repo `README.md` table
lists 5309804 as E-Z15196; YFull v14.07 places it on **E-FT402211**, with 5309851.

## 3. Evidence matrix (private variants)

Cells: `CALL ref/alt/other (filtered depth) MAPQ-median BQ-median`. Filters: primary,
non-duplicate, MAPQ ≥ 20, BQ ≥ 20. Full per-cell metrics (raw depth, excluded reads, strand
counts, low-MAPQ fraction, soft clipping, read-end proximity) are in the CSV.

| Site | GRCh38 | Ref>Der | Prelim. mappability | 5309804 (FT402211) | 5309851 (FT402211) | 5309829 (Z15196) |
|---|---|---|---|---|---|---|
| PV01 | 3339724 | G>C | **fail** (XTR) | BADMAP 3/0/0 MQ27 | BADMAP 3/0/0 MQ27 | BADMAP 1/0/0 |
| PV02 | 5039793 | C>T | caution (XTR) | ANC(low) 6/0/0 | ANC(low) 4/0/0 | ANC(low) 5/0/0 |
| PV03 | 5650649 | G>T | caution (XTR) | ANC(low) 6/0/0 | ANC(low) 2/0/0 | LOWCOV 1/0/0 |
| PV04 | 6220233 | T>C | caution (XTR) | ANC(low) 6/0/0 | LOWCOV 1/0/0 | ANC(low) 4/0/0 |
| PV05 | 6926993 | A>G | pass | ANC(high) 7/0/0 | ANC(mod) 3/0/0 | ANC(low) 2/0/0 |
| PV06 | 7013476 | T>C | pass | ANC(high) 17/0/0 | ANC(high) 7/0/0 | ANC(high) 7/0/0 |
| PV07 | 7758704 | C>A | pass | ANC(high) 7/0/0 | ANC(high) 6/0/0 | ANC(high) 11/0/0 |
| PV08 | 9023319 | C>A | caution (chrY paralog) | ANC(low) 9/0/0 | ANC(low) 6/0/0 | ANC(low) 5/0/0 |
| PV09 | 10820475 | A>T | pass* | ANC(mod) 6/0/0 | ANC(mod) 5/0/0 | ANC(mod) 5/0/0 |
| PV10 | 10973240 | C>T | pass* (satellite edge) | BADMAP 10/0/0 | BADMAP 3/0/0 | BADMAP 11/0/0 |
| PV11 | 14952809 | C>A | pass | ANC(high) 10/0/0 | ANC(mod) 3/0/0 | ANC(high) 5/0/0 |
| PV12 | 15252463 | T>C | pass | ANC(high) 5/0/0 | LOWCOV 1/0/0 | ANC(mod) 3/0/0 |
| PV13 | 15391518 | G>A | pass | ANC(high) 5/0/0 | ANC(mod) 3/0/0 | ANC(high) 5/0/0 |
| PV14 | 15936416 | G>A | pass | ANC(high) 9/0/0 | ANC(high) 9/0/0 | ANC(high) 7/0/0 |
| PV15 | 15990294 | C>A | pass | ANC(high) 8/0/0 | ANC(high) 4/0/0 | ANC(mod) 3/0/0 |
| PV16 | 18949120 | C>T | pass | ANC(high) 4/0/0 | ANC(high) 5/0/0 | LOWCOV 1/0/0 |
| PV17 | 21349919 | G>A | pass | ANC(mod) 3/0/0 | ANC(high) 4/0/0 | ANC(mod) 3/0/0 |
| PV18 | 21680474 | C>G | pass (segdup 0.965) | ANC(high) 4/0/0 | LOWCOV 1/0/0 | ANC(low) 2/0/0 |
| PV19 | 26572350 | C>T | pass (segdup 0.962) | ANC(high) 7/0/0 | ANC(low) 2/0/0 | ANC(high) 4/0/0 |

\* The preliminary test passed, but the empirical CRAM evidence shows many low-MAPQ reads (§5).

| Sample | Derived | Ancestral | Mixed | Insufficient cov. | Unreliable mapping |
|---|---|---|---|---|---|
| SC_GMJOL5309804 | 0 | 17 | 0 | 0 | 2 |
| SC_GMJOL5309851 | 0 | 14 | 0 | 3 | 2 |
| SC_GMJOL5309829 | 0 | 15 | 0 | 2 | 2 |

**Branch-marker control:** 5309804 is derived at 4/4 markers. 5309851 is derived at 2/2 covered
markers (no usable reads at Z15174/Z15234). 5309829 is derived at 3/3 covered markers (no reads at
Z15174). There were zero ancestral reads at any marker.

## 4. Site reliability (preliminary)

> **The mappability assessment is preliminary.** 100-mers tiled across each site were aligned
> with `bwa mem -a` to **chrX + chrY only**. This detects near-identical X/Y paralogs. It does
> **not** test autosomes, decoys or alt contigs, so it **does not establish genome-wide uniqueness.**

| Tier | Sites | Basis |
|---|---|---|
| **A: well-behaved** | PV05, PV06, PV07, PV11, PV12, PV13, PV14, PV15, PV16, PV17 | Preliminary pass; Jola reads MAPQ 60 with no low-MAPQ pile-up |
| **B: usable with caution** | PV09, PV18, PV19 | PV09: 38–50% of Jola reads at the site are MAPQ<20, all showing reference. PV18/PV19: inside segmental duplications (~96% identity) |
| **C: paralog risk** | PV02, PV03, PV04, PV08 | PV02–PV04 lie in the **X-transposed region** (98.8% identical to Xq21); PV08 has a near-identical chrY copy (~21.8 Mb). This matches your Big Y MAPQ<30 at PV02/PV04 |
| **D: unsuitable for phylogenetics** | PV01, PV10 | PV01: majority of 100-mers ambiguous between X and Y (Big Y MAPQ<30 too). PV10: edge of an `(ATTCC)n` satellite; 86–91% of Jola reads MAPQ<20 and mostly heavily clipped |

Your Big Y flagged three sites with MAPQ<30 (PV01, PV02, PV04). This analysis independently places
all three in the X-transposed region and identifies three additional weaker sites (PV03, PV08, PV10).
On a genome-wide basis, roughly **10–13 of the 19 variants** look suitable as stable phylogenetic markers.

## 5. Phylogenetic interpretation

| Category | Finding |
|---|---|
| Shared ancestral (branch-defining) mutations | E-Z15174 block: all three Jola men are derived where covered, which the tree predicts. |
| Your private mutations | Absent from all three men at every callable site. That is consistent with arising on your line after the E-Z15174* / E-Z5991 split. |
| Possible recurrent mutations | None detected. No private site was derived in any Jola man. |
| Sequencing / mapping artifacts | PV10: 7 low-MAPQ, heavily clipped (75 bp soft-clip), reverse-strand-only reads carry T (4 in 5309804, 3 in 5309851, 0 in 5309829). This is best explained by satellite-derived mis-placement. It is **not** evidence of sharing (see §6). |
| Missing data | 5 cells with <2 usable reads, and 6 cells with unreliable mapping. These are no-calls, not reference. |

**What this does and does not show**

* It **does** show that three Gambian Jola men belong to the sister branch of your paternal
  lineage. Their shared ancestor with you lived around 7,000 years ago, and none of them sits
  downstream of E-FTF75935.
* It **does not** show that your 18th-century paternal ancestor was Jola. A ~7,000-year TMRCA
  predates any modern ethnic group. A haplogroup found today among some Jola and Mandinka men
  (E-Z15196 also includes 1000 Genomes Mandinka HG03024/HG03027) indicates a Senegambian-region
  lineage, but it cannot assign an ethnic identity to a specific ancestor. Several statements in
  the repo `README.md` go beyond what the genetic data support: "confirmed genetic relatives",
  "4x great-grandfather Jola-Fonyi ancestor", and "the three matched families are almost certainly
  from the same community" (that last one fits 5309804/5309851, TMRCA ~375 ybp, but not 5309829).
* Of the earlier repo screen: its "3/47 carriers" figure assumed, without re-testing, that these
  three men were E-Z15174. This run confirms them at read level. The 40 other men were assessed at
  only 4 markers, often with 1–3 reads, and several were classed "ref" with a single read, so that
  frequency estimate should be treated as provisional.

## 6. Findings that warrant deeper validation

1. **PV10 (chrY:10,973,240).** Derived-allele reads appear only in the two E-FT402211 men. They
   are very likely a satellite artifact (low MAPQ, one strand, clipped fragments within a pile of
   ~100 similar reference-allele fragments), but the lineage-specific pattern is the only signal in
   the dataset that is not fully explained. To resolve it: examine your Big Y BAM at this position
   (read MAPQ and clipping), test whether YFull counts the site (it may sit in a region YFull
   excludes), and do genome-wide k-mer / full-reference alignment of the flanks. **Until then,
   exclude PV10 from phylogenetic use.**
2. **PV01–PV04 and PV08 (paralog-risk sites).** Ancestral calls here are concordant but low
   confidence. If any of these is ever cited as a branch marker (for example, when comparing with
   kit IN147690), first confirm it with long reads or high-coverage data, plus genome-wide
   mappability (UCSC Umap k100 / multi-read mappability, or alignment of tiles to the full GRCh38+decoy reference).
3. **Genome-wide mappability for all 19 sites.** The preliminary test cannot rule out autosomal
   or decoy homologs; PV09 and PV10 show that empirically. Recommended next step: align the tiles
   to the full analysis-set reference (needs the full 3.1 GB FASTA plus a ~5 GB bwa index; ask
   before downloading) or use the published Umap/Bismap k100 tracks.
4. **Branch-marker coordinates.** These are empirically supported here, but they should still be
   cross-checked against YBrowse or your Big Y VCF before being reused for screening other
   populations.

## 7. Limitations

* GGVP depth at these sites was 1–17 filtered reads (typically 3–7). Low-confidence ancestral calls
  rest on 2 reads. On the haploid Y, two concordant MAPQ-60 reads make a missed derived allele
  unlikely, but not impossible, especially at paralogous sites.
* One CRAM per man. There is no independent replicate data, and no high-coverage data for these individuals.
* Your own reads were not re-examined. The 19 positions and alleles are as you supplied them,
  with reference alleles verified against GRCh38. Please supply the Big Y VCF/BED (and ideally the
  BAM region around PV10) for the next stage.
* YFull SNP coordinates for Z-series names could not be retrieved programmatically.
