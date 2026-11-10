# Y-chromosome comparison of YF133256 with three Gambian Jola men: research report (v2)

**Subject:** YF133256 (FTDNA terminal E-FTF75935; YFull E-Z15174*), direct paternal line of an enslaved
ancestor, Woodley Drumgoole (b. ~1809–1820, Virginia)
**Comparison samples:** Gambian Genome Variation Project (GGVP) low-coverage GRCh38 CRAMs: SC_GMJOL5309804,
SC_GMJOL5309851, SC_GMJOL5309829
**Tree reference:** YFull YTree v14.07.00 (02 Nov 2026). Marker coordinates from YBrowse `snps_hg38.vcf.gz` (09 Nov 2026)
**Date:** 2026-11-10 · v1 of this report (same day) is superseded

---

## Summary

1. **None of the three Jola men shares any of your 19 private variants.** Across 57 sample × site cells,
   all 284 strictly filtered reads show the ancestral (reference) allele, and none shows your derived allele.
2. **Their paternal lineages sit exactly where YFull places them.** At 13 independently named tree
   markers from E-M132 down to each man's terminal branch, **34 of 34 callable genotypes agree with the
   published tree, 0 disagree**. That includes the five branches where a man must be *ancestral*.
3. **Your closest relationship to these men is through E-Z15174**, a shared ancestor that YFull estimates
   at ~7,000 years ago (95% CI 8,400–5,700). Nothing in the data suggests a closer relationship.
4. **The only derived-allele signal (PV10) has been traced to satellite DNA from elsewhere in the genome.**
   It is not a shared mutation.
5. **About 11 of your 19 private SNPs are robust phylogenetic markers.** Four are unsuitable (PV01, PV02,
   PV04, PV10), and four need caution.
6. The earlier repo scripts **mislabelled one marker**: chrY:8,414,659 is Z15197, not Z15234.
   The true Z15234 is at chrY:15,274,693.
7. The most informative next comparison is **SNP-level data from FTDNA kit IN147690**, plus a **public-data
   screen for FTF75935** (your terminal SNP) in Gambian men (§7).

---

## 1. Established findings

These are supported by multiple independent lines of evidence and are unlikely to change with more data
of the same kind.

### 1.1 The three Jola men do not carry your private variants

Strict filters: primary, non-duplicate, MAPQ ≥ 20, base quality ≥ 20, ≥ 50% of the read aligned, and no
better alignment elsewhere (AS ≥ XS). Haploid calling rules: ≥ 2 reads and ≥ 90% agreement are needed for
any call. A site with zero or one read is a no-call, never "reference".

| | Derived | Ancestral | Mixed | Too few reads | Unreliable mapping |
|---|---|---|---|---|---|
| SC_GMJOL5309804 (E-FT402211) | 0 | 17 | 0 | 0 | 2 |
| SC_GMJOL5309851 (E-FT402211) | 0 | 14 | 0 | 3 | 2 |
| SC_GMJOL5309829 (E-Z15196) | 0 | 15 | 0 | 2 | 2 |

* Ancestral calls: 23 high, 10 moderate and 13 low confidence. The full per-cell matrix is in
  `results/evidence_matrix.md` and `results/jola_site_evidence_strict.csv`.
* At the 11 Tier-A sites (§3), all three men are ancestral wherever covered (31/33 cells callable).
  The strict filter changed no call compared with the standard filter.

### 1.2 The three men's positions in the tree are confirmed at read level

| Marker (YBrowse) | GRCh38 | Branch | 5309804 | 5309851 | 5309829 |
|---|---|---|---|---|---|
| M132 | 19,734,375 G>T | E1a | D ✔ (12) | D ✔ (3) | D ✔ (5) |
| Z958 | 20,558,142 C>T | E-Z958 | D ✔ (6) | D ✔ (2) | D ✔ (8) |
| CTS246 | 2,972,260 G>A | E-CTS246 | D ✔ (6) | no-call (1 read) | D ✔ (2) |
| Z906 | 8,360,264 T>C | E-CTS246 block | D ✔ (8) | D ✔ (5) | D ✔ (4) |
| Z15174 | 4,202,358 G>T | E-Z15174 | D ✔ (5) | no-call (0) | no-call (0) |
| Z15197 | 8,414,659 C>T | E-Z15174 | D ✔ (4) | no-call (0) | D ✔ (3) |
| Z15231 | 14,719,978 G>A | E-Z15174 | D ✔ (6) | D ✔ (7) | D ✔ (3) |
| Z15232 | 14,765,245 T>C | E-Z15174 | D ✔ (9) | D ✔ (3) | D ✔ (5) |
| Z15234 | 15,274,693 C>G | E-Z15174 | D ✔ (12) | D ✔ (4) | D ✔ (4) |
| **Z5991** | 2,849,335 C>T | Jola sister clade | D ✔ (8) | D ✔ (4) | no-call (1) |
| **FT402211** | 9,176,816 A>G | E-FT402211 | D ✔ (10) | D ✔ (10) | **A ✔ (5)** |
| **Z15196** | 8,336,326 A>C | E-Z15196 | **A ✔ (8)** | **A ✔ (5)** | D ✔ (7) |
| **FTF75935** | 2,895,562 A>C | **your terminal** | **A ✔ (3)** | **A ✔ (5)** | **A ✔ (3)** |

D = derived, A = ancestral, (n) = supporting reads. Every supporting read agreed: 167 derived and 29
ancestral reads at the positions where those alleles were expected, with no contrary read anywhere.

The ancestral calls at FT402211, Z15196 and FTF75935 matter most. They show the pipeline is not simply
calling "derived" everywhere. They also place all three men off your E-FTF75935 branch, using a marker
that is well-mapped genome-wide (Umap k100 = 1.0).

### 1.3 The CRAM decoding problem is solved

The original timeouts were caused by htslib fetching the 57 MB chrY reference from the ENA server at
~85 KB/s. The pipeline now range-fetches chrY directly from the exact 1000 Genomes reference file and
checks its MD5 (`ce3e3110…`) against the CRAM header. Each site query takes 2–8 s. See `TECHNICAL_LOG.md`.

### 1.4 Coordinate audit of the positive-control markers

The original four "E-Z15174" markers came from the repo's `run_all.sh`. Two separate questions apply to them:

| Question | Z15174 (4,202,358) | "Z15234" (8,414,659) | Z15231 (14,719,978) | Z15232 (14,765,245) |
|---|---|---|---|---|
| **A. Empirical:** are the Jola men derived at this position? | yes | yes | yes | yes |
| **B. Independent:** does a public registry assign the claimed name to this coordinate? | **yes** (YBrowse) | **no: YBrowse says Z15197** | **yes** | **yes** |
| Is the named SNP in YFull's E-Z15174 block? | yes | Z15197: yes | yes | yes |

Conclusion: question A alone could never confirm a SNP's *name*. A derived allele at a coordinate shows
only that the men carry *some* variant there that you also carry. The YBrowse audit (B) corrects one
label. Because Z15197 is also an E-Z15174-block SNP, the earlier positive-control conclusion still holds.
The true Z15234 (15,274,693 C>G) has now been tested directly, and all three men are derived. The repo's
`run_all.sh`, `check_jola_males.py` and `README.md` still carry the wrong label.

Remaining caveat: YBrowse is a curated aggregator, not a primary source. The strongest confirmation would
be your own Big Y VCF, which lists FTDNA's SNP names at these positions.

---

## 2. Provisional observations

These are supported by current evidence but could change with new data.

* **PV13 (15,391,518) may be a mutable site.** YBrowse also lists FT365779 at the same position with a
  *different* derived allele (G>C; yours is G>A). That suggests the position has mutated independently
  in another lineage. It does not mean anyone else shares your mutation.
* **Names on your private SNPs.** 17 of 19 positions already have YBrowse names (FTG…, FTC76522,
  Y621321–Y621323), and PV09/PV10 have none. These are probably names assigned to your own novel variants
  by FTDNA and YFull. A name alone does **not** show that anyone else carries the variant. This can only be
  checked against the FTDNA Block Tree or YFull's SNP pages for those names.
* **PV10's satellite signal appears only in the E-FT402211 pair** (§4). The most economical explanation is
  a heritable satellite-array sequence variant on that lineage. It is untested and does not bear on your SNPs.
* **Frequency estimate from the earlier repo screen.** "3/47 Jola males carry E-Z15174" depends on 4 markers
  (one mislabelled), many with 1–3 reads, and some men classed as "ref" on one read. Treat it as provisional.
  A re-screen using the better-mapped markers in §1.2 would be quick (§7).

---

## 3. Ranking of the 19 private SNPs by phylogenetic reliability

Lines of evidence:
* Big Y MAPQ, which you supplied
* **Published genome-wide mappability**: Umap k24/36/50/100 (Karimzadeh et al. 2018, NAR 46:e120), queried
  remotely from UCSC
* ENCODE blacklist v2 and GRC exclusions
* UCSC segmental duplications and RepeatMasker
* The preliminary chrX/chrY 100-mer test
* Jola read behaviour: fraction of low-MAPQ reads, and whether the strict filter changes calls
* Read-level forensics, used for PV10 only

Scoring is additive and auditable: `scripts/05_rank_private_snps.py`, output
`results/private_snp_reliability_ranking.csv`.

| Rank | SNP | GRCh38 | Tier | Main evidence |
|---|---|---|---|---|
| 1 | PV05 | 6,926,993 A>G | **A** | Umap k100 1.0; no flags; ANC in 3/3 |
| 2 | PV06 | 7,013,476 T>C | **A** | Umap k24–k100 all 1.0; ANC in 3/3 (17/7/7 reads) |
| 3 | PV07 | 7,758,704 C>A | **A** | Umap all 1.0; ANC 3/3 |
| 4 | PV11 | 14,952,809 C>A | **A** | Umap k100 1.0; ANC 3/3 |
| 5 | PV14 | 15,936,416 G>A | **A** | Umap all 1.0; ANC 3/3 |
| 6 | PV15 | 15,990,294 C>A | **A** | Umap k100 1.0; ANC 3/3 |
| 7 | PV17 | 21,349,919 G>A | **A** | Umap all 1.0; ANC 3/3 |
| 8 | PV13 | 15,391,518 G>A | **A** | Umap k100 1.0; ANC 3/3; other allele named at site (possible mutable site) |
| 9 | PV12 | 15,252,463 T>C | **A** | Umap all 1.0; one Jola man has a single read |
| 10 | PV16 | 18,949,120 C>T | **A** | Umap all 1.0; one Jola man has a single read |
| 11 | PV19 | 26,572,350 C>T | **A** | Umap k100 1.0; inside a 96.2% segdup; ANC 3/3 |
| 12 | PV08 | 9,023,319 C>A | B | chrY paralog (~21.8 Mb) near-identical over 100 bp; Umap k50 0.8 |
| 13 | PV18 | 21,680,474 C>G | B | segdup 96.5%; Umap k36 0.11, k50 0.56 (unique only at ≥100 bp) |
| 14 | PV09 | 10,820,475 A>T | B | **ENCODE blacklist** (pericentromeric high-signal region); 46% of Jola reads MAPQ<20 |
| 15 | PV03 | 5,650,649 G>T | C | X-transposed region (segdup 98.8% with Xq21); Umap k50 min 0.44 |
| 16 | PV04 | 6,220,233 T>C | D | X-transposed region; **Big Y MAPQ<30**; Umap k50 min 0.18 |
| 17 | PV02 | 5,039,793 C>T | D | X-transposed region; **Big Y MAPQ<30**; **Umap k100 = 0.74** |
| 18 | PV10 | 10,973,240 C>T | D | **ENCODE blacklist**; satellite edge; 88% low-MAPQ reads; forensics (§4) |
| 19 | PV01 | 3,339,724 G>C | D | **Big Y MAPQ<30**; **Umap k100 = 0.48**; X/Y 100-mer test fails |

**Interpretation:** Tier A (11 SNPs) is suitable for phylogenetic use and age estimation. Tier B (3) is
usable, but each call should be checked against raw reads. Tiers C and D (5) should not be relied on as
branch markers without long-read or high-coverage confirmation. Your Big Y's three low-MAPQ sites all land
in tier D, and the published Umap tracks independently confirm PV01 and PV02 are not uniquely mappable at 100 bp.

---

## 4. PV10 (chrY:10,973,240 C>T): forensic read analysis

**Question:** seven reads carry your derived T. They come only from the two E-FT402211 men: 4 in
SC_GMJOL5309804 and 3 in SC_GMJOL5309851, and none in SC_GMJOL5309829. Could they come from another
genomic location?

**Answer: yes. Every line of evidence points to another location.** The full per-read table is
`results/pv10/pv10_read_table.tsv`, and the raw SAM records are `results/pv10/*.sam`.

| Feature | 7 T reads | 24 C reads with MAPQ ≥ 20 |
|---|---|---|
| Aligned length | **25 bp** (44–75 bp soft-clipped) | mostly 100 bp (full length) |
| MAPQ | 0, 0, 1, 1, 6, 6, 12 | mostly 60 |
| bwa found a *better* alignment elsewhere (XS > AS) | **7/7** (AS 20 vs XS 30–35) | 3/24 |
| Alternative hits listed (XA) | **49–72 per read**; 30–49 exact (NM=0) | 0 for the full-length reads |
| Where the alternatives are | 122 primary autosome/X hits (96 exact), 129 chrY, 101 unlocalized `*_random`, 65 decoy, 28 `chrUn`, 14 alt contigs | – |
| Strand and start | all **reverse** strand, **all start at 10,973,222** | both strands, varied starts |
| Mate MAPQ | 1–19 (mates also clipped) | mostly 60 |
| Sequence | shifted copies of one satellite-like source sequence | unique flanking sequence |

Two direct tests:
1. The 25-mer the T reads align by, `TTCCATTCCATTCCATTCTATTCCA`, occurs **exactly 47 times elsewhere on
   reference chrX/chrY**. The reference-allele version occurs 731 times. A read carrying the T is therefore
   explained perfectly by dozens of other loci.
2. The rest of these reads contains a sequence (31-mer probe) with **zero exact matches on reference chrX
   or chrY**. It is found in exactly the 7 T reads and in none of 5309829's 89 reads. So the reads
   originate from satellite DNA that is either unassembled or lies on another chromosome. bwa placed a
   25 bp fragment here because the mate rescued the pairing.

**Why the published track misses this:** Umap scores PV10 as uniquely mappable (k100 = 1.0) because the
*reference* 100-mer there is unique. The problem is sequence present in the men's genomes but missing
from, or collapsed in, the reference. ENCODE blacklists the whole region (chrY:10.25–11.04 Mb) for exactly
this behaviour.

**Conclusion:** the PV10 alternate reads are mapping artifacts from satellite DNA. They are **not**
evidence that the E-FT402211 men share your PV10 variant. PV10 should be excluded from phylogenetic use.
Your own Big Y call at PV10 deserves the same scrutiny, using read-level MAPQ and clipping in your BAM.

---

## 5. Unresolved questions

1. **Is PV10 real in *your* genome?** PV10 is in a blacklisted satellite region, so your own derived call
   may be the same kind of artifact. Resolving it needs your Big Y BAM at chrY:10,973,100–10,973,400.
2. **Are your 19 "private" variants private relative to IN147690?** This depends on whether FTDNA computed
   them before or after he was placed on E-FTF75935 (§6).
3. **What is E-Z15174's distribution beyond these samples?** No Fulani, Wolof or Mandinka men have been
   tested yet, apart from the two 1000 Genomes Mandinka men already on YFull's E-Z15196.
4. **Are the tier B–D sites reliable genome-wide?** Published tracks are now used, and a full-genome
   re-alignment would add little for tiers A and D. It would only help resolve PV08, PV18, PV09 and PV03.

---

## 6. Gambian Fulani tester IN147690: what evidence would establish a recent shared ancestor?

**Current status:** association with your branch only. The SNP-level relationship is unverified and his
data are not public. YFull still shows you alone on E-Z15174*, so his data are not in YFull.

**Structural point.** FTDNA creates a named branch such as E-FTF75935 only when at least two testers share
its defining variants. Your terminal SNP FTF75935 is **not** among your 19 private variants. That implies
FTF75935, and any equivalent SNPs in the same block, are shared with at least one other tester, plausibly
IN147690. This is an inference, not a confirmed fact.

**Evidence needed, from most to least decisive:**

| # | Evidence | What it would show | How to obtain it |
|---|---|---|---|
| 1 | FTDNA Block Tree for E-FTF75935: number of SNPs in the branch block, and the **private-variant count for each tester** on it | Whether he is on your branch, and how deep your shared ancestor is | Visible to you on FTDNA (Block Tree / Discover), or via the E-M132 or Gambia project administrators |
| 2 | His position relative to your 19 private SNPs: derived or ancestral at each | If he is derived at some, they are not private but shared and recent, which is the strongest signal of a recent relationship. If he is ancestral at all, your shared ancestor predates all 19 mutations | FTDNA block tree (shared variants move into the branch), or his VCF/BAM with consent |
| 3 | His Big Y BAM or VCF, analysed with this pipeline restricted to **tier A sites** and matched callable regions | Independent, artifact-controlled confirmation | His written consent; an upload to YFull would also place both of you on one YFull branch |
| 4 | Y-111 STR genetic distance | Weak supporting evidence only; STRs saturate over long time spans | FTDNA match list |

**Calibration with your own data.** Your line carries 19 private SNPs below E-FTF75935, 11 of them tier A.
SNP-clock calibrations depend on the callable region and are uncertain; FTDNA and YFull use different
scales. As a rough guide, though, **11–19 accumulated SNPs on one line corresponds to roughly a thousand
years or more, not a few centuries**. If IN147690 shares only the FTF75935 block and none of your 19
private SNPs, your common ancestor with him most likely lived **well before the transatlantic slave trade
era**. A genealogically recent connection, within the last ~300 years, would require him to share most of
your private SNPs. The block-tree counts in row 1 would settle this immediately. No public dataset can
substitute for them.

---

## 7. Most informative next comparisons (no full-genome download needed)

1. **IN147690 block-tree counts (§6, row 1).** Zero computational cost, and the highest information value.
2. **Public screen at FTF75935 (chrY:2,895,562 A>C)** across all GGVP males (Jola, Fula, Mandinka, Wolof)
   and 1000 Genomes GWD (Gambian Mandinka, high coverage). FTF75935 is well-mapped (Umap k100 1.0, no
   blacklist or segdup), and all three Jola men are confirmed ancestral. **Any man derived at FTF75935
   would be on your own branch**: the first public-data relative closer than 7,000 years.
   Cost: ~1 query per man, ~5 s each, ≈ 35–45 min for ~400 men, a few hundred MB transferred.
   Run with a sex/coverage check first, and confirm any positive at Z15174-block markers and tier-A private SNPs.
3. **Your Big Y VCF/BED (and the BAM region around PV10).** This confirms the marker names at the source,
   lets me add FTF75935's block equivalents as extra screening markers, and settles question 5.1.
4. **Negative controls.** Run the 19 private sites on 2–3 non-E-Z15174 Jola men (`ROLES=control`).
   Expected: ancestral throughout. This guards against reference-specific artifacts.

---

## 8. What these findings do and do not mean

* They **do** show that your paternal lineage and these three Jola men's lineages share an ancestor on
  E-Z15174, around 7,000 years ago, and that no tested Jola man is on your E-FTF75935 branch.
* They **do not** identify the ethnicity of your 18th-century African ancestor. A lineage that is ~7,000
  years old predates modern ethnic groups. E-Z15174 is currently known from Jola and Mandinka men in
  western Gambia (YFull), which points to Senegambian regional ancestry, but it does not assign a people
  or community to a particular ancestor. Several statements in the repo's top-level `README.md` go beyond
  the genetic evidence: "confirmed genetic relatives", "4x great-grandfather Jola-Fonyi ancestor", and
  "almost certainly from the same community".
* The finding closest to a "recent" relationship is still open, and depends on IN147690 (§6) and on the
  FTF75935 screen (§7.2).

---

## 9. Limitations

* GGVP coverage at these sites was 1–17 filtered reads (typically 3–7). Low-confidence calls rest on 2 reads.
* There is one CRAM per man and no replicate or high-coverage data for these individuals.
* Your own reads were not examined. Positions and alleles are as you supplied them, with reference alleles
  verified against GRCh38.
* Published mappability (Umap) describes the reference genome, not each person's genome. PV10 shows that
  real reads can behave worse than the track predicts.
* The reliability score uses explicit but subjective weights. The CSV lists each component so it can be re-weighted.

## Files

| File | Content |
|---|---|
| `results/jola_site_evidence_strict.csv` | primary evidence table (strict filter; 19 private + 13 tree markers × 3 men) |
| `results/jola_site_evidence.csv` | standard-filter run (v1) |
| `results/evidence_matrix.md` | 57-cell private-variant matrix |
| `results/tree_concordance.md` | expected vs observed at the 13 YBrowse-named markers |
| `results/private_snp_reliability_ranking.csv` | ranking with all evidence components |
| `results/genomewide_mappability.csv` | Umap k24–k100, ENCODE blacklist, GRC exclusions |
| `results/site_annotation.csv` | reference check, repeats, segdups, preliminary chrX/chrY test |
| `results/pv10/` | PV10 SAM records and per-read forensic table |
| `results/audit/` | YBrowse position→name lookups |
