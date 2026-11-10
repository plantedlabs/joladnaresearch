# Searching public African genomes for upstream, sister and intermediate branches of E-FTF75935

**Phase 1 (audit) and Phase 2 (design) complete. Phase 3 pilot proposed and awaiting approval.**
Date: 2026-11-10 · Builds on `REPORT.md` (v2) and `FTF75935_SCREEN_REPORT.md`

**Deliverables:**
* this report
* `results/topology/marker_matrix_wide.csv` and `marker_matrix_long.csv` (sample × marker matrix)
* `results/topology/provisional_tree.md` (provisional tree)
* `results/topology/nesting_check.md` (empirical topology test)
* `scripts/10_build_marker_matrix.py`

---

## 1. Audit of completed work

| Completed | Result | Where |
|---|---|---|
| CRAM decoding fix (MD5-verified chrY/chrX range fetch) | 2–8 s per site query | `TECHNICAL_LOG.md` |
| 19 private SNPs × 3 GGVP Jola men | 0 derived; 46/57 ancestral | `REPORT.md` §1 |
| Private-SNP reliability (Umap, ENCODE, segdup, read behaviour) | 11 tier A, 3 B, 1 C, 4 D | `results/private_snp_reliability_ranking.csv` |
| PV10 forensics | satellite artifact | `REPORT.md` §4 |
| YBrowse coordinate audit | repo "Z15234" was Z15197; true Z15234 = 15,274,693 | `results/audit/` |
| 13 tree markers × 3 Jola men | 34/34 callable genotypes match YFull | `results/tree_concordance.md` |
| FTF75935 × 394 GGVP samples | 203 males; 171 ancestral, 0 derived | `FTF75935_SCREEN_REPORT.md` |
| FTF75935 × 1000G (NYGC 3202 chrY calls) | no record; 0/93 GWD males | same |

**Not done yet** (not a repeat of earlier work): no screen at **E-CTS246 or E-Z15174 resolution** outside the
three Jola men and 1000 Genomes. The September 2026 repo screen tested 40 GGVP Jola men at four markers. One of
those markers was mislabelled, and many calls rested on a single read.

## 2. Marker inventory (GRCh38; coordinates from YBrowse; alleles verified against the MD5-checked reference)

| Node | SNP(s) tested | GRCh38 | Anc>Der | Mappability (Umap k100 / blacklist) |
|---|---|---|---|---|
| E-M132 | M132 | 19,734,375 | G>T | 1.0 / – |
| E-Z958 | Z958 (= CTS10713) | 20,558,142 | C>T | 1.0 / – |
| E-CTS246 | CTS246 (= V1119) | 2,972,260 | G>A | 1.0 / – |
| E-CTS246 | Z906 | 8,360,264 | T>C | 1.0 / – |
| E-Z15174 (FTDNA: E-Z15234) | Z15174 | 4,202,358 | G>T | 1.0 / – (X-transposed region; caution) |
| E-Z15174 | Z15197 | 8,414,659 | C>T | 1.0 / – |
| E-Z15174 | Z15231 | 14,719,978 | G>A | 1.0 / – |
| E-Z15174 | Z15232 | 14,765,245 | T>C | 1.0 / – |
| E-Z15174 | Z15234 | 15,274,693 | C>G | 1.0 / – |
| E-Z5991 | Z5991 | 2,849,335 | C>T | 1.0 / – |
| E-FT402211 | FT402211 | 9,176,816 | A>G | 1.0 / – (segdup 0.965; caution) |
| E-Z15196 | Z15196 | 8,336,326 | A>C | 1.0 / – |
| E-FTF75935 | FTF75935 | 2,895,562 | A>C | 1.0 / – |
| YF133256 private (tier A) | PV05 PV06 PV07 PV11 PV12 PV13 PV14 PV15 PV16 PV17 PV19 | `config/tierA_private_variants.tsv` | | 1.0 / – |

**Build conversion:** none needed. All data used or proposed (GGVP, 1000G NYGC, HGDP) are GRCh38 against the
same analysis-set reference. For HGDP the chrY header M5 is `ce3e3110…`, identical to GGVP. If GRCh37 data are
ever added, use the YBrowse hg19 records together with an MD5-verified hg19 chrY, never coordinate arithmetic.

## 3. Topology verified from the evidence, not assumed

`results/topology/nesting_check.md` tests every parent→child pair on all observed genotypes: 3 GGVP Jola men from
CRAM, and 38 1000G E-M132 men from VCF.

| Relationship | Observed | Verdict |
|---|---|---|
| Z958 within M132 | 41/41 Z958 carriers are M132-derived | nested |
| CTS246 within Z958 | 22/22 | nested |
| CTS246 ≡ Z906 | identical carrier sets | same node (equivalent SNPs) |
| Z15174 block within CTS246 | 4/4 (Z15174); 0 conflicts | nested |
| Z15174, Z15197, Z15231, Z15232, Z15234 | identical carriers in 1000G (3); no conflicts in Jola | same node (block) |
| Z5991 within Z15174 block | 5/5 | nested |
| FT402211 and Z15196 | **mutually exclusive** (0 men carry both; 2 vs 4) | **sister** branches under Z5991 |
| FTF75935 | 0 carriers among observed public men | position below E-Z15174 only from FTDNA/YFull |

**Topology is not linear.** E-Z15174 splits into E-Z5991 (five public men) and your own line (E-Z15174* on
YFull, E-FTF75935 on FTDNA). E-Z5991 splits again into E-FT402211 and E-Z15196. Your ancestral state at Z5991
comes from YFull's placement. It is **not observed** here, because your VCF has not been supplied.

In 1000G alone, Z5991 and Z15196 cannot be separated, because the same 3 men carry both. Separating them depends
on the two GGVP Jola men (Z5991+ Z15196−) and on YFull's placement of you.

## 4. Upstream and sister clades already present in public data

The E-CTS246 sister clades of E-Z15174 are already placed by YFull and confirmed here as **CTS246+ but
ancestral at all five E-Z15174 block SNPs** (NYGC VCF):

| Sister clade (YFull) | Public men | Sampling label |
|---|---|---|
| E-CTS10935 | HG03096 (+ son) | 1000G MSL (Mende, Sierra Leone) |
| E-CTS3380 > E-Z5993 | HG03048 (+ son) | 1000G GWD (Gambian, Western Division) |
| E-CTS3380 > E-FT296631 | HG01094 (+ son) | 1000G PUR |
| E-CTS3380 > E-P110 | NA19239, NA19189 (+ son) | 1000G YRI |
| E-L133 (> E-CTS736 / E-Z5992) | NA18910, NA18856 (+ sons); HG02449, HG01882 (+ son) | 1000G YRI; ACB |

The non-CTS246 E-Z958 branches (mostly E-Z5985) contain 18 1000G men (GWD, MSL), 3 GGVP Fula men and ≥ 20 GGVP
Jola men. Those placements are reported by YFull.

**What this means:** every E-CTS246 man in 1000 Genomes is already classified. **No new E-Z15174 sister or
intermediate lineage can come from 1000G.** The only E-Z15174 men there are HG03024, HG03027 and HG03029.

## 5. Where new branches could still be found (gap analysis)

YFull's E-M132 tree lists 26 GGVP Jola, 3 GGVP Fula, **0 GGVP Mandinka, 0 GGVP Wolof**, 20 1000G men and
1 HGDP man (HGDP01284, Senegal, 'mnk', on a non-Z958 E-M132 branch). 1000G GWD Mandinka have 11 E-M132 men
among 55 unrelated males (20%), so the absence of GGVP Mandinka and Wolof most plausibly means **YFull has not
processed them**. It is not evidence that these lineages are absent.

| Dataset | Males | E-CTS246 / E-Z15174 status | Priority |
|---|---|---|---|
| GGVP Mandinka | 49 | **untested**; not on YFull | **1** |
| GGVP Wolof | 56 | **untested**; not on YFull | **1** |
| GGVP Fula | 51 | untested (3 men on YFull, in E-Z5985) | **1** |
| HGDP Mandenka (Senegal), EBI high coverage | ≤ 20 (sex not yet inferred) | untested; reference identical | **2** |
| SGDP West Africa (EBI) | few | untested | 3 |
| GGVP Jola | 47 | 26 on YFull + 40 in the Sept 2026 screen (weak markers) | 4 (re-check only) |
| 1000G (all) | ~1,233 | fully classified (§4) | done |

## 6. Multi-marker classification strategy

Hierarchical design, so each man gets only the queries his result justifies:

| Step | Markers | Who | Decision rule |
|---|---|---|---|
| 1 | CTS246, Z906 (equivalent) | every male | E-CTS246 if ≥ 1 marker D with ≥ 2 fragments and no A at the other; NOT CTS246 if both A (or 1 A plus 1 no-call); otherwise **unresolved** |
| 2 | M132, Z958, plus the 5 Z15174-block SNPs | Step-1 derived or unresolved | **E-Z15174 if ≥ 2 block SNPs D and 0 A**. Block D and A mixed → **conflicting**, so stop and inspect reads. 0 D and ≥ 2 A → **CTS246 sister lineage** |
| 3 | Z5991, FT402211, Z15196, FTF75935 + 11 tier-A private SNPs | E-Z15174 men | assign E-Z5991 / E-FT402211 / E-Z15196 / E-FTF75935, or **E-Z15174 x (new)** if Z5991 A and FTF75935 A |

Call vocabulary for each man × marker:
* D derived · A ancestral (≥ 2 fragments, ≥ 90% agreement)
* N missing (0–1 fragments) · M ambiguous (mixed)
* X conflicting (call contradicts the established nesting)
* nr no record (VCF sources only; not proof of ancestral)

A man is assigned to a node only through concordant evidence at **≥ 2 independent SNPs** of that node or below.

**What would count as an "intermediate" branch with you:** a man who is E-Z15174 derived, Z5991 ancestral and
FTF75935 ancestral, **and** shares ≥ 2 tier-A private SNPs of YF133256. That pattern would place him on a branch
between E-Z15174 and E-FTF75935, a genuinely new node. If he is E-Z15174 derived, Z5991 ancestral and shares
none of your SNPs, he is a **new E-Z15174 sister lineage**, which is still informative about the clade's age
and distribution. Neither would be claimed from one SNP or one low-quality read. Every candidate gets
read-level review, and ideally a high-coverage re-test.

## 7. Proposed pilot (smallest informative)

| | Samples | Queries | Purpose |
|---|---|---|---|
| Positive control | SC_GMJOL5309829 (E-Z15196) | 2 + 7 + 15 | must run through to E-Z15196 |
| Negative/sister control | SC_GMFUL5306400 (YFull E-Z5985: Z958+, CTS246−) | 2 (+7 if needed) | must stop at "not CTS246" |
| Unscreened | 5 GGVP Mandinka + 5 GGVP Wolof males with the highest FTF75935 depth (proxy for chrY coverage) | 20, plus follow-ups for any CTS246+ | tests the classifier on the untested gap |
| High-coverage | 2 HGDP Mandenka (sex check first) | 2–24 | tests 30x CRAM cost and gives definitive calls |

**Estimate:**
* 14 samples
* ~30–60 site queries (most men stop after step 1)
* **~5–10 min**
* **~100–300 MB** of indexed byte ranges (GGVP indexes are already cached; 2 new HGDP indexes)

No full CRAMs and no population-scale screen.

**What would follow, after your approval:** step 1 on all 156 untested GGVP Fula, Mandinka and Wolof males
(~312 queries, ~25 min, ~0.5–1 GB), then steps 2–3 for positives only. In GWD, 3 of 55 unrelated men are
CTS246+ (~5%). That suggests roughly 5–10 CTS246+ men and perhaps 1–4 E-Z15174 men among the 156, if GGVP
resembles GWD. That expectation is uncertain, not a prediction.

## 8. Interpretation so far

**Where E-Z15174 has been observed** (public, independent lineages):
* 3 GGVP Jola men, two tagged 'dyo' by YFull
* 2 1000G GWD Mandinka men ('mnk'), plus one son
* 1 African American (you)

All public carriers were sampled in western Gambia. The E-Z5991 subclade is shared across Jola- and
Mandinka-labelled men within its ~2,300-year TMRCA: E-Z15196 contains both, with a TMRCA of ~1,000 years.

Distinguishing the three concepts:
* **Sample ethnicity:** what the dataset records, i.e. GGVP sampling population or YFull language tag. Even
  this varies: two GGVP "Jola" men in E-Z5985 carry YFull's 'ful' tag.
* **Geographic origin:** all public carriers come from Western Gambia (GGVP and GWD sampling regions).
* **Inferred ancestral population of the lineage:** **unresolved.** A clade with a TMRCA of about 7,000 years,
  known from six independent men, cannot be assigned to any present-day ethnic group.

**Assessment against your four options:**

| Hypothesis | Current evidence |
|---|---|
| Jola-associated | Partially consistent: 3 of the 5 public African carriers are GGVP Jola. But Jola are also the most heavily sampled and YFull-processed GGVP group (26 placed), so this is confounded by sampling |
| Fulani-associated | **Not supported by observed data**: 0/37 callable GGVP Fula at FTF75935; Fula untested at E-Z15174; IN147690's placement unverified |
| Broader Senegambian | **Best supported**: carriers in two ethno-linguistic groups (Jola, Mandinka), all from western Gambia |
| Unresolved | Appropriate overall: n = 6 independent lineages, uneven sampling, Mandinka/Wolof/Fula GGVP not yet tested at E-Z15174 |

No ethnic probabilities are given. With carriers this rare and sampling this uneven across groups (Jola heavily
processed by YFull; Mandinka and Wolof GGVP not at all), no defensible sampling model exists yet. The pilot and
step-1 screen would fix the most important part of that imbalance.

## 9. Evidence-strength assessment

| Claim | Strength | Basis |
|---|---|---|
| The three Jola men are E-Z15174 > E-Z5991 (> FT402211 / Z15196) | **Strong** | 34/34 concordant read-level genotypes; matches YFull |
| E-Z15174 block SNPs are one node; FT402211 and Z15196 are sisters | **Strong** | identical carrier sets; mutual exclusivity in observed data |
| 1000G holds no further E-Z15174 or new E-CTS246 lineages | **Strong** | 30x joint calls; every CTS246 carrier matches a YFull placement |
| No public Gambian man carries FTF75935 | **Moderate–strong** for sampled men; **weak** for populations | 0/171 callable GGVP; no 1000G record; 95% UB 1.7% |
| You are ancestral at Z5991 | **Moderate** (reported) | YFull placement; not observed in your VCF here |
| IN147690 is on E-FTF75935 | **Unverified** | no data |
| E-Z15174 is "Jola" or "Fulani" | **Not supported**; unresolved | see §8 |

## 10. Three most informative next analyses

1. **Pilot, then step-1 CTS246/Z906 screen of the 156 GGVP Fula, Mandinka and Wolof males** (§7). This is the
   only public data where new E-Z15174 sister or intermediate lineages can still be found at scale, and it
   corrects the Jola-heavy sampling of current carriers.
2. **Your FTDNA E-FTF75935 Block Tree, plus your Big Y VCF.** These give IN147690's branch position and the
   private-SNP counts, confirm your Z5991 ancestral state from your own data, and supply the E-FTF75935 block
   SNPs for multi-SNP testing of any candidate.
3. **HGDP Mandenka (Senegal) at high coverage.** These are 20 30x genomes with an identical reference, from a
   different country than all current carriers. Any E-Z15174 man there would extend the known distribution
   beyond western Gambia, with definitive genotypes.
