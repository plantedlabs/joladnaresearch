# Branch-classification pilot: 14 men, E-M132 → E-Z15174 → E-FTF75935 path

**Date:** 2026-11-10 · Builds on `BRANCH_DISCOVERY_REPORT.md` · **Stopped after 14 samples, as instructed.**

**Deliverables:**
* `results/branch_pilot/pilot_sample_by_marker.csv`: sample × marker matrix
* `results/branch_pilot/pilot_marker_evidence_long.csv`: every query, with read counts, strands, MAPQ, BQ and
  no-call reason
* `results/branch_pilot/pilot_placements.csv`: per-node status and placement
* `results/audit/marker_verification.csv`: pre-query verification
* Scripts: `scripts/11_verify_markers.py`, `scripts/12_branch_pilot.py`

---

## Summary

* **Both controls reproduced their known placements.** Positive: SC_GMJOL5309829 → E-Z15196. Negative/sister:
  SC_GMFUL5306400 → E-Z958, not E-CTS246, matching YFull's E-Z5985.
* **Among 12 test men:**
  * 2 GGVP Mandinka men are **E-CTS246 but not E-Z15174**: SC_GMMAN5482185 and SC_GMMAN5482207. These are
    sister-lineage men at the E-CTS246 level, not on YFull, and their CTS246 subclade is not yet determined.
  * 1 Mandinka man is E-Z958, not E-CTS246.
  * 9 men are not E-M132.
* **No man is on E-Z15174**, so no new E-Z15174 sister or intermediate lineage was found, and no man shares
  any of your Tier-A candidate SNPs.
* **No conflicting calls; 0 query errors; 6 no-calls out of 106 marker queries**, each with its reason recorded.
* Recommendation (§8): **a full screen is technically justified and cheap**, but its expected yield for
  *your* branch is low. It is worth running for the E-CTS246 population structure, preferably after adding
  sister-clade markers. Your decision.

## 1. Reconciling the 156 GGVP males

| GGVP population | Males (inferred) | FTF75935 screen (done) | Callable at FTF75935 | Classified at E-M132 → E-Z15174 before this pilot |
|---|---|---|---|---|
| Fula | 51 | yes | 37 | 0 (3 men placed by YFull; one is this pilot's negative control) |
| Mandinka | 49 | yes | 46 | 0 (none on YFull) |
| Wolof | 56 | yes | 50 | 0 (none on YFull) |
| **Non-Jola total** | **156** | yes | 133 | **0 independently genotyped** |
| Jola | 47 | yes | 39 | 3 genotyped here (13 markers); 26 placed by YFull |

**Distinction:** the FTF75935 screen asked only one question: is a man on *your terminal branch*? It cannot
detect men on E-Z15174 sister branches, E-CTS246 sister clades, or intermediate nodes, because all of those
are ancestral at FTF75935. Upstream branch classification is a separate question. Before this pilot it had
been done for no GGVP Fula, Mandinka or Wolof man. This pilot classified 10 of the 156 (5 Mandinka, 5 Wolof)
and 1 Fula control. **146 remain unclassified.**

## 2. Pre-query verification (all 24 markers PASS)

`scripts/11_verify_markers.py` checks, for each marker:
1. the GRCh38 analysis-set reference base (chrY MD5 `ce3e3110…`; the same reference is declared in GGVP and
   HGDP CRAM headers) equals REF;
2. YBrowse `snps_hg38.vcf.gz` has a record at that exact coordinate with the claimed name and identical REF/ALT;
3. YBrowse AA (ancestral allele) = REF, so derived = ALT;
4. for tree markers, the SNP name appears in the claimed node's block on the live YFull page.

| Markers | Result |
|---|---|
| 12 YFull tree markers (M132 … Z15196) | PASS on all four checks |
| FTF75935 | PASS on 1–3. It is an FTDNA-only node (YF133256 is E-Z15174* on YFull), so equivalence cannot be checked on YFull |
| 11 Tier-A private SNPs | PASS on 1–3. Treated as **candidate markers only**, never used to define a branch |

No GRCh37 data were used and no coordinate conversion was performed.

## 3. Sample selection (not a population sample)

| Role | Samples | Selection |
|---|---|---|
| Positive control | SC_GMJOL5309829 | known E-Z15196 (YFull; 34/34 concordant earlier) |
| Negative/sister control | SC_GMFUL5306400 | YFull E-Z5985 (E-Z958, not E-CTS246) |
| Test | SC_GMMAN5482228, 5482185, 5482178, 5482207, 5482240 | top-5 GGVP Mandinka males by FTF75935 depth |
| Test | SC_GMWOF5428831, 5428829, 5428855, 5428877, 5428879 | top-5 GGVP Wolof males by FTF75935 depth |
| Test | HGDP00904, HGDP00905 | first 2 HGDP Mandenka (Senegal) IDs found male by a sex-only check (HGDP00909 not used: 102 reads, female-like at 30×) |

The test men were **selected for coverage, not at random.** Nothing in this pilot estimates the frequency of
any branch in any population, and population labels are dataset sampling labels only.

## 4. Controls: observed evidence

Cells show ref/alt fragments and strand split (fwd/rev) of the supporting allele. All MAPQ medians are 60
except PV19 in the negative control (56). BQ medians are 22–33.

| Marker | Expected 5309829 | **Observed 5309829** | Expected 5306400 | **Observed 5306400** |
|---|---|---|---|---|
| M132 | D | **D** 0/5 (3/2) | D | **D** 0/7 (3/4) |
| Z958 | D | **D** 0/8 (5/3) | D | **D** 0/2 (1/1) |
| CTS246 | D | **D** 0/2 (1/1) | A | **A** 2/0 (1/1) |
| Z906 | D | **D** 0/4 (3/1) | A | **A** 6/0 (3/3) |
| Z15174 | D | N: zero coverage | A | **A** 5/0 (4/1) |
| Z15197 | D | **D** 0/3 (1/2) | A | **A** 3/0 (1/2) |
| Z15231 | D | **D** 0/3 (3/0) | A | **A** 3/0 (3/0) |
| Z15232 | D | **D** 0/5 (3/2) | A | N: single fragment (ref) |
| Z15234 | D | **D** 0/4 (1/3) | A | **A** 9/0 (3/6) |
| Z5991 | D | N: single fragment (alt) | A | **A** 9/0 (5/4) |
| FT402211 | A | **A** 5/0 (1/4) | A | **A** 7/0 (1/6) |
| Z15196 | D | **D** 0/6 (2/4) | A | **A** 8/0 (7/1) |
| FTF75935 | A | **A** 3/0 (0/3) | A | **A** 6/0 (4/2) |
| Tier-A candidates (11) | A | 10 A, 1 N (PV16 single fragment) | A | 10 A, 1 N (PV05 single fragment) |

**Every called genotype in both controls matches expectation, with 0 discordant.** The four no-calls stay
unknown. For the positive control, E-Z5991 membership is inferred from derived Z15196, its verified child; the
Z5991 site itself was not called. That is recorded as an implication, not an observation.

## 5. Sample-by-marker matrix (tree markers; Tier-A tested only where E-Z15174 was not excluded)

D = derived · A = ancestral (≥ 2 fragments, ≥ 90%) · N = no-call · . = not tested (gated out by the hierarchy)

| Sample | Label | M132 | Z958 | CTS246 | Z906 | Z15174 | Z15197 | Z15231 | Z15232 | Z15234 | Z5991 | FT402211 | Z15196 | FTF75935 | Tier-A (11) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SC_GMJOL5309829 | GGVP Jola (ctrl+) | D | D | D | D | N | D | D | D | D | N | A | D | A | 10A 1N |
| SC_GMFUL5306400 | GGVP Fula (ctrl−) | D | D | A | A | A | A | A | N | A | A | A | A | A | 10A 1N |
| SC_GMMAN5482185 | GGVP Mandinka | D | N | D | D | A | A | A | A | A | . | . | . | . | . |
| SC_GMMAN5482207 | GGVP Mandinka | D | D | D | D | A | A | A | A | A | . | . | . | . | . |
| SC_GMMAN5482240 | GGVP Mandinka | D | D | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMMAN5482228 | GGVP Mandinka | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMMAN5482178 | GGVP Mandinka | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMWOF5428831 | GGVP Wolof | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMWOF5428829 | GGVP Wolof | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMWOF5428855 | GGVP Wolof | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMWOF5428877 | GGVP Wolof | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| SC_GMWOF5428879 | GGVP Wolof | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| HGDP00904 | HGDP Mandenka (SEN) | A | A | A | A | . | . | . | . | . | . | . | . | . | . |
| HGDP00905 | HGDP Mandenka (SEN) | A | A | A | A | . | . | . | . | . | . | . | . | . | . |

Read support at the defining sites of the test men (ref/alt fragments; fwd/rev of the called allele):
* **SC_GMMAN5482185:**
  * CTS246 0/5 (3/2) and Z906 0/5 (2/3): E-CTS246 by 2 SNPs, both strands.
  * Z15174-block A at all five (2–4 fragments each).
  * M132 0/2, one strand only (2/0). Z958: single fragment (no-call).
* **SC_GMMAN5482207:**
  * M132 0/4, Z958 0/7, CTS246 0/9, Z906 0/4, all on both strands.
  * Z15174-block A at all five (3–10 fragments).
* **SC_GMMAN5482240:** M132 0/7 and Z958 0/8, then CTS246 14/0 and Z906 9/0.
* **The nine men classed as not E-M132:** M132 ancestral with 4–18 fragments, every one on both strands
  except SC_GMWOF5428829 (0/4, reverse only).

Full metrics for all 106 queries are in `pilot_marker_evidence_long.csv`.

## 6. Provisional placements

| Sample | Placement | Evidence strength |
|---|---|---|
| SC_GMJOL5309829 (ctrl+) | **E-Z15196** (E-Z5991 implied by derived child) | Strong: 9 derived tree SNPs, 0 conflicts |
| SC_GMFUL5306400 (ctrl−) | **E-Z958, not E-CTS246** | Strong: CTS246 and Z906 ancestral, 2/2 SNPs |
| SC_GMMAN5482207 | **E-CTS246, not E-Z15174** | Strong: 4 derived upstream SNPs; 5/5 block SNPs ancestral |
| SC_GMMAN5482185 | **E-CTS246, not E-Z15174** (E-Z958 implied) | Strong for E-CTS246 (2 SNPs × 5 fragments, both strands) and for not E-Z15174 (5/5 A). The M132 call itself is weak (2 fragments, 1 strand), but E-M132 is implied by E-CTS246 |
| SC_GMMAN5482240 | **E-Z958, not E-CTS246** | Strong: 2/2 CTS246-node SNPs ancestral (14 and 9 fragments) |
| 7 GGVP men + 2 HGDP men | **not E-M132** | Moderate: rests on **one defining SNP** (M132, 4–18 fragments). Ancestral at Z958/CTS246/Z906 is consistent but not independent evidence, since those lie downstream |

Updated provisional tree (pilot men in **bold**; other branches as in `results/topology/provisional_tree.md`):

```
E-M132
└── E-Z958
    ├── (not E-CTS246)  1000G GWD/MSL men; GGVP Fula 5306392/5306398/5306400 [YFull E-Z5985];
    │                   **SC_GMMAN5482240** (subclade undetermined)
    └── E-CTS246
        ├── known sister clades (E-CTS10935, E-CTS3380, E-L133): 1000G MSL, GWD, YRI, ACB, PUR men [YFull + VCF]
        ├── **SC_GMMAN5482185**, **SC_GMMAN5482207**  E-CTS246, not E-Z15174; which sister clade is UNDETERMINED
        │     (could be a known clade such as E-CTS3380 > E-Z5993 like GWD HG03048, or an unsampled one)
        └── E-Z15174
            ├── E-Z15174* / E-FTF75935: YF133256 [reported]
            └── E-Z5991 > E-FT402211 (Jola 5309804, 5309851) | E-Z15196 (Jola 5309829; GWD HG03024, HG03027)
```

## 7. Evidence assessment

| Finding | Strength | Caveat |
|---|---|---|
| Classifier and markers behave correctly | **Strong** | 2 controls: 23 called tree-marker + 20 Tier-A genotypes, 0 discordant; 24/24 markers verified |
| Two GGVP Mandinka men are E-CTS246 outside E-Z15174 | **Strong** for the node; **unresolved** for subclade | needs CTS246 sister-clade SNPs (E-CTS10935, E-CTS3380/Z5993, E-L133) |
| No pilot man is E-Z15174 / intermediate / on FTF75935 | **Strong** for these 12 men | coverage-selected; says nothing about population frequency |
| 9 men not E-M132 | **Moderate** | single defining SNP |
| HGDP at 30× works with the same reference | **Strong** | sex thresholds must scale with coverage (female HGDP00909 = 102 reads) |

**On ethnicity:** finding E-CTS246 in 2 of 5 coverage-selected Mandinka men, and E-Z15174 in none, supports no
statement about Mandinka, Wolof or Mandenka prevalence or "ownership" of any branch. It only shows that
E-CTS246 lineages outside E-Z15174 exist among GGVP Mandinka-labelled men, as they do in 1000G GWD (HG03048).

## 8. Recommendation on the full screen

**Technically justified and cheap.**
* Step 1 (4 SNPs) for the remaining 146 GGVP Fula, Mandinka and Wolof males: ~584 queries.
* Step 2 for E-CTS246-positive or unresolved men: perhaps 50–150 more.
* Total ~4.3 s/query sequential ≈ 45–55 min, or ~25 min with 2 workers. About 1–2 GB of indexed transfer,
  with GGVP indexes already cached.

**Scientific value, stated honestly:**
* **For your own branch (E-Z15174 / intermediate / FTF75935): low expected yield.** All 156 men were already
  screened at FTF75935 with zero carriers. Among the 12 pilot men, and among 55 unrelated 1000G GWD men,
  E-Z15174 appears only as HG03024/HG03027. A new E-Z15174 sister or intermediate lineage is possible but
  unlikely to be found here.
* **For the structure of E-CTS246: moderate.** The pilot already found 2 unplaced E-CTS246 men. A full screen
  would show how many GGVP men sit in E-CTS246 sister clades. That bears on the clade's regional history,
  not on your recent ancestry.

**Recommended order** (each a separate approval):
1. **Before any full screen**, add 3–5 verified sister-clade SNPs (E-CTS10935, E-CTS3380 / E-Z5993, E-L133)
   and test them on the **two pilot E-CTS246 men only**: ~10–20 queries, ~2 min. That shows whether they fall
   in known YFull clades or need a new node. It is the cheapest way to learn whether the pool holds unplaced
   E-CTS246 diversity.
2. Then, if you want the population-structure picture, run the 146-man step-1 screen with E-Z15174
   follow-up, using coverage-aware sex thresholds and 2 workers.
3. For *your* branch specifically, the higher-value evidence is still the **FTDNA E-FTF75935 Block Tree
   (IN147690)** and your **Big Y VCF**, which no public screen can replace.

Nothing beyond the 14 pilot samples has been run.
