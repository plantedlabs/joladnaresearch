# Upstream branch screen of all GGVP Fula, Mandinka and Wolof males (E-M132 → E-Z15174)

**Date:** 2026-11-10 · Builds on `BRANCH_PILOT_REPORT.md` · Approved scope only: no new datasets or analyses

**Deliverables:**
* `results/upstream_screen/ggvp_nonjola_upstream_placements.csv`: all 156 men, one row each, with evidence
  strings and the basis for excluding E-Z15174
* `results/upstream_screen/ggvp_nonjola_sample_by_marker.csv`: sample × marker matrix
* `results/upstream_screen/pilot_marker_evidence_long.csv`: all 709 screen queries (fragments, ref/alt/other,
  strands, MAPQ, BQ, confidence, no-call reason)
* `results/upstream_screen/upstream_summary.md`: count tables and caveated proportions
* `results/topology/provisional_tree.md`: updated tree
* Scripts: `scripts/12_branch_pilot.py` (classifier), `scripts/14_summarize_upstream.py`

---

## Answer

**No new E-Z15174 carriers.** None of the 156 GGVP Fula, Mandinka and Wolof males carries a derived allele at
any of the five E-Z15174 block SNPs. Not a single derived read was seen at those sites in any screened man,
apart from one low-quality read that the filters removed (SC_GMFUL5306420, Z15197).

| Status at E-Z15174 | Men | Basis |
|---|---|---|
| Excluded | 23 | ≥ 2 of the 5 block SNPs called ancestral (2–10 fragments each) |
| Excluded | 132 | parent node E-CTS246 ancestral at **both** CTS246 and Z906 |
| **Unknown** | 1 | SC_GMWOF5428755, the lowest-coverage male; 0–1 fragment at nearly every site, kept as unknown |
| Derived | **0** | |

Because no man was E-Z15174, step 3 (Z5991, FT402211, Z15196, FTF75935 and the 11 Tier-A candidates) was not
triggered for anyone except the one unknown man, whose step-3 sites were also almost all no-calls.

**Relationship to your branch:** no GGVP Fula, Mandinka or Wolof man shares your E-Z15174 branch. That means
no screened man is closer to you on the paternal line than the ~11,000-year-old E-CTS246 split (YFull
formation estimate for E-Z15174; TMRCA of E-CTS246 ~11,000 ybp). The closest men found are **7 E-CTS246
sister-lineage men** (below). They share a paternal ancestor with you at or above E-CTS246, deep in prehistory,
not in a genealogical time frame.

## 1. Audit and scope

| Item | Value |
|---|---|
| GGVP non-Jola males (inferred from sex window) | 156 (Fula 51, Mandinka 49, Wolof 56) |
| Classified in the pilot (coverage-selected test men) | 10 (5 Mandinka, 5 Wolof): reused, not re-queried |
| Screened here | 146, of which 145 were queried fresh and 1 was reused from cache (SC_GMFUL5306400, the pilot negative control) |
| Markers | the 24 verified markers (`results/audit/marker_verification.csv`: 24/24 PASS); hierarchy M132/Z958/CTS246/Z906 → 5-SNP E-Z15174 block → Z5991/FT402211/Z15196/FTF75935 + Tier-A |
| Read rules | primary, non-duplicate, MAPQ ≥ 20, BQ ≥ 20, ≥ 50% aligned, AS ≥ XS; depth in **independent fragments** (overlapping mates count once) |
| Calls | ≥ 2 fragments and ≥ 90% agreement; 0–1 fragment = **no-call (unknown)**; never imputed as ancestral |
| Node rules | derived if ≥ 1 SNP D and 0 A; **outside only if ≥ 2 defining SNPs A** (single-SNP nodes flagged); D plus A = conflict |
| Queries | 709; **0 errors; 0 mixed; 0 conflicting nodes**; 62 no-calls (39 single fragment, 21 zero coverage, 2 one unusable read) |
| Transfer | **~0.22 GB** of CRAM slice data (computed from `.crai` slice sizes for every executed query); indexes were already cached; no complete CRAMs |
| Runtime | 72 min, sequential (~6 s per query incl. network), per-query JSON checkpoints |

## 2. Observed placements: sample counts (all 156 men)

| Placement | Fula | Mandinka | Wolof | Total |
|---|---|---|---|---|
| E-Z15174 | 0 | 0 | 0 | **0** |
| E-CTS246, not E-Z15174 (sister lineages) | 0 | 5 | 2 | **7** |
| E-Z958, not E-CTS246 | 2 | 4 | 0 | 6 |
| E-Z958, E-CTS246 unresolved (E-Z15174 excluded by 4–5 block SNPs) | 1 | 1 | 0 | 2 |
| E-M132, not E-Z958 (Z958 ancestral, single SNP; CTS246/Z906 also ancestral) | 0 | 2 | 0 | 2 |
| Not E-M132 (M132 ancestral; single defining SNP) | 48 | 37 | 52 | 137 |
| Unresolved at E-M132 (E-Z15174 excluded for SC_GMWOF5428851; unknown for SC_GMWOF5428755) | 0 | 0 | 2 | 2 |
| **Total** | 51 | 49 | 56 | 156 |

These are **counts of sampled men**, nothing more. 10 Mandinka and Wolof men came from the coverage-selected
pilot; the other 146 are the complete remaining set.

## 3. Sample proportions (separate from counts, and not population frequencies)

| Node | Fula | Mandinka | Wolof |
|---|---|---|---|
| E-M132 | 3/51 = 5.9% (2.0–15.9) | 12/49 = 24.5% (14.6–38.1) | 2/54 = 3.7% (1.0–12.5) |
| E-CTS246 | 0/50 = 0% (0–7.1) | 5/48 = 10.4% (4.5–22.2) | 2/54 = 3.7% (1.0–12.5) |
| E-Z15174 | 0/50 = 0% (0–7.1) | 0/48 = 0% (0–7.4) | 0/54 = 0% (0–6.6) |

Wilson 95% intervals; denominators are men with a resolved status at that node. **These are not population
frequencies:**
* GGVP sampled families and particular communities in the Western Division, so relatives may share a Y and
  the effective number of independent lineages is smaller than n, making the intervals too narrow;
* labels are the dataset's sampling labels;
* 10 men were coverage-selected.

They describe this dataset only.

## 4. The 7 E-CTS246 sister-lineage men (flagged: subclade undetermined)

| Sample | Label | M132 | Z958 | CTS246 | Z906 | E-Z15174 block | Run |
|---|---|---|---|---|---|---|---|
| SC_GMMAN5482181 | GGVP Mandinka | D 9 (6/3) | D 9 (4/5) | D 3 (2/1) | D 6 (6/0) | 5/5 A | screen |
| SC_GMMAN5482185 | GGVP Mandinka | D 2 (2/0) | N (1 frag) | D 5 (3/2) | D 5 (2/3) | 5/5 A | pilot |
| SC_GMMAN5482191 | GGVP Mandinka | D 8 (4/4) | D 8 (2/6) | D 4 (4/0) | D 12 (7/5) | 5/5 A | screen |
| SC_GMMAN5482207 | GGVP Mandinka | D 4 (2/2) | D 7 (2/5) | D 9 (5/4) | D 4 (3/1) | 5/5 A | pilot |
| SC_GMMAN5482299 | GGVP Mandinka | D 8 (7/1) | D 5 (4/1) | D 8 (5/3) | D 7 (4/3) | 5/5 A | screen |
| SC_GMWOF5428789 | GGVP Wolof | D 9 (5/4) | D 5 (0/5) | D 3 (2/1) | D 6 (4/2) | 5/5 A | screen |
| SC_GMWOF5428805 | GGVP Wolof | D 7 (2/5) | D 5 (2/3) | D 6 (3/3) | D 7 (4/3) | 5/5 A | screen |

Cells: call, derived fragments, (fwd/rev). Every man has two concordant E-CTS246 SNPs and five ancestral
E-Z15174 SNPs.

**None is on YFull's tree**, so their placement within E-CTS246 is **unresolved**. They could belong to known
sister clades (E-CTS3380 > E-Z5993, like 1000G GWD HG03048; E-CTS10935; E-L133) or to unsampled ones.
Some may be close relatives of each other, since GGVP includes families. They are **flagged for verification**
with sister-clade SNPs; the five Mandinka men are the larger set.

**Other flags:**
* SC_GMWOF5428755: E-Z15174 unknown (coverage).
* SC_GMWOF5428851: E-M132 unresolved (M132 single fragment); E-Z15174 excluded.
* SC_GMFUL5306398 and SC_GMMAN5482307: E-CTS246 unresolved (1 A + 1 no-call); E-Z15174 excluded.

## 5. Consistency checks

* **YFull concordance:** the three GGVP Fula men placed by YFull in E-Z5985 (E-Z958, not E-CTS246) are
  classified consistently. SC_GMFUL5306392 and SC_GMFUL5306400 are E-Z958 and not E-CTS246 (2/2 A).
  SC_GMFUL5306398 is E-Z958 with E-CTS246 unresolved (Z906 A, CTS246 1 fragment) and E-Z15174 excluded.
* **FTF75935 screen concordance:** every man here was already ancestral or a no-call at FTF75935. None is
  E-Z15174, so none could carry FTF75935 consistently with the tree. No contradiction.
* **Nesting:** no man is derived at a child node while ancestral at its parent (0 conflicts).

## 6. Updated provisional tree (GGVP non-Jola men in **bold**)

```
E-M132                                    **2 GGVP Mandinka: E-M132, not E-Z958** (Z958 single-SNP A)
└── E-Z958                                **6 men (2 Fula, 4 Mandinka) not E-CTS246**; **2 unresolved at CTS246**
    │                                     1000G GWD/MSL; GGVP Fula 5306392/398/400 = YFull E-Z5985
    └── E-CTS246
        ├── known sister clades           1000G MSL, GWD (HG03048), YRI, ACB, PUR men
        ├── **7 GGVP men, subclade undetermined**: Mandinka 5482181, 5482185, 5482191, 5482207,
        │                                  5482299; Wolof 5428789, 5428805
        └── E-Z15174                      no GGVP Fula/Mandinka/Wolof man (0/155 resolved; 1 unknown)
            ├── E-Z15174* / E-FTF75935    YF133256 (+ IN147690 reported, unverified)
            └── E-Z5991 > E-FT402211 (GGVP Jola 5309804, 5309851) | E-Z15196 (GGVP Jola 5309829; GWD HG03024, HG03027)
```

## 7. What this does and does not mean

* **Does:** in the GGVP collection, E-Z15174 has now been tested in every Fula, Mandinka and Wolof male, and
  it appears only among Jola-sampled men (3, placed by YFull and genotyped here). In 1000G, it appears in two
  GWD Mandinka men.
* **Does not:**
  * show that E-Z15174 is "Jola" or absent from other groups. Zero carriers in about 50 men per group is
    compatible with a lineage carried by a few percent of men, and the observed carriers are confounded with
    where sampling and YFull processing were densest;
  * say anything about Fula men outside the GGVP sampling, or about IN147690;
  * assign an ethnicity to your 18th-century ancestor. A lineage whose nearest public relatives share an
    ancestor ~7,000 years ago (E-Z5991 men) or ~11,000 years ago (the E-CTS246 sister men here) does not map
    onto present-day ethnic groups.

## 8. Three most informative remaining research questions

1. **Where exactly is IN147690, and are any of your 19 "private" SNPs actually shared with him?** Only the
   FTDNA E-FTF75935 Block Tree (or his Big Y data, with consent) can answer this. Together with your Big Y VCF
   (which also confirms your Z5991-ancestral state from your own reads), it is the only evidence that bears on
   a *recent* paternal relative. No public screen can replace it.
2. **Which E-CTS246 sister clades do the 7 newly found GGVP men belong to: known YFull clades, or unplaced
   lineages?** Testing 3–5 verified SNPs for E-CTS10935, E-CTS3380/E-Z5993 and E-L133 on these 7 men
   (~35 queries, a few minutes) would show whether western Gambia holds E-CTS246 diversity not yet on the
   tree. That bears on where and when E-Z15174's parent lineage diversified.
3. **Is E-Z15174 really restricted to western-Gambian Jola- and Mandinka-sampled men?**
   * Within GGVP, the 44 Jola men not placed by YFull at E-Z15174 have only been screened with the
     September 2026 four-marker scan (one marker mislabelled, many single-read calls). Re-screening them with
     this pipeline would firm up the only within-dataset comparison.
   * Beyond GGVP, the HGDP Mandenka (Senegal; 18 not yet tested) and other West African datasets would test
     the distribution outside the Gambia.

Each item needs your approval before any data are accessed.
