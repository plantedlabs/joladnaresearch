# Phase 1: GGVP Jola upstream re-screen, with all 203 GGVP males combined

**Date:** 2026-11-10 · **Phase 2 (HGDP Mandenka) planned only, not executed** (§7)

**Deliverables:**
* `results/audit/jola_male_reconciliation.csv`: all 100 Jola samples: sex, YFull status, earlier screens
* `results/jola_screen/pilot_marker_evidence_long.csv`: 196 queries, with fragments, strands, MAPQ, BQ and
  no-call reasons
* `results/jola_screen/raw_reads/`: raw-read validation SAMs
* `results/ggvp_all/ggvp_all_sample_by_marker.csv`: **sample × marker matrix for all 203 GGVP males**
* `results/ggvp_all/ggvp_all_males_placements.csv`: placement and E-Z15174 basis for each man
* `results/ggvp_all/ggvp_branch_distribution.md`: counts and caveated proportions
* `results/topology/provisional_tree.md`: updated tree
* Scripts: `scripts/12_branch_pilot.py`, `scripts/16_summarize_ggvp_all.py`

---

## Answer

* **No additional Jola man carries E-Z15174.** Of the 47 GGVP Jola-sampled males, exactly **3** are E-Z15174:
  the three already known (SC_GMJOL5309804 and 5309851 on E-FT402211; 5309829 on E-Z15196). The other
  **44 are excluded on ≥ 2 independent ancestral SNPs each**: 4 at the E-Z15174 block, and 40 at both E-CTS246
  SNPs. No Jola man is unknown.
* **No Jola man sits closer to your lineage.** All three carriers are on E-Z5991, the sister branch of your
  E-Z15174* / E-FTF75935 line, and all are ancestral at FTF75935. None of them shares any of your Tier-A
  candidate SNPs (earlier work). No man anywhere in GGVP is on an E-Z15174 branch other than E-Z5991.
* Across **all 203 GGVP males**, E-Z15174 is derived in 3, excluded in 199 (each on ≥ 2 SNPs), and unknown in
  1 (a Wolof man with very low coverage).

## 1. Reconciling the 47 Jola males

| Group | n | Basis |
|---|---|---|
| Jola samples / inferred male | 100 / **47** | sex window (females ≤ 56 reads, males ≥ 890) |
| Already genotyped here (13 tree markers, strict) | 3 | 5309804, 5309851, 5309829 |
| On YFull's E-M132 tree, never genotyped by this pipeline | **23** | YFull places all of them in E-Z958 > E-Z5985 sub-branches (E-FT402321, E-MF28713, E-Y318152, …) |
| Not on YFull's E-M132 tree, never genotyped | **21** | either non-E-M132 (YFull does not list) or never processed by YFull |
| September 2026 repo screen | 46 IDs listed | 4 markers, one mislabelled (Z15197 called "Z15234"), many 1-read calls; **missed SC_GMJOL5309891** (male, on YFull E-MF28713) |

All 26 YFull-listed Jola IDs are male in our inference, and all 46 repo IDs are male. **44 men were screened**:
the 23 YFull-placed men serve as an external concordance test, and the 21 others are the genuinely
unclassified men.

## 2. Screen results (44 Jola men)

| Placement | YFull-placed (23) | Not on YFull (21) |
|---|---|---|
| E-Z958, not E-CTS246 | 21 | **2 (SC_GMJOL5309950, 5309952)** |
| E-Z958, E-CTS246 unresolved (E-Z15174 excluded by block) | 1 (5309844) | 0 |
| E-M132, Z958 unresolved (E-CTS246 ancestral 2/2) | 1 (5309973) | 0 |
| Not E-M132 | 0 | 18 |
| Unresolved at E-M132 (E-CTS246 ancestral 2/2) | 0 | 1 (5309898) |
| **E-Z15174** | **0** | **0** |

* 196 queries; 0 errors; 10 no-calls (single fragment or zero coverage); 1 mixed call.
* The mixed call is SC_GMJOL5309926 at Z958: 4 ref and 1 alt (one strand). He is M132-ancestral, so this is
  noise, not a placement.
* At the 20 block-SNP queries, there was **no derived read at all**, usable or filtered.
* **YFull concordance: 0 contradictions among 23.** 21 match exactly. The other 2 are only partially resolved
  because of a single-fragment no-call; their observed calls are compatible with YFull.
* **Raw-read validation of the unexpected placements.** SC_GMJOL5309950 and SC_GMJOL5309952 are E-M132 men
  absent from YFull's tree.

  | Man | Site | Reads | Alignment | AS vs XS | Both strands |
  |---|---|---|---|---|---|
  | 5309950 | M132 | 4 | all 100M, MAPQ 60 | 90–98 vs 19–20 | yes |
  | 5309950 | Z958 | 6 | all 100M, MAPQ 60 | 95–96 vs 0–24 | yes |
  | 5309952 | M132 | 7 | all 100M, MAPQ 60 | 90–95 vs ≤ 20 | yes |
  | 5309952 | Z958 | 5 | all 100M, MAPQ 60 | 95 vs ≤ 24 | yes |

  Reads come from several sequencing runs. CTS246 and Z906 are ancestral in both (2–7 fragments). These are
  reliable E-Z958 (not E-CTS246) placements, which suggests YFull has not processed every GGVP Jola sample.

## 3. Branch distribution across all 203 GGVP males: sample counts

| Placement | Fula | Jola | Mandinka | Wolof | Total |
|---|---|---|---|---|---|
| E-Z15174 > E-Z5991 > E-FT402211 | 0 | 2 | 0 | 0 | 2 |
| E-Z15174 > E-Z5991 > E-Z15196 | 0 | 1 | 0 | 0 | 1 |
| E-CTS246, not E-Z15174 | 0 | 0 | 5 | 2 | 7 |
| E-Z958, not E-CTS246 | 2 | 23 | 4 | 0 | 29 |
| E-Z958, E-CTS246 unresolved | 1 | 1 | 1 | 0 | 3 |
| E-M132, not E-Z958 | 0 | 0 | 2 | 0 | 2 |
| E-M132, E-Z958 unresolved | 0 | 1 | 0 | 0 | 1 |
| Not E-M132 | 48 | 18 | 37 | 52 | 155 |
| Unresolved | 0 | 1 | 0 | 2 | 3 |
| **Men** | **51** | **47** | **49** | **56** | **203** |

E-Z15174 status basis:

| Basis | Fula | Jola | Mandinka | Wolof |
|---|---|---|---|---|
| Derived | 0 | 3 | 0 | 0 |
| Excluded, ≥ 2 block SNPs ancestral | 11 | 4 | 6 | 6 |
| Excluded, E-CTS246 ancestral at 2/2 | 40 | 40 | 43 | 49 |
| Single-SNP exclusion only | 0 | 0 | 0 | 0 |
| Unknown | 0 | 0 | 0 | 1 |

### Sample proportions (separate from counts; **not** population frequencies)

| Node | Fula | Jola | Mandinka | Wolof |
|---|---|---|---|---|
| E-M132 | 3/51 = 5.9% (2.0–15.9) | 28/46 = 60.9% (46.5–73.6) | 12/49 = 24.5% (14.6–38.1) | 2/54 = 3.7% (1.0–12.5) |
| E-CTS246 | 0/50 = 0% (0–7.1) | 3/45 = 6.7% (2.3–17.9) | 5/48 = 10.4% (4.5–22.2) | 2/54 = 3.7% (1.0–12.5) |
| E-Z15174 | 0/51 = 0% (0–7.0) | 3/47 = 6.4% (2.2–17.2) | 0/49 = 0% (0–7.3) | 0/55 = 0% (0–6.5) |

Wilson 95% intervals, describing **this dataset only**:
* GGVP sampled families and particular Western Division communities, so relatives can share a Y chromosome
  and the effective number of lineages is smaller than *n*. The intervals are therefore too narrow.
* SC_GMJOL5309804 and 5309851 share an ancestor only ~375 years ago, and the 23 Jola E-Z5985 men include
  YFull clusters with TMRCAs of 50–150 years. These are close paternal relatives, not independent draws.

The earlier repo figure of "3/47 Jola males" is now confirmed with read-level evidence for all 47 men. It
remains a count in this collection, not a frequency for any people.

## 4. Updated provisional tree (GGVP men, observed; YFull-reported where marked)

```
E-M132
├── (not E-Z958)            2 GGVP Mandinka
└── E-Z958
    ├── (not E-CTS246)       29 GGVP men: Jola 23 (21 = YFull E-Z5985 sub-branches, + 5309950, 5309952 new),
    │                        Mandinka 4, Fula 2 (5306392/5306400 = YFull E-Z5985); + 3 with CTS246 unresolved
    └── E-CTS246
        ├── sister lineages  7 GGVP men: Mandinka 5482181/185/191/207/299; Wolof 5428789/805 (subclade undetermined)
        │                    + 1000G MSL/GWD/YRI/ACB/PUR men in E-CTS10935 / E-CTS3380 / E-L133 [YFull + VCF]
        └── E-Z15174
            ├── E-Z15174* / E-FTF75935   YF133256 [reported]; 0/171 callable GGVP males and no 1000G carrier at FTF75935
            └── E-Z5991
                ├── E-FT402211           GGVP Jola SC_GMJOL5309804, 5309851
                └── E-Z15196             GGVP Jola SC_GMJOL5309829; 1000G GWD HG03024, HG03027 (+HG03029)
```

## 5. Interpretation: sampling labels vs inferred ethnic origin

* **Sampling labels:** "GGVP Jola", "GGVP Mandinka" and so on are the cohort's recruitment categories.
  YFull's language tags ('dyo', 'mnk', 'ful') are another layer, and they do not always agree. Two
  GGVP-"Jola" men in E-Z5985 are tagged 'ful' by YFull.
* **Observed fact:**
  * Within GGVP, E-Z15174 occurs only among Jola-sampled men (3 of 47), all on E-Z5991.
  * In 1000G it occurs in two Gambian Mandinka men, also on E-Z5991 (E-Z15196). E-Z15196 itself joins a
    Jola-sampled man and Mandinka-sampled men within ~1,000 years.
* **What can be inferred:**
  * The known E-Z15174 diversity in public data is concentrated in western-Gambian men of at least two
    sampling labels.
  * Your line split from all of them ~7,000 years ago (YFull TMRCA of E-Z15174). That split predates any
    present-day ethnic group, so **no ethnic origin can be inferred for your ancestor** from these data.
  * "Jola-associated" overstates the data. Carriers are few, partly closely related (two of them share an
    ancestor ~375 years ago), and drawn from a cohort where Jola men were sampled and YFull-processed
    alongside Mandinka men.
  * **"Broader western Senegambian, unresolved"** remains the defensible summary.

## 6. Evidence strength

| Finding | Strength |
|---|---|
| No additional GGVP Jola E-Z15174 carrier | **Strong**: 44/44 excluded on ≥ 2 SNPs |
| No GGVP man of any label on an E-Z15174 branch other than E-Z5991 | **Strong**: 199/203 excluded, 3 on E-Z5991, 1 unknown |
| Pipeline agrees with YFull | **Strong**: 0 contradictions in 29 placements |
| SC_GMJOL5309950 and 5309952 are E-Z958, not E-CTS246 | **Strong**: raw reads validated |
| 7 E-CTS246 sister-lineage men (Mandinka, Wolof) | **Strong** for the node; subclade unresolved |

## 7. Phase 2 plan: HGDP Mandenka (not executed; needs approval)

**Samples.** 20 HGDP Mandenka (Senegal) high-coverage CRAMs on EBI; the reference is identical (chrY M5
verified). Already done: HGDP00904 and HGDP00905 are classified "not E-M132". HGDP00906–00908 were found
male by the sex-window count only, and HGDP00909 is female-like. **Remaining for classification: 18**
(00906–00919 range plus 01200–01202, 01283, 01285).

**Sex inference adapted to 30× coverage.** Fixed read-count thresholds set at GGVP's 4× do not transfer:
female HGDP00909 has 102 Y-window reads, above the GGVP female cut-off of 60. Plan:
1. Count primary, non-dup, MAPQ ≥ 20 reads in the Y window chrY:7,000,001–7,100,000 **and** in two
   100%-unique X windows (chrX:70,000,001–70,100,000 and chrX:100,000,001–100,100,000; Umap k100 = 1.0).
   chrX is already MD5-verified locally.
2. Use **R = Y / mean(X)**. Expected R ≈ 0.5 for XY (one Y versus one X), and ≈ 0 for XX, with residual
   mismapping well below 0.05.
3. Thresholds: male R ≥ 0.25, female R ≤ 0.05, and anything else ambiguous (not classified). Before use,
   calibrate on the 6 HGDP samples already counted, plus a handful of GGVP controls, and report each
   sample's R.
4. As an independent check, compare against the published HGDP sex where available (HGDP sample sheets;
   it is not in the EBI index). If that would need a new host, ask first.

**Screening.** The same 24 verified markers and hierarchy (`12_branch_pilot.py`), with coverage-appropriate
call rules. At 30× a call needs ≥ 2 fragments as before. Report confidence, because most sites will have
10–20 fragments. Validate any E-CTS246 or E-Z15174 hit at raw-read level.

**Cost, estimated from the 6 cached HGDP indexes:** per sample ≈ 1.5 MB index + 1–2 MB sex windows + 1.3–3.7 MB
step 1 (+1.6–4.7 MB block if needed) ≈ 6–12 MB. For 18 samples that is **≈ 0.1–0.2 GB and ~15–20 min**.
No full CRAMs.

**What it can and cannot show.** With about 10–15 males expected, a zero result would bound E-Z15174 in this
Senegalese Mandenka sample only weakly (95% upper bound ≈ 20–25%). A positive would be the first public
E-Z15174 carrier outside the Gambia, so it is informative mainly if positive.
