# E-Z15174 among Fulani and neighbouring West Africans: dataset discovery and feasibility

**Date:** 2026-11-12 · **Status:** Phase 1 (discovery) and Phase 2 (feasibility) complete. Phase 3 is proposed
only; **no data were downloaded or screened** in this pass beyond metadata and published supplementary tables.

**Lineage:** E-M132 → E-Z958 → E-CTS246 → E-Z15174 → (FTDNA) E-FTF75935. On YFull you are E-Z15174*.

**Markers and builds (already verified; `results/audit/marker_verification.csv`):**

| Marker | GRCh38 | Ref>Alt | GRCh37 cross-check |
|---|---|---|---|
| M132 | 19,734,375 | G>T | 21,896,261 G>T |
| CTS246 (= V1119) | 2,972,260 | G>A | 2,840,301 G>A |
| Z906 | 8,360,264 | T>C | – |
| Z15174 / Z15197 / Z15231 / Z15232 / Z15234 (E-Z15174 block) | 4,202,358 / 8,414,659 / 14,719,978 / 14,765,245 / 15,274,693 | (see audit file) | – |
| FTF75935 | 2,895,562 | A>C | – |

The GRCh37 positions and alleles for M132 and CTS246 were found independently in D'Atanasio et al. 2018,
Table S2. They agree with YBrowse.

**Already done (not repeated here):**
* all 203 GGVP males (Fula, Jola, Mandinka, Wolof) at E-M132 → E-Z15174, plus FTF75935
* 1000 Genomes high-coverage chrY calls for all 3,202 samples, including **MSL** (Mende, Sierra Leone) and GWD

---

## 1. What was checked (discovery log)

| Source | What was verified | Finding |
|---|---|---|
| **Fortes-Lima et al. 2025, AJHG 112:261** ("Population history and admixture of the Fulani people from the Sahel"); GitHub `Schlebusch-lab/Sahel_study`; EGA metadata API | Populations, sample counts, data type, access | 452 Fulani in 18 groups: Mauritania 16, Senegal 18 + Halpularen 23, **Guinea 25**, Mali 22 + 16, Burkina Faso 19 + 24 + 54, Niger 5 groups (91), Chad 23 + 17, Cameroon 30, Gambia 74 (= GGVP). New data: **EGAD50000000653** (90 saliva samples, Senegal / Mauritania / Mali) and **EGAD50000000654** (329 WGA samples). Both are **Illumina H3Africa array** (H3Africa_2019_20037295_B1, ~2.27 M SNPs), **controlled access** (Sahel DAC EGAC50000000290; policy EGAP50000000253). The Guinea and Senegal-Halpularen samples come from earlier published data; the repo does not say which. **Not sequencing.** Whether the array has probes at CTS246, the Z15174 block or FTF75935 is **unverified** (needs the Illumina manifest). Probes for SNPs this rare are unlikely |
| **D'Atanasio et al. 2023** ("Echoes from the last Green Sahara", bioRxiv 2023.04.06.535569; EGA **EGAS00001007499 / EGAD00001011812**) | Sample table (preprint Table 1), access | 43 genomes at 30×, 23 of them Fulani: Benin 1, Burkina Faso 2, Cameroon 3, Ghana 5, Niger 3, Nigeria 4, **Senegal 3**, Chad 1, USA-sampled Fulani 1. Also Senegal Wolof 2 and others. **No Guinea, Sierra Leone or Mali.** Sex not listed. **Controlled access** (EGAP00001003524). Ideal Y coverage. The paper's Y discussion concerns E-M2/Z15939, not E1a |
| **D'Atanasio et al. 2018**, Genome Biol 19:20 (Y capture of 104 men; 142 markers in ~8,000 men from 145 populations); ENA **PRJEB24071** (public BAMs) | Supplementary Tables S1, S2, S5 (via Europe PMC) | **No E1a men in the 104-sample lab set.** The only E1a is the public YRI NA19239. The 145-population frequency survey (Table S5), which includes Fulbe of Niger, Nigeria and Cameroon/Chad, Mandenka of Senegal, GWD and MSL, covers **only A3-M13, E-M2, E-M78 and R-V88**. It has **no E1a/E-M132 data**. CTS246 (GRCh37 2,840,301) is inside capture block 1, but there are no E1a samples to look at. **No value for E-Z15174** |
| **Trombetta et al. 2015**, GBE 7:1940 | Scope | E-M35 (E1b1b) only. Not relevant |
| **Bučková, Černý & Novelletto 2013**, AJPA; **Cruciani et al. 2002 / Rosa et al. 2007** (as summarised) | E1a in Fulani | Low-resolution Y-SNP surveys. E-M132 is reported in some Fulani/Fulbe samples (e.g. Sudan "12%", Burkina Faso Fulbe "10%" of n = 20, Cameroon Fulbe "53%"), but these are second-hand figures not checked against the primary tables. **None resolves CTS246 or Z15174** |
| **SGDP** (Mallick 2016; EBI GRCh38DH CRAMs, public) | Population list from EBI sequence index | West Africa: Mandenka 3, Mende 2, "Gambian" 2, Yoruba 3, Esan 2. **No Fulani.** Sex not yet checked. Possible overlap with HGDP Mandenka and 1000G GWD/MSL individuals (not verified) |
| **HGDP high-coverage** (Bergström 2020; EBI GRCh38DH, public) | Alignment index; reference M5 | **20 Mandenka (Senegal)** CRAMs, identical chrY reference (M5 verified). No Fulani. YFull lists one HGDP Mandenka (HGDP01284) on a non-Z958 E-M132 branch. A screening plan already exists (`JOLA_RESCREEN_REPORT.md` §7) |
| **1000 Genomes MSL** | NYGC 3202 chrY VCF (already queried) | All MSL males covered. E-CTS246: HG03096 (+ son), placed by YFull in E-CTS10935. E-Z958: HG03439. **E-Z15174: none.** Nothing new to gain |
| **"Sudan Fulani WGS study"** | Web and EGA search | **Not located.** Found instead: a 2026 PNAS WGS of 125 Sudanese from 5 groups (group list not confirmed to include Fulani; no accession seen); Hollfelder 2017 arrays (Sudan, includes Hausa); SudanMitoSeq (mtDNA only). No Sudan Fulani WGS with an accession could be confirmed |
| **H3Africa high-depth WGS** (Choudhury 2020, 426 genomes, 50 groups, 13 countries) | Literature | Includes Guinea (Soussou). Whether any Fulani are included is not confirmed. **Controlled** (H3Africa DAC) |
| EGA free-text search | API test | Not usable: the query parameter is ignored. Accessions were found through the literature |

## 2. Feasibility table

Males are counted from published or metadata information where available. "?" means not stated.

| Dataset | Population(s) | Males | Data type | Y coverage for our SNPs | Access | Est. cost | Expected research value |
|---|---|---|---|---|---|---|---|
| **D'Atanasio 2023 Fulani WGS** (EGAD00001011812) | Fulani of 9 countries (23); Wolof, Songhai etc. | ? (≤ 43; sex not listed) | WGS 30× | **Excellent** (all markers, including FTF75935 and the E-Z15174 block) | **Controlled** (EGA DAC; institutional application) | Trivial once accessed (~10 sites × 43 samples) | **Moderate–high per sample**: multi-country Fulani at full Y resolution. No Guinea/SL/Mali |
| **Fortes-Lima 2025 Fulani arrays** (EGAD50000000653/654 + published Guinea/Senegal sets) | Fulani, 18 groups incl. **Guinea 25**, Senegal, Mauritania, Mali | ~half of 452 (?) | H3Africa array | **Unverified, probably poor**: only SNPs on the array; Z15174 block and FTF75935 unlikely to be probed; M132/CTS246 unknown | **Controlled** (Sahel DAC) | Trivial once accessed | **High N and the right geography**, but likely only E-M132-level resolution. **Check the array manifest first** |
| **HGDP Mandenka** | Mandenka (Senegal) | ~10–15 expected of 20 (5 of 6 checked were male) | WGS ~30× | Excellent | **Public** (EBI) | ~0.1–0.2 GB, ~15–20 min | **Low–moderate**: only public, high-coverage, non-Gambian Senegambian set left. Not Fulani |
| **SGDP West Africa** | Mandenka 3, Mende 2, Gambian 2 (+ Yoruba 3, Esan 2) | ? | WGS ~30× | Excellent | Public (EBI) | < 0.1 GB | **Low**: tiny, and likely overlaps HGDP and 1000G |
| **H3Africa high-depth WGS** | 50 groups incl. Guinea Soussou | ? | WGS 30× | Excellent | Controlled (H3Africa DAC) | Trivial once accessed | Moderate if Guinea/Senegal men are present; Fulani content unconfirmed |
| D'Atanasio 2018 Y capture (PRJEB24071) | 104 men, Sahel/Africa | 104 | Y capture | Covers CTS246, but no E1a men | Public | – | **None** (ruled out) |
| 1000G MSL / GWD | Mende, Mandinka | ~40 / 93 | WGS 30× | Excellent | Public | – | **None new** (already screened) |
| GGVP Fula / Mandinka / Wolof / Jola | Gambia | 203 | WGS ~4× | Adequate | Public | – | **None new** (already screened) |
| Published Y-SNP frequency surveys | Fulani, Sahel | hundreds | SNP typing | Only to E-M132/M33 at most | Public (papers) | none | Context only. Cannot resolve E-Z15174 |

## 3. The three strongest opportunities

1. **D'Atanasio 2023 Fulani WGS (EGAD00001011812), through the authors (Cruciani lab, Sapienza University of
   Rome) rather than a data request.** It has the right population, full Y resolution, and a lab that
   specialises in African Y phylogeny. Checking whether any of 23 Fulani (and 2 Senegal Wolof) are E-M132,
   E-CTS246 or E-Z15174 is a few-minute task for anyone with access. Limitation: no Guinea, Sierra Leone or
   Mali samples.
2. **Fortes-Lima 2025 Fulani array data (Sahel DAC), but only after confirming the H3Africa array's chrY
   content.**
   * It is the largest and most relevant sample: about 450 Fulani across the western Sahel, including 25
     from Guinea.
   * If the array carries M132/CTS246 probes, it would give the first Fulani-wide E-M132/E-CTS246
     distribution. It is unlikely to resolve E-Z15174 itself.
   * First step: check the public Illumina manifest for the GRCh37/38 positions above. This costs nothing
     and needs no data access.
3. **HGDP Mandenka (public), via the pilot already designed.** It is the only remaining public, freely
   screenable, high-coverage Senegambian set. Low expected yield for the Fulani question, but it is the one
   thing that can be done now without any access application.

Outside these datasets, the most informative evidence on your lineage's *recent* origin is still the
**FTDNA E-FTF75935 Block Tree, kit IN147690 (reported Gambian Fulani), and FTDNA/YFull testers from Guinea,
Sierra Leone and Senegal**. Consumer Big Y/WGS data resolve FTF75935 directly, and that is where any recent
branch-mates would appear. No academic dataset found here resolves your branch below E-Z15174.

## 4. Conditional pilot (proposal only; not run)

Only one promising dataset is accessible without an application: **HGDP Mandenka**. Pilot as designed in
`JOLA_RESCREEN_REPORT.md` §7:
* coverage-normalised Y/X sex inference
* hierarchical screen M132/Z958 → CTS246/Z906 → E-Z15174 block → Z5991/FT402211/Z15196/FTF75935, plus
  Tier-A candidates if E-Z15174 is found
* ≈ 0.1–0.2 GB, ~15–20 min

Optionally add the 7 SGDP Mandenka/Mende/Gambian samples (< 0.1 GB), after checking for duplicate
individuals.

Calls follow the established rules:
* derived or ancestral needs ≥ 2 independent fragments and ≥ 90% agreement
* 0–1 fragment is a **no-call**, never ancestral
* mixed reads are **ambiguous**
* "outside a branch" requires ≥ 2 defining SNPs called ancestral

**Expected outcome (a judgement, not a calculation from data):** E-Z15174 is absent from 152 resolved
non-Jola GGVP men and 93 GWD men. With ~10–15 Mandenka males, the chance of finding a carrier is low. A
negative result would add little.

## 5. Answer to the final question

**Public datasets have essentially reached diminishing returns for E-Z15174.**
* Every large public West African genome set that resolves your branch has been screened:
  * GGVP (all four groups)
  * 1000 Genomes GWD and MSL
* What remains public is small:
  * HGDP Mandenka (~20)
  * SGDP (~7 West Africans), probably overlapping
* The published Fulani Y surveys stop at E-M132-level resolution or study other haplogroups entirely
  (E-M2, E-M78, R-V88). None reports E-CTS246 or E-Z15174 in any Fulani population.

**A realistic opportunity for new evidence does exist, but it is not in open data:**
* (a) **Collaboration requests** to the two groups holding controlled Fulani data:
  * Cruciani/D'Atanasio, Rome: 30× WGS of 23 Fulani
  * Schlebusch/Černý, Uppsala and Prague: about 450 Fulani on H3Africa arrays; check the array manifest
    first
* (b) **Consumer genetic-genealogy data**: the FTDNA E-FTF75935 Block Tree (IN147690), and recruiting Big Y
  testers from Guinea, Sierra Leone and Senegal Fulani communities through FTDNA projects.

These are the routes that could show whether E-Z15174, or your E-FTF75935 sub-branch, is present among
Fulani beyond the single reported Gambian Fulani tester.

Two things would make the HGDP Mandenka pilot worth running anyway:
* you want every free public option closed off;
* it is a prerequisite for a collaboration request ("we screened all public data").
