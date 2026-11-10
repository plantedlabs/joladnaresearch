# Population screen for FTF75935 (terminal SNP of YF133256) in public Gambian genomes

**Date:** 2026-11-10
**Target:** FTF75935, chrY:2,895,562 A>C (GRCh38; coordinate and allele from YBrowse; Umap k100 = 1.0, not
blacklisted, no segmental duplication)
**Question:** does any publicly sequenced Gambian man carry your terminal-branch mutation? Such a man would be a
paternal relative on your branch, much closer than the ~7,000-year E-Z15174 link to the Jola men already studied.

**Files:**
* `results/ftf75935_screen/ggvp_ftf75935_screen.csv`: every one of the 394 GGVP samples (individual level)
* `results/ftf75935_screen/gwd_ftf75935_from_nygc_calls.csv`: the 93 GWD males (from published calls; depth not assessed)
* `results/ftf75935_screen/population_summary.md`: population tables
* `results/ftf75935_screen/cache/`: per-sample raw results (JSON)
* Code: `scripts/07_screen_site.py` (screen), `scripts/08_compile_screen.py` (summary),
  `scripts/09_confirm_candidates.sh` (confirmation, not triggered)

---

## Answer

**No Fula, Jola, Mandinka or Wolof man in the GGVP, and no Gambian Mandinka man in 1000 Genomes, shows any
evidence of FTF75935.**

* **GGVP:** 203 confidently male samples screened. 172 are callable at FTF75935, and **171 carry the ancestral
  A**. One call is mixed (A plus one G read, no C). **0 are derived.** Not a single read carrying C was observed
  in any man, including reads removed by the quality filters.
* **1000 Genomes GWD (Gambian Mandinka), high coverage:** the published NYGC chrY calls for all 3,202 samples
  have **no variant record at FTF75935**. That covers all 93 GWD males and ~1,233 males worldwide in the
  2,504-sample set. The surrounding ±5 kb is densely called (79–98 records, including singletons), so a carrier
  would very likely have appeared. The only E-Z15174 men in the whole collection are HG03024, HG03027 and
  HG03029 (HG03027's son), all on E-Z15196.
* **No candidates**, so no confirmation testing at the tree markers and Tier-A private SNPs was needed.
  **No previously unrecognised paternal relationship was found.**

---

## 1. Samples and sex inference (GGVP)

GGVP publishes no sex metadata, so sex was inferred from primary, non-duplicate, MAPQ ≥ 20 reads in a 100 kb
chrY window that is entirely uniquely mappable (chrY:7,000,001–7,100,000, Umap k100). The groups separated
sharply: **females 0–56 reads, males 890–10,000+ reads**.

| Population (GGVP label) | Samples | Male | Female | Ambiguous |
|---|---|---|---|---|
| Fula | 95 | 51 | 44 | 0 |
| Jola | 100 | 47 | 53 | 0 |
| Mandinka | 99 | 49 | 50 | 0 |
| Wolof | 100 | 56 | 43 | 1 |
| **Total** | **394** | **203** | **190** | **1** |

The one ambiguous sample, SC_GMWOF5428740, had 63 reads. That is just above the female cut-off (≤ 60) and 14×
below the lowest male (890), so it is almost certainly female. It was not counted as a male. The 47 Jola males
match the count in the earlier repo screen.

## 2. FTF75935 genotypes in GGVP males

Read filters: primary, non-duplicate, MAPQ ≥ 20, base quality ≥ 20, ≥ 50% of the read aligned, no better
alignment elsewhere (AS ≥ XS). **Depth is counted in DNA fragments**, so overlapping mates of one pair count once.

Haploid calling rules:
* ≥ 2 fragments and ≥ 90% agreement are needed for a call.
* A man with 0 or 1 usable fragment is a **no-call, never "ancestral"**.

| Population | Males | Callable (≥ 2) | Ancestral (A) | Derived (C) | Mixed | Insufficient coverage | Median depth |
|---|---|---|---|---|---|---|---|
| Fula | 51 | 37 | 37 | **0** | 0 | 14 | 2 |
| Jola | 47 | 39 | 39 | **0** | 0 | 8 | 3 |
| Mandinka | 49 | 46 | 46 | **0** | 0 | 3 | 4 |
| Wolof | 56 | 50 | 49 | **0** | 1 | 6 | 4 |
| **Total** | **203** | **172** | **171** | **0** | **1** | **31** | **4** |

Depth distribution among the 203 males:

| Depth | 0 | 1 | 2 | 3 | 4–5 | 6–9 | ≥ 10 |
|---|---|---|---|---|---|---|---|
| Males | 5 | 26 | 42 | 24 | 67 | 35 | 4 |

* **Confidence of the 171 ancestral calls:** 96 high (≥ 4 fragments on both strands), 33 moderate, and
  42 low (2 fragments).
* **The 31 uncallable males:** 5 had no reads. 26 had a single read, and all 26 of those reads showed A, none C.
  These 31 men are **not** classed as ancestral.
* **The mixed call** (SC_GMWOF5428859, Wolof): 4 independent A fragments (MAPQ 60) plus one fragment carrying
  **G**, with both mates showing G at MAPQ 45–49. G is not the FTF75935 allele. This is most likely an error on
  one molecule, and it is not evidence of FTF75935.
* **Controls:** the three Jola men studied earlier (SC_GMJOL5309804, 5309851, 5309829) are again ancestral,
  with 3, 5 and 3 fragments.
* **Quality-control note:** the first pass counted overlapping mates of a pair as two reads. All 203 males
  were re-screened with fragment-level counting. Six men each lost one unit of depth, and no call changed.
  The first-pass results are kept in `cache_v1_superseded/`.

## 3. Additional variants shared with YF133256

None to report, because no man carried FTF75935. Of the Gambian men examined at your private SNPs (the three
Jola men in earlier work), none shares any of them.

## 4. What the negative result can and cannot establish

**It can establish:**
* Among 171 GGVP men with a confident call and 93 GWD men in the high-coverage call set, **none is on your
  E-FTF75935 branch**.
* With zero carriers among *n* independent callable men, FTF75935 frequency in the sampled populations is
  unlikely to exceed the exact one-sided 95% upper bound:

| Group | Callable men (n) | 95% upper bound on FTF75935 frequency |
|---|---|---|
| GGVP Fula | 37 | 7.8% |
| GGVP Jola | 39 | 7.4% |
| GGVP Mandinka | 46 | 6.3% |
| GGVP Wolof | 49 | 5.9% |
| GGVP all | 171 | 1.7% |
| GGVP + 1000G GWD unrelated males (55) | 226 | 1.3% |

These bounds assume the men are independent paternal lineages. GGVP sampled families, so some men may be
father and son, and the true number of independent lineages, and hence the power, is somewhat lower. The GWD
figure uses only the 55 unrelated GWD males, because the other 38 are relatives of men already counted.

**It cannot establish:**
* **That the lineage is absent from these peoples.** A lineage carried by about 1–2% of men, or concentrated in
  particular villages, lineages or regions not sampled, could easily be missed with ~40–50 men per group.
  GGVP and GWD sampled mainly the Western Division of The Gambia.
* **Anything about neighbouring populations** that were not sampled: other Senegambian and Guinea-Bissau
  peoples (for example Serer, Balanta, Bainuk, Manjak, Papel), Jola communities in Casamance, or Fula
  populations elsewhere.
* **Anything about the 31 GGVP men without a call**, or about the 1000 Genomes men's per-person depth. For the
  latter, the inference rests on the absence of a joint-called variant record, not on direct read counts. The
  GWD CRAM depth check was deliberately not run.
* **Anything about IN147690.** If he is on E-FTF75935, his absence from these public samples is unremarkable.
  This screen neither supports nor weakens a relationship with him.
* **Ethnic origin.** A single SNP, present or absent in about 50 sampled men per group, cannot assign your
  ancestor to any people. Carriers found in one population would show a shared paternal lineage, not a shared
  ethnicity.

## 5. Limitations

* GGVP is low coverage (median 4 fragments at this site). 42 ancestral calls rest on 2 fragments.
  42 + 33 calls are low or moderate confidence.
* Only one SNP was screened. A man on your branch would carry FTF75935 *and* its equivalent SNPs. Screening
  2–4 equivalents would rescue some of the 31 uncalled men and guard against a site-specific artifact. That
  needs the E-FTF75935 block SNP list from your FTDNA Block Tree or Big Y VCF.
* Sex inference is empirical. Thresholds of ≤ 60 and ≥ 300 reads were set from the pilot, and the observed gap
  (56 vs 890) is wide.
* GGVP alignments are from 2015 (bwakit 0.7.12, GRCh38DH). Re-alignment is not expected to matter at a site
  this well mapped.

## 6. Recommended follow-up

1. **FTDNA Block Tree for E-FTF75935.** Get the number of branch SNPs and each tester's private-variant count,
   including IN147690. This is the most informative step and costs nothing computationally. See `REPORT.md` §6.
2. **Your Big Y VCF/BED.** It supplies the FTF75935 equivalents, which allows a multi-SNP re-screen of the
   31 uncalled GGVP men (cached indexes make this ~5 min).
3. **Other public West African genomes on the same EBI server** (indexed access, no full downloads):
   * **HGDP Mandenka (Senegal):** 20 high-coverage CRAMs.
   * **Simons Genome Diversity Project** West African samples.
   * Each is a single-site screen of a few minutes, after the same sex check.
4. **Optional: GWD CRAM depth check** (93 men, ~15 min, ~1 GB) to replace "no record in the joint call set"
   with direct per-man depth. Deferred, as you requested.
5. **Controlled-access West African cohorts** (for example H3Africa / EGA datasets from Senegal and
   Guinea-Bissau) would need a data-access application.

---
*Data: Gambian Genome Variation Project (MRC Unit The Gambia, Wellcome Sanger Institute, MalariaGEN), accessed
through IGSR/EBI under the Fort Lauderdale principles. 1000 Genomes high-coverage chrY calls: New York Genome
Center.*
