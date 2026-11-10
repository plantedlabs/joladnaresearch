# Provisional E-M132 > E-Z958 > E-CTS246 > E-Z15174 tree (public samples)

Evidence layers:
* **[obs-CRAM]**: genotyped by this project from reads
* **[obs-VCF]**: genotyped from NYGC 3202 chrY joint calls
* **[YFull]**: placement reported by YFull v14.07 (02 Nov 2026), not re-genotyped here at that node's own SNPs
* **[FTDNA]**: reported by FamilyTreeDNA, not verified

Population labels are the **sampling labels** of each dataset: GGVP population, 1000G population code, or
YFull's language tag (dyo = Jola-Fonyi, mnk = Mandinka, ful = Fula). They are not inferred ancestry.

```
E-M132 (M132 = 19,734,375 G>T)
└── E-Z958 (Z958 = 20,558,142 C>T)                          1000G: 38/38 E-M132 men are E-Z958 [obs-VCF]
    ├── E-Z5985 and other non-CTS246 E-Z958 branches          [YFull]; Z958+ CTS246− [obs-VCF] for 18 1000G men
    │     1000G GWD: HG02461, HG02570, HG02623, HG02634, HG02810, HG02819, HG03246, HG03258 (+ sons); MSL HG03439
    │     GGVP Fula: SC_GMFUL5306392, 5306398, 5306400 [YFull, 'ful']
    │     GGVP Jola: ≥ 20 men incl. 5309827, 5309828 (YFull tag 'ful'), 5309844, 5309852, 5309903, 5309922 ... [YFull]
    └── E-CTS246 (CTS246 = 2,972,260 G>A; Z906 = 8,360,264 T>C)   carriers nested in Z958 [obs-VCF, obs-CRAM]
        ├── E-CTS10935 > E-FT327560 / E-Y610327     MSL HG03096 (+HG03098)          [YFull]; CTS246+ Z15174− [obs-VCF]
        ├── E-CTS3380 > E-Z5993                     GWD HG03048 (+HG03050)          [YFull]; CTS246+ Z15174− [obs-VCF]
        │             > E-FT296631                  PUR HG01094 (+HG01096)          [YFull]; CTS246+ Z15174− [obs-VCF]
        │             > E-P110                      YRI NA19239, NA19189 (+NA19191) [YFull]; CTS246+ Z15174− [obs-VCF]
        ├── E-L133 / E-CTS736 / E-Z5992             YRI NA18910 (+NA18911), NA18856 (+NA18857); ACB HG02449, HG01882 (+HG01888)
        │                                           [YFull]; CTS246+ Z15174− [obs-VCF]
        └── E-Z15174 (39 SNPs; tested: Z15174 4,202,358; Z15197 8,414,659; Z15231 14,719,978; Z15232 14,765,245; Z15234 15,274,693)
            │   all 5 tested block SNPs: same carriers in 1000G (3/3202) [obs-VCF]; derived in 3 Jola [obs-CRAM]
            ├── E-Z15174* ─ YF133256 (USA, African American) ── E-FTF75935 (2,895,562 A>C) [FTDNA/YFull]
            │       FTF75935: 0/171 callable GGVP males; no record in 3202 1000G [obs]
            │       IN147690 (reported Gambian Fulani): branch position unverified [FTDNA, not observed]
            └── E-Z5991 (Z5991 = 2,849,335 C>T)      YF133256 ancestral [YFull]; 5 public men derived [obs]
                ├── E-FT402211 (9,176,816 A>G)        GGVP Jola SC_GMJOL5309804 ('dyo'), SC_GMJOL5309851 [obs-CRAM]
                │                                     TMRCA ~375 ybp [YFull]
                └── E-Z15196 (8,336,326 A>C)          GGVP Jola SC_GMJOL5309829 ('dyo') [obs-CRAM]
                                                      1000G GWD HG03024, HG03027 (+ son HG03029) ('mnk') [obs-VCF]
                                                      TMRCA ~1,000 ybp [YFull]
```

Dates are YFull estimates (formed / TMRCA): E-CTS246 formed 12,000 / TMRCA 11,000; E-Z15174 formed 11,000 /
TMRCA 7,000; E-Z5991 formed 7,000 / TMRCA 2,300.

**Not yet tested at the E-CTS246 or E-Z15174 level:** GGVP Mandinka (49 males), Wolof (56) and most Fula (51)
males; HGDP Mandenka (20, high coverage); SGDP West Africans. YFull's tree contains no GGVP Mandinka or Wolof
samples, so it offers no placement for them.

---

## Update after the full GGVP upstream screen (2026-11-10)

All 156 GGVP Fula, Mandinka and Wolof males are now classified (`docs/UPSTREAM_SCREEN_REPORT.md`):

* **E-Z15174: 0 men.** 155 are excluded with ≥ 2 ancestral SNPs (23 at the block, 132 at E-CTS246);
  1 is unknown (SC_GMWOF5428755, coverage).
* **E-CTS246 sister lineages, subclade undetermined [obs-CRAM]:**
  * GGVP Mandinka SC_GMMAN5482181, 5482185, 5482191, 5482207, 5482299
  * GGVP Wolof SC_GMWOF5428789, 5428805
* **E-Z958, not E-CTS246 [obs-CRAM]:** 2 Fula (incl. YFull-placed 5306392, 5306400) and 4 Mandinka.
  2 more (Fula 5306398, Mandinka 5482307) are E-Z958 with CTS246 unresolved.
* **E-M132, not E-Z958:** 2 Mandinka (Z958 ancestral; single SNP).
* **Not E-M132:** 137 men (single defining SNP).
