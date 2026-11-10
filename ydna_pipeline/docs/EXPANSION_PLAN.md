# Prioritized plan for additional African samples

Principle: screen cheaply on **branch-defining markers first**, and run the 19 private sites only on
men who are independently shown to be E-Z15174 (or deeper). Haplogroup is not ethnicity. Population
labels below are the data provider's sampling labels.

## Update (v2): top priority
**Screen FTF75935 (chrY:2,895,562 A>C; your terminal SNP, well-mapped genome-wide) in all GGVP + 1000G GWD males.**
Any derived man is on your own branch. All three Jola men tested are ancestral. Confirm positives with
E-Z15174-block markers (Z15231, Z15232, Z15234) and tier-A private SNPs. Then do the IN147690 block-tree
comparison (REPORT §6).

Coordinate note: the repo's original "Z15234" position (8,414,659) is Z15197. Use `config/tree_markers.tsv`.

## Tier 0: before expanding (small, high value)
1. **Verify the branch-marker coordinates.** The four Z-SNP positions came from the repo script. The
   CRAM run (all three Jola men derived at all four) supports them empirically, but they should also be
   checked against YBrowse (`ybrowse.org`; that host would have to be allowed) or your Big Y VCF.
2. **Send your Big Y VCF + BED.** That would let us (a) confirm your genotype at the four markers,
   (b) take more E-Z15174-block SNPs (39 on YFull) and E-Z5991 / E-FT402211 / E-Z15196 SNPs
   (all ancestral in you) as additional, better-mapped screening markers, and (c) re-assess the 19
   sites with your read-level depth.
3. **Negative controls.** Run the 19 sites on 2–3 non-E-Z15174 Jola men (configured as
   `control` in `config/samples.tsv`; run with `ROLES=control`). A "private" site that turns up derived
   in unrelated men is not private: it would be a recurrent site, a reference artifact or a mapping artifact.

## Tier 1: 1000 Genomes GWD (Mandinka, Western Division), high coverage
* **HG03024 and HG03027** are already placed by YFull on **E-Z15196**, the sister clade of the
  E-FT402211 Jola pair. 1000G NYGC 30x CRAMs exist, about 7x the coverage of GGVP.
  Expected for the 19 sites: ancestral. A derived call here at a well-mapped site would contradict the tree.
* Other GWD men: screen on the E-Z15174 markers (~113 GWD samples, about half male).
* Cost: indexed queries only, ~4 sites × ~60 men × ~5 s ≈ 20 min, a few hundred MB transferred.
  Needs `ftp.1000genomes.ebi.ac.uk` (already allowed).

## Tier 2: GGVP Fula, Mandinka and Wolof (~4x low coverage)
* Available: 95 FULA, 99 MANDINKA and 100 WOLOFF CRAMs (from the alignment index).
* Step A: sex/coverage check, using `samtools view -c` in a 1 Mb X-degenerate chrY window (no
  reference decode, ~2 s per sample).
* Step B: males only, decode 4–8 E-Z15174 markers. Expect ~150 men × 6 sites × ~5 s ≈ 1.5 h.
* Step C: run the 19 private sites only on E-Z15174-positive men.
* Low coverage means many single-read sites. Require ≥2 concordant markers for a positive call.

## Tier 3: other public resources (eligibility to verify)
* **HGDP Mandenka (Senegal)**, high coverage. Check access route and whether any carry E-M132 / E-Z15174.
* **SGDP** West African samples.
* **Published E-M132 phylogenies**, for whether any Bainuk / Balanta / Manjak data exist.
  The repo README raises this, but I have not verified it.
* **FTDNA kit IN147690 (Gambian Fulani)**: not public. The SNP-level relationship can only be checked
  through the FTDNA Block Tree / group admins, or with the tester's consent to share BAM/VCF. If he is
  truly downstream of E-FTF75935, he should carry some of the 19 sites. That is the single most
  informative comparison available, and no public dataset substitutes for it.

## What would change the interpretation
* Any well-mapped (pass) private site **derived** in a man whom the tree places outside E-FTF75935
  → treat it as recurrent or artifact and remove it from the private set.
* A man who is **derived at several private sites** (≥3, well-mapped, ≥2 reads each) and ancestral
  for E-Z5991 → a candidate for a new shared branch with YF133256. That would be a genuine
  phylogenetic discovery requiring high-coverage confirmation.
