#!/usr/bin/env python3
"""
05_rank_private_snps.py: rank the 19 private SNPs by phylogenetic reliability.

Transparent additive score (start at 10; every deduction is listed in the
output so the ranking can be audited and re-weighted):

  Evidence line                                   Deduction
  Big Y MAPQ < 30 (user-supplied)                 -3
  Umap k100 at site < 1.0 (published, genome-wide)-3   (<0.5: -4)
  Umap k50 min over +/-50 bp < 0.5                -1
  ENCODE blacklist region                         -3
  chrX/chrY 100-mer test: fail / caution(X-Y hit) -3 / -2  (caution w/o X/Y hit: -1)
  Segmental duplication fracMatch >= 0.98 / 0.95  -2 / -1
  Jola reads: mean fraction MAPQ<20 > 0.5 / > 0.25-3 / -1
  Jola strict-filter calls change vs standard     -1
  Another derived allele named at same position   -0.5 (mutable site; informational)
  Jola callable (ANC) in 3/3 men                  +1  (site genotypes cleanly in short reads)

Tiers: A >= 10, B 7-9.5, C 4-6.5, D < 4.
"""
import csv, sys, statistics as st
from collections import defaultdict

ann = {r['site_id']: r for r in csv.DictReader(open('results/site_annotation.csv'))}
gw = {r['site_id']: r for r in csv.DictReader(open('results/genomewide_mappability.csv'))}
std = defaultdict(dict); strict = defaultdict(dict)
for r in csv.DictReader(open('results/jola_site_evidence.csv')):
    std[r['site_id']][r['sample_id']] = r
for r in csv.DictReader(open('results/jola_site_evidence_strict.csv')):
    strict[r['site_id']][r['sample_id']] = r
yb = {}
for l in open('results/audit/ybrowse_position_lookup.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0].startswith('PV'):
        yb[f[0]] = f[3] if len(f) > 3 else ''
SAME_POS_OTHER_ALLELE = {'PV13': 'FT365779 (G>C) also named at this position'}
pv = [x for x in (l.rstrip('\n').split('\t') for l in open('config/private_variants.tsv')) if not x[0].startswith('#')]

rows = []
for sid, chrom, pos, ref, alt, mq30 in pv:
    sc, why = 10.0, []
    if mq30 == 'yes':
        sc -= 3; why.append('BigY MAPQ<30 (-3)')
    g = gw[sid]; k100 = float(g['umap_k100_site'])
    if k100 < 0.5:
        sc -= 4; why.append(f'Umap k100={k100} (-4)')
    elif k100 < 1:
        sc -= 3; why.append(f'Umap k100={k100} (-3)')
    if float(g['umap_k50_min_pm50']) < 0.5:
        sc -= 1; why.append(f"Umap k50 min±50={g['umap_k50_min_pm50']} (-1)")
    if g['encode_blacklist']:
        sc -= 3; why.append('ENCODE blacklist (-3)')
    a = ann[sid]
    if a['mappability_flag'] == 'fail':
        sc -= 3; why.append('chrX/Y 100-mer FAIL (-3)')
    elif a['mappability_flag'] == 'caution':
        d = 2 if int(a['tiles_with_alt_hits_nm_le2'] or 0) > 0 else 1
        sc -= d; why.append(f'chrX/Y 100-mer caution (-{d})')
    sd = float(a['segdup_best_fracMatch']) if a['segdup_best_fracMatch'] else 0
    if sd >= 0.98:
        sc -= 2; why.append(f'segdup {sd:.3f} (-2)')
    elif sd >= 0.95:
        sc -= 1; why.append(f'segdup {sd:.3f} (-1)')
    fl = [float(r['frac_primary_lowmapq'] or 0) for r in std[sid].values()]
    mfl = st.mean(fl) if fl else 0
    if mfl > 0.5:
        sc -= 3; why.append(f'Jola low-MAPQ frac {mfl:.2f} (-3)')
    elif mfl > 0.25:
        sc -= 1; why.append(f'Jola low-MAPQ frac {mfl:.2f} (-1)')
    calls_std = {k: v['call_classification'] for k, v in std[sid].items()}
    calls_str = {k: v['call_classification'] for k, v in strict[sid].items()}
    if calls_str and calls_std != calls_str:
        sc -= 1; why.append('strict filter changes a call (-1)')
    if sid in SAME_POS_OTHER_ALLELE:
        sc -= 0.5; why.append(SAME_POS_OTHER_ALLELE[sid] + ' (-0.5)')
    n_anc = sum(1 for v in calls_str.values() if v == 'ANCESTRAL_SUPPORTED')
    if n_anc == 3:
        sc += 1; why.append('ANC in 3/3 Jola (+1)')
    tier = 'A' if sc >= 10 else 'B' if sc >= 7 else 'C' if sc >= 4 else 'D'
    rows.append(dict(rank=0, site_id=sid, grch38_pos=pos, ref=ref, alt=alt, ybrowse_names=yb.get(sid, '').split('\t')[0] if yb.get(sid) else '',
                     score=sc, tier=tier, bigy_mapq_lt30=mq30, umap_k100=k100, umap_k50_min_pm50=g['umap_k50_min_pm50'],
                     umap_k36_site=g['umap_k36_site'], encode_blacklist=bool(g['encode_blacklist']),
                     chrXY_prelim=a['mappability_flag'], segdup_fracMatch=sd or '', repeatmasked=a['repeatmasked_site'],
                     jola_mean_lowmapq_frac=round(mfl, 2), jola_strict_calls='/'.join(calls_str.get(s, '?')[:3] for s in
                     ('SC_GMJOL5309804', 'SC_GMJOL5309851', 'SC_GMJOL5309829')), evidence=' ; '.join(why)))
rows.sort(key=lambda r: (-r['score'], int(r['grch38_pos'])))
for i, r in enumerate(rows, 1):
    r['rank'] = i
with open('results/private_snp_reliability_ranking.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['rank']:2d} {r['site_id']} {r['grch38_pos']:>9} {r['tier']} {r['score']:5.1f}  {r['evidence']}")
