#!/usr/bin/env python3
"""
08_compile_screen.py: compile per-sample JSON cache from 07_screen_site.py into
  (1) an individual-level CSV (every sample, incl. females, marked not screened)
  (2) a population-level Markdown summary with callability and detection power.

Usage: 08_compile_screen.py CACHE_DIR OUT_CSV OUT_MD
"""
import csv, glob, json, sys, math
from collections import defaultdict, Counter

cache, out_csv, out_md = sys.argv[1:4]
recs = [json.load(open(f)) for f in sorted(glob.glob(f'{cache}/*.json'))]
cols = ['sample_id', 'population', 'source_dataset', 'inferred_sex', 'sex_window_reads', 'site_id', 'grch38_pos', 'ref', 'alt',
        'call_classification', 'confidence', 'callable', 'usable_depth', 'ref_reads', 'alt_reads', 'other_reads',
        'alt_strands', 'alt_reads_in_filtered_out', 'raw_primary_reads', 'excluded_dup_secondary', 'low_mapq_reads',
        'clipped_reads', 'suboptimal_reads', 'low_bq_reads', 'frac_lowmapq', 'mapq_median', 'bq_median',
        'confirmation_status', 'cram_url', 'error']
for r in recs:
    c = r.get('call_classification', '')
    r['confirmation_status'] = ('CANDIDATE: confirm with tree markers + Tier-A' if c == 'DERIVED_SUPPORTED'
                                else 'review (mixed evidence)' if c == 'MIXED_AMBIGUOUS'
                                else 'n/a')
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
    for r in sorted(recs, key=lambda x: (x['population'], x['sample_id'])):
        w.writerow({k: r.get(k, '') for k in cols})

pops = ['FULA', 'JOLA', 'MANDINKA', 'WOLOFF']
by = defaultdict(list)
for r in recs:
    by[r['population']].append(r)
L = []
L.append('| Population | Samples | Male | Female | Ambiguous | Error | Males screened | Callable (≥2 reads) | Ancestral | Derived | Mixed | Insufficient cov. | Median usable depth (males) |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
tot = Counter()
for p in pops + ['ALL']:
    rs = recs if p == 'ALL' else by[p]
    sx = Counter(r['inferred_sex'] for r in rs)
    m = [r for r in rs if r['inferred_sex'] == 'male']
    cl = Counter(r.get('call_classification') for r in m)
    callable_n = sum(1 for r in m if r.get('callable') == 'yes')
    deps = sorted(int(r.get('usable_depth') or 0) for r in m)
    med = deps[len(deps) // 2] if deps else '-'
    L.append(f"| {p} | {len(rs)} | {sx['male']} | {sx['female']} | {sx['ambiguous']} | {sx['error']} | {len(m)} | {callable_n} | "
             f"{cl['ANCESTRAL_SUPPORTED']} | {cl['DERIVED_SUPPORTED']} | {cl['MIXED_AMBIGUOUS']} | "
             f"{cl['INSUFFICIENT_COVERAGE'] + cl['UNRELIABLE_MAPPING']} | {med} |")
L.append('')
# depth distribution among males
L.append('Usable-depth distribution at FTF75935 (males):\n')
L.append('| Population | 0 | 1 | 2 | 3 | 4–5 | 6–9 | ≥10 |\n|---|---|---|---|---|---|---|---|')
bins = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 5), (6, 9), (10, 10**6)]
for p in pops + ['ALL']:
    m = [r for r in (recs if p == 'ALL' else by[p]) if r['inferred_sex'] == 'male']
    ds = [int(r.get('usable_depth') or 0) for r in m]
    L.append(f'| {p} | ' + ' | '.join(str(sum(1 for d in ds if lo <= d <= hi)) for lo, hi in bins) + ' |')
L.append('')
sexy = sorted(int(r['sex_window_reads']) for r in recs if r['inferred_sex'] in ('male', 'female', 'ambiguous'))
fem = [x for x in sexy if x <= 60]; mal = [x for x in sexy if x >= 300]
L.append(f"Sex-window read counts: females max {max(fem) if fem else '-'}, males min {min(mal) if mal else '-'} "
         f"(thresholds ≤60 female, ≥300 male).\n")
# alt reads anywhere
alt_any = [r for r in recs if int(r.get('alt_reads') or 0) > 0 or int(r.get('alt_reads_in_filtered_out') or 0) > 0]
L.append(f"Males with ANY read carrying C at FTF75935 (incl. filtered-out reads): {len(alt_any)}"
         + (': ' + ', '.join(f"{r['sample_id']} ({r['alt_reads']} usable, {r.get('alt_reads_in_filtered_out')} filtered)" for r in alt_any) if alt_any else ''))
open(out_md, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
