#!/usr/bin/env python3
"""
14_summarize_upstream.py: combine pilot (10 GGVP test men) + full upstream screen (146 men)
into one classification of all 156 GGVP non-Jola males, and summarise.

Outputs:
  results/upstream_screen/ggvp_nonjola_upstream_placements.csv   one row per man (156)
  results/upstream_screen/ggvp_nonjola_sample_by_marker.csv      sample x marker matrix (codes)
  results/upstream_screen/upstream_summary.md                    counts, proportions (separate), flags

Count tables are observations. Proportions are given separately, with Wilson 95% intervals, and carry the
caveats that GGVP sampled families (relatives share Y) and specific communities, so they are not
population frequencies.
"""
import csv, math, os
from collections import Counter, defaultdict
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

P = list(csv.DictReader(open('results/branch_pilot/pilot_placements.csv')))
U = list(csv.DictReader(open('results/upstream_screen/pilot_placements.csv')))
Mp = list(csv.DictReader(open('results/branch_pilot/pilot_sample_by_marker.csv')))
Mu = list(csv.DictReader(open('results/upstream_screen/pilot_sample_by_marker.csv')))
EVp = list(csv.DictReader(open('results/branch_pilot/pilot_marker_evidence_long.csv')))
EVu = list(csv.DictReader(open('results/upstream_screen/pilot_marker_evidence_long.csv')))
keep = lambda r: r['sample_id'].startswith(('SC_GMFUL', 'SC_GMMAN', 'SC_GMWOF'))
pilot_test = [r for r in P if r['role'] == 'test' and keep(r)]
rows = pilot_test + U
ids = [r['sample_id'] for r in rows]
assert len(ids) == len(set(ids)), 'duplicate sample'
mat = {r['sample_id']: r for r in Mp + Mu if keep(r)}
ev = defaultdict(dict)
for r in EVp + EVu:
    if keep(r):
        ev[r['sample_id']][r['marker']] = r


def category(pl):
    if pl.startswith('CONFLICT'): return 'conflicting'
    if 'candidate NEW' in pl or 'E-Z15174, not E-Z5991' in pl: return 'E-Z15174 (not E-Z5991): flag'
    if pl.startswith(('E-Z15196', 'E-FT402211', 'E-Z5991')): return 'E-Z15174 > E-Z5991'
    if pl.startswith('E-Z15174'): return 'E-Z15174 (unresolved below)'
    if pl.startswith('E-CTS246, not E-Z15174'): return 'E-CTS246, not E-Z15174'
    if pl.startswith('E-CTS246'): return 'E-CTS246 (E-Z15174 unresolved)'
    if pl.startswith('E-Z958, not E-CTS246'): return 'E-Z958, not E-CTS246'
    if pl.startswith('E-Z958'): return 'E-Z958 (E-CTS246 unresolved)'
    if pl.startswith('E-M132, not E-Z958'): return 'E-M132, not E-Z958'
    if pl.startswith('E-M132'): return 'E-M132 (E-Z958 unresolved)'
    if pl.startswith('not E-M132'): return 'not E-M132'
    if pl.startswith('not classified'): return 'not classified (not male)'
    return 'unresolved'


cats = ['E-Z15174 (not E-Z5991): flag', 'E-Z15174 > E-Z5991', 'E-Z15174 (unresolved below)', 'E-CTS246, not E-Z15174',
        'E-CTS246 (E-Z15174 unresolved)', 'E-Z958, not E-CTS246', 'E-Z958 (E-CTS246 unresolved)', 'E-M132, not E-Z958',
        'E-M132 (E-Z958 unresolved)', 'not E-M132', 'unresolved', 'conflicting', 'not classified (not male)']
pops = ['GGVP Fula', 'GGVP Mandinka', 'GGVP Wolof']
for r in rows:
    r['category'] = category(r['placement'])
    r['source_run'] = 'pilot (coverage-selected)' if r in pilot_test else ('screen (control reused)' if r['role'] == 'screen_control_reused' else 'full screen')
    e = ev[r['sample_id']]
    def sup(m):
        x = e.get(m)
        if not x: return ''
        s = x['alt_strands_fwd_rev'] if x['call'] == 'D' else x['ref_strands_fwd_rev']
        return f"{x['call']} {x['ref_count']}/{x['alt_count']}/{x['other_count']} str{s}" + (f" ({x['no_call_reason']})" if x['no_call_reason'] else '')
    for m in ('M132', 'Z958', 'CTS246', 'Z906'):
        r[f'{m}_evidence'] = sup(m)
    blk = [e.get(m, {}).get('call', '.') for m in ('Z15174', 'Z15197', 'Z15231', 'Z15232', 'Z15234')]
    r['Z15174_block_calls'] = ''.join(blk) if any(b != '.' for b in blk) else 'not tested'

cols = ['sample_id', 'population_label', 'source_run', 'sex_window_reads', 'category', 'placement', 'E-M132', 'E-Z958',
        'E-CTS246', 'E-Z15174', 'E-Z5991', 'E-FTF75935', 'M132_evidence', 'Z958_evidence', 'CTS246_evidence',
        'Z906_evidence', 'Z15174_block_calls', 'tierA_candidates']
with open('results/upstream_screen/ggvp_nonjola_upstream_placements.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
    for r in sorted(rows, key=lambda x: (x['population_label'], cats.index(x['category']), x['sample_id'])):
        w.writerow(r)
mcols = [c for c in Mu[0] if c not in ('role', 'dataset')]
with open('results/upstream_screen/ggvp_nonjola_sample_by_marker.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(mcols)
    for s in sorted(ids):
        w.writerow([mat[s].get(c, '.') for c in mcols])

L = ['## Observed branch placements: SAMPLE COUNTS (all 156 GGVP Fula, Mandinka and Wolof males)\n',
     '| Placement | ' + ' | '.join(pops) + ' | Total |', '|---|' + '---|' * (len(pops) + 1)]
cnt = Counter((r['population_label'], r['category']) for r in rows)
for c in cats:
    if sum(cnt[(p, c)] for p in pops):
        L.append(f'| {c} | ' + ' | '.join(str(cnt[(p, c)]) for p in pops) + f' | {sum(cnt[(p, c)] for p in pops)} |')
L.append('| **Total men** | ' + ' | '.join(str(sum(1 for r in rows if r['population_label'] == p)) for p in pops) + f' | {len(rows)} |')


def wilson(k, n):
    if n == 0: return (float('nan'),) * 3
    z = 1.96; p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0, c - h), min(1, c + h)


L += ['', '## Sample proportions (NOT population frequencies)\n',
      'Denominator = men with a **resolved** status at that node (unknowns excluded). Wilson 95% intervals assume '
      'independent, randomly sampled paternal lineages. GGVP sampled families and specific communities, so relatives '
      'may share a Y chromosome and the intervals are too narrow. 10 of these men (pilot) were coverage-selected.\n',
      '| Node | Population | Derived / resolved | Proportion (95% CI) |', '|---|---|---|---|']
for node, field in (('E-M132', 'E-M132'), ('E-CTS246', 'E-CTS246'), ('E-Z15174', 'E-Z15174')):
    for p in pops:
        rs = [r for r in rows if r['population_label'] == p]
        der = sum(1 for r in rs if r.get(field, '').startswith('DERIVED') or
                  (node == 'E-M132' and not r['category'].startswith(('not E-M132', 'unresolved', 'conflicting'))) or
                  (node == 'E-CTS246' and r['category'].startswith(('E-CTS246', 'E-Z15174'))) or
                  (node == 'E-Z15174' and r['category'].startswith('E-Z15174')))
        if node == 'E-M132':
            res = der + sum(1 for r in rs if r['category'] == 'not E-M132')
        elif node == 'E-CTS246':
            res = der + sum(1 for r in rs if r['category'] in ('not E-M132', 'E-Z958, not E-CTS246', 'E-M132, not E-Z958'))
        else:
            res = der + sum(1 for r in rs if r['category'] in ('not E-M132', 'E-Z958, not E-CTS246', 'E-M132, not E-Z958', 'E-CTS246, not E-Z15174'))
        pr, lo, hi = wilson(der, res)
        L.append(f'| {node} | {p} | {der} / {res} | {pr*100:.1f}% ({lo*100:.1f}–{hi*100:.1f}%) |' if res else f'| {node} | {p} | 0 / 0 | – |')
flags = [r for r in rows if r['category'] in ('E-Z15174 (not E-Z5991): flag', 'E-Z15174 > E-Z5991', 'E-Z15174 (unresolved below)',
                                            'E-CTS246 (E-Z15174 unresolved)', 'conflicting', 'unresolved', 'E-CTS246, not E-Z15174')]
L += ['', '## Flagged for follow-up\n', '| Sample | Population | Placement | Why flagged |', '|---|---|---|---|']
for r in flags:
    why = {'E-CTS246, not E-Z15174': 'E-CTS246 sister lineage; subclade undetermined (not on YFull)',
           'conflicting': 'conflicting calls: inspect reads', 'unresolved': 'insufficient coverage at defining sites'}.get(
        r['category'], 'E-Z15174-related: verify with reads + more SNPs')
    L.append(f"| {r['sample_id']} | {r['population_label']} | {r['placement']} | {why} |")
open('results/upstream_screen/upstream_summary.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
