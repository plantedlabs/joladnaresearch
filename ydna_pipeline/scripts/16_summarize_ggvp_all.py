#!/usr/bin/env python3
"""
16_summarize_ggvp_all.py: one classification for ALL 203 GGVP males (Fula, Jola, Mandinka, Wolof)
using identical node rules (imported from 12_branch_pilot.py).

Sources: pilot (10 test men), upstream screen (146 non-Jola), Jola screen (44), and the 3 Jola target men
from the strict 13-marker run (results/jola_site_evidence_strict.csv).
YFull placements (reconciliation table) are attached as REPORTED context, never as genotypes.

Outputs (results/ggvp_all/):
  ggvp_all_males_placements.csv     one row per man
  ggvp_all_sample_by_marker.csv     sample x 13 tree markers (D/A/N/M/.)
  ggvp_branch_distribution.md       counts by population; proportions in a separate, caveated table
"""
import csv, math, os, importlib.util
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(os.path.join(HERE, '..'))
spec = importlib.util.spec_from_file_location('bp', os.path.join(HERE, '12_branch_pilot.py'))
bp = importlib.util.module_from_spec(spec); spec.loader.exec_module(bp)
os.chdir(os.path.join(HERE, '..'))
TREE = ['M132', 'Z958', 'CTS246', 'Z906', 'Z15174', 'Z15197', 'Z15231', 'Z15232', 'Z15234', 'Z5991', 'FT402211', 'Z15196', 'FTF75935']
LAB = {'SC_GMFUL': 'GGVP Fula', 'SC_GMJOL': 'GGVP Jola', 'SC_GMMAN': 'GGVP Mandinka', 'SC_GMWOF': 'GGVP Wolof'}
calls = defaultdict(dict); run = {}
for path, tag in (('results/branch_pilot/pilot_marker_evidence_long.csv', 'pilot'),
                  ('results/upstream_screen/pilot_marker_evidence_long.csv', 'upstream screen'),
                  ('results/jola_screen/pilot_marker_evidence_long.csv', 'Jola screen')):
    for r in csv.DictReader(open(path)):
        s = r['sample_id']
        if not s.startswith('SC_GM') or r['marker'] not in TREE:
            continue
        if tag == 'pilot' and r['role'] != 'test':
            continue
        calls[s][r['marker']] = r['call']; run[s] = tag
CODE = {'DERIVED_SUPPORTED': 'D', 'ANCESTRAL_SUPPORTED': 'A', 'MIXED_AMBIGUOUS': 'M'}
for r in csv.DictReader(open('results/jola_site_evidence_strict.csv')):
    if r['site_set'] == 'tree_marker' and r['site_id'] in TREE:
        calls[r['sample_id']][r['site_id']] = CODE.get(r['call_classification'], 'N'); run[r['sample_id']] = 'Jola targets (strict 13-marker run)'
yf = {r['sample_id']: r['yfull_branch_context'] for r in csv.DictReader(open('results/audit/jola_male_reconciliation.csv')) if r['yfull_E-M132_page'] == 'True'}
for s in ('SC_GMFUL5306392', 'SC_GMFUL5306398', 'SC_GMFUL5306400'):
    yf[s] = 'E-Z5985 subtree'
yf.update({'SC_GMJOL5309804': 'E-FT402211', 'SC_GMJOL5309851': 'E-FT402211', 'SC_GMJOL5309829': 'E-Z15196'})


def cat(pl):
    for k, v in (('CONFLICT', 'conflicting'), ('E-Z15196', 'E-Z15174 > E-Z5991 > E-Z15196'), ('E-FT402211', 'E-Z15174 > E-Z5991 > E-FT402211'),
                 ('E-Z5991', 'E-Z15174 > E-Z5991'), ('E-Z15174', 'E-Z15174 (other)'), ('E-CTS246, not E-Z15174', 'E-CTS246, not E-Z15174'),
                 ('E-CTS246', 'E-CTS246 (E-Z15174 unresolved)'), ('E-Z958, not E-CTS246', 'E-Z958, not E-CTS246'),
                 ('E-Z958', 'E-Z958 (E-CTS246 unresolved)'), ('E-M132, not E-Z958', 'E-M132, not E-Z958'),
                 ('E-M132', 'E-M132 (E-Z958 unresolved)'), ('not E-M132', 'not E-M132'), ('unresolved', 'unresolved')):
        if pl.startswith(k):
            return v
    return 'unresolved'


rows = []
for s, c in calls.items():
    st = {n: bp.node_status(c, snps) for n, snps in bp.NODES}
    pl = bp.placement(st)
    z, ct = st['E-Z15174'][0], st['E-CTS246'][0]
    basis = ('DERIVED' if z == 'DERIVED' else 'excluded: >=2 block SNPs A' if z == 'NOT_DERIVED' else
             'excluded: E-CTS246 2/2 A' if ct == 'NOT_DERIVED' and '2/2' in st['E-CTS246'][1] else
             'weak: single-SNP exclusion' if ct == 'NOT_DERIVED' or st['E-M132'][0] == 'NOT_DERIVED' else 'UNKNOWN')
    rows.append(dict(sample_id=s, population_label=LAB[s[:8]], source=run[s], placement=pl, category=cat(pl),
                     **{n: st[n][0] for n, _ in bp.NODES}, E_Z15174_basis=basis, yfull_reported=yf.get(s, ''),
                     **{m: c.get(m, '.') for m in TREE}))
rows.sort(key=lambda r: (r['population_label'], r['sample_id']))
os.makedirs('results/ggvp_all', exist_ok=True)
keys = list(rows[0])
with open('results/ggvp_all/ggvp_all_males_placements.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=[k for k in keys if k not in TREE]); w.writeheader()
    w.writerows([{k: v for k, v in r.items() if k not in TREE} for r in rows])
with open('results/ggvp_all/ggvp_all_sample_by_marker.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['sample_id', 'population_label', 'source'] + TREE + ['placement', 'yfull_reported'])
    for r in rows:
        w.writerow([r['sample_id'], r['population_label'], r['source']] + [r[m] for m in TREE] + [r['placement'], r['yfull_reported']])

pops = ['GGVP Fula', 'GGVP Jola', 'GGVP Mandinka', 'GGVP Wolof']
cats = ['E-Z15174 > E-Z5991 > E-FT402211', 'E-Z15174 > E-Z5991 > E-Z15196', 'E-Z15174 > E-Z5991', 'E-Z15174 (other)',
        'E-CTS246, not E-Z15174', 'E-CTS246 (E-Z15174 unresolved)', 'E-Z958, not E-CTS246', 'E-Z958 (E-CTS246 unresolved)',
        'E-M132, not E-Z958', 'E-M132 (E-Z958 unresolved)', 'not E-M132', 'unresolved', 'conflicting']
cnt = Counter((r['population_label'], r['category']) for r in rows)
L = ['## Branch placements: SAMPLE COUNTS, all 203 GGVP males\n', '| Placement | ' + ' | '.join(p.replace('GGVP ', '') for p in pops) + ' | Total |',
     '|---|' + '---|' * (len(pops) + 1)]
for c in cats:
    t = sum(cnt[(p, c)] for p in pops)
    if t:
        L.append(f'| {c} | ' + ' | '.join(str(cnt[(p, c)]) for p in pops) + f' | {t} |')
L.append('| **Men** | ' + ' | '.join(str(sum(1 for r in rows if r['population_label'] == p)) for p in pops) + f' | {len(rows)} |')
b = Counter((r['population_label'], r['E_Z15174_basis']) for r in rows)
L += ['', '### E-Z15174 status basis\n', '| Basis | ' + ' | '.join(p.replace('GGVP ', '') for p in pops) + ' |', '|---|' + '---|' * len(pops)]
for k in ('DERIVED', 'excluded: >=2 block SNPs A', 'excluded: E-CTS246 2/2 A', 'weak: single-SNP exclusion', 'UNKNOWN'):
    L.append(f'| {k} | ' + ' | '.join(str(b[(p, k)]) for p in pops) + ' |')


def wilson(k, n):
    z = 1.96; p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0, c - h), min(1, c + h)


L += ['', '## Sample proportions (separate from counts; NOT population frequencies)\n',
      'Denominator = men whose status at the node is resolved (derived, or excluded by >=2 SNPs at the node or its parent). '
      'GGVP sampled families and specific Western Division communities, so relatives may share a Y and the Wilson '
      'intervals are too narrow; 10 non-Jola men were coverage-selected. Labels are sampling labels.\n',
      '| Node | ' + ' | '.join(p.replace('GGVP ', '') for p in pops) + ' |', '|---|' + '---|' * len(pops)]
for node in ('E-M132', 'E-CTS246', 'E-Z15174'):
    cells = []
    for p in pops:
        rs = [r for r in rows if r['population_label'] == p]
        if node == 'E-Z15174':
            d = sum(1 for r in rs if r['E_Z15174_basis'] == 'DERIVED'); n = d + sum(1 for r in rs if r['E_Z15174_basis'].startswith('excluded'))
        elif node == 'E-CTS246':
            d = sum(1 for r in rs if r['category'].startswith(('E-CTS246', 'E-Z15174')))
            n = d + sum(1 for r in rs if r['E-CTS246'] == 'NOT_DERIVED' and r['category'] != 'unresolved' or r['category'] == 'not E-M132')
        else:
            d = sum(1 for r in rs if not r['category'].startswith(('not E-M132', 'unresolved', 'conflicting'))); n = d + sum(1 for r in rs if r['category'] == 'not E-M132')
        pr, lo, hi = wilson(d, n)
        cells.append(f'{d}/{n} = {pr*100:.1f}% ({lo*100:.1f}–{hi*100:.1f})')
    L.append(f'| {node} | ' + ' | '.join(cells) + ' |')
conc = [(r['sample_id'], r['yfull_reported'], r['placement']) for r in rows if r['yfull_reported']]
L += ['', f'## YFull concordance ({len(conc)} GGVP men with a YFull placement)\n']
bad = [c for c in conc if ('E-Z5985' in c[1] or c[1].startswith(('E-FT4023', 'E-MF', 'E-Y31', 'E-Z1497', 'E-FT402336', 'E-FT402356'))) and c[2].startswith(('E-CTS246', 'E-Z15174', 'not E-M132'))]
L.append(f'Contradictions with YFull: {len(bad)}' + (': ' + '; '.join(f'{a} ({b} vs {c})' for a, b, c in bad) if bad else ''))
open('results/ggvp_all/ggvp_branch_distribution.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L)); print(len(rows), 'men')
