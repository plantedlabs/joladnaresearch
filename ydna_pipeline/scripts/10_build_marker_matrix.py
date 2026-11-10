#!/usr/bin/env python3
"""
10_build_marker_matrix.py: sample-by-marker matrix from ALL existing observed evidence
(no new data access), plus an empirical nesting (topology) check.

Sources (observed, not reported):
  CRAM   results/jola_site_evidence_strict.csv           3 GGVP Jola men x 13 tree markers
  CRAM   results/ftf75935_screen/cache/*.json            203 GGVP males x FTF75935
  VCF    results/topology/nygc3202_tree_markers_raw.tsv  38 1000G E-Z958 men x 13 markers
         (NYGC GATK 3202 joint calls). Marker with no VCF record => 'nr'
         (no derived call in any of 3,202 samples; NOT verified ancestral per person).
Reported (not observed) placements are written separately (tree only).

Call codes: D derived | A ancestral | N no-call (coverage) | M mixed | X conflicting
            nr no VCF record (no derived call) | . not tested
Outputs:
  results/topology/marker_matrix_long.csv
  results/topology/marker_matrix_wide.csv
  results/topology/nesting_check.md
"""
import csv, glob, json, os, sys
from collections import defaultdict

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
M = [l.rstrip('\n').split('\t') for l in open('config/tree_markers.tsv') if l.strip() and not l.startswith('#')]
markers = [m[0] for m in M]; pos2m = {m[2]: m[0] for m in M}
ped = {r[1]: r for r in (l.split() for l in open('config/1kgp_3202_samples_ped_population.txt')) if r[0] != 'FamilyID'}
CODE = {'DERIVED_SUPPORTED': 'D', 'ANCESTRAL_SUPPORTED': 'A', 'INSUFFICIENT_COVERAGE': 'N',
        'UNRELIABLE_MAPPING': 'N', 'MIXED_AMBIGUOUS': 'M', 'NO_DATA_ERROR': 'N'}
long = []  # sample, population, dataset, marker, call, ref, alt, source


def add(s, pop, ds, mk, call, r='', a='', src=''):
    long.append(dict(sample_id=s, population_label=pop, dataset=ds, marker=mk, call=call, ref_count=r, alt_count=a, source=src))


# 1. GGVP Jola, CRAM strict
for r in csv.DictReader(open('results/jola_site_evidence_strict.csv')):
    if r['site_set'] != 'tree_marker':
        continue
    add(r['sample_id'], 'GGVP Jola', 'GGVP CRAM (low cov)', r['site_id'], CODE[r['call_classification']],
        r['filtered_ref_reads'], r['filtered_alt_reads'], 'observed: strict read evidence')
# 2. GGVP males, FTF75935
for f in glob.glob('results/ftf75935_screen/cache/*.json'):
    j = json.load(open(f))
    if j['inferred_sex'] != 'male' or j['sample_id'] in ('SC_GMJOL5309804', 'SC_GMJOL5309851', 'SC_GMJOL5309829'):
        continue
    add(j['sample_id'], f"GGVP {j['population'].title().replace('Woloff','Wolof')}", 'GGVP CRAM (low cov)', 'FTF75935',
        CODE.get(j['call_classification'], 'N'), j.get('ref_reads', ''), j.get('alt_reads', ''), 'observed: FTF75935 screen')
# 3. 1000G NYGC VCF
seen = defaultdict(set)
for l in open('results/topology/nygc3202_tree_markers_raw.tsv'):
    f = l.rstrip('\n').split('\t'); mk = pos2m[f[0]]
    for cell in f[4:]:
        s, rest = cell.split('='); gt, ad = rest.split(':')
        ads = ad.split(',') if ad != '.' else ['', '']
        call = 'D' if gt not in ('0', '.', './.') and '1' in gt else 'A' if gt in ('0', '0/0') else 'N'
        add(s, f"1000G {ped[s][5]}", 'NYGC 3202 chrY joint calls (30x)', mk, call, ads[0], ads[1] if len(ads) > 1 else '',
            'observed: VCF genotype'); seen[s].add(mk)
for s in list(seen):
    for mk in markers:
        if mk not in seen[s]:
            add(s, f"1000G {ped[s][5]}", 'NYGC 3202 chrY joint calls (30x)', mk, 'nr', '', '', 'no VCF record at this position')

with open('results/topology/marker_matrix_long.csv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, fieldnames=list(long[0])); w.writeheader(); w.writerows(long)
wide = defaultdict(dict); meta = {}
for r in long:
    wide[r['sample_id']][r['marker']] = r['call']; meta[r['sample_id']] = (r['population_label'], r['dataset'])
with open('results/topology/marker_matrix_wide.csv', 'w', newline='') as fo:
    w = csv.writer(fo); w.writerow(['sample_id', 'population_label', 'dataset'] + markers + ['deepest_supported_node'])
    order = ['M132', 'Z958', 'CTS246', 'Z15174', 'Z5991', 'FT402211', 'Z15196', 'FTF75935']

    def deepest(c):
        blk = [c.get(k, '.') for k in ('Z15174', 'Z15197', 'Z15231', 'Z15232', 'Z15234')]
        node = 'not determined'
        if c.get('M132') == 'D': node = 'E-M132'
        if c.get('Z958') == 'D': node = 'E-Z958'
        if c.get('CTS246') == 'D' or c.get('Z906') == 'D': node = 'E-CTS246'
        if blk.count('D') >= 2 and 'A' not in blk: node = 'E-Z15174'
        if node == 'E-Z15174' and c.get('Z5991') == 'D': node = 'E-Z5991'
        if node == 'E-Z5991' and c.get('FT402211') == 'D': node = 'E-FT402211'
        if node == 'E-Z5991' and c.get('Z15196') == 'D': node = 'E-Z15196'
        if c.get('FTF75935') == 'D': node += ' + FTF75935'
        if node.startswith('E-CTS246') and blk.count('A') + blk.count('nr') >= 2: node = 'E-CTS246 (not E-Z15174)'
        if node == 'E-Z958' and c.get('CTS246') in ('A', 'nr') and c.get('Z906') in ('A', 'nr'): node = 'E-Z958 (not E-CTS246)'
        return node
    for s in sorted(wide, key=lambda x: (meta[x][1], meta[x][0], x)):
        c = wide[s]; w.writerow([s, *meta[s]] + [c.get(k, '.') for k in markers] + [deepest(c)])

# Nesting check: for each marker pair, is carrier set of child a subset of parent's?
car = defaultdict(set); tested = defaultdict(set)
for r in long:
    tested[r['marker']].add(r['sample_id'])
    if r['call'] == 'D':
        car[r['marker']].add(r['sample_id'])
pairs = [('M132', 'Z958'), ('Z958', 'CTS246'), ('CTS246', 'Z906'), ('Z906', 'CTS246'), ('CTS246', 'Z15174'), ('Z15174', 'Z15197'),
         ('Z15174', 'Z15231'), ('Z15174', 'Z15232'), ('Z15174', 'Z15234'), ('Z15231', 'Z5991'), ('Z5991', 'FT402211'), ('Z5991', 'Z15196'),
         ('Z15231', 'FTF75935'), ('Z5991', 'FTF75935')]
L = ['| Parent | Child | child carriers | of which also parent-derived | parent-ancestral among child carriers | Verdict |', '|---|---|---|---|---|---|']
anc = defaultdict(set)
for r in long:
    if r['call'] in ('A', 'nr'):
        anc[r['marker']].add(r['sample_id'])
for p, c in pairs:
    cc = car[c]; both = cc & car[p]; bad = cc & anc[p]
    verdict = ('no carriers' if not cc else 'CONFLICT' if bad else 'nested (child ⊆ parent)' if both == cc else 'nested where tested')
    L.append(f'| {p} | {c} | {len(cc)} | {len(both)} | {len(bad)} | {verdict} |')
mut = [('FT402211', 'Z15196')]
for a_, b_ in mut:
    L.append(f'\nMutual exclusivity {a_} vs {b_}: carriers of both = {len(car[a_] & car[b_])}; '
             f'{a_}-only {len(car[a_] - car[b_])}, {b_}-only {len(car[b_] - car[a_])}')
open('results/topology/nesting_check.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
print(f"\nmatrix: {len(wide)} samples x {len(markers)} markers; long rows {len(long)}")
