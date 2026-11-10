#!/usr/bin/env python3
"""06_tree_concordance.py: observed vs tree-expected genotypes at YBrowse-named markers.
Usage: 06_tree_concordance.py results/jola_site_evidence_strict.csv config/tree_markers.tsv > results/tree_concordance.md
"""
import csv, sys
ev = {(r['sample_id'], r['site_id']): r for r in csv.DictReader(open(sys.argv[1]))}
lines = [l.rstrip('\n').split('\t') for l in open(sys.argv[2]) if l.strip()]
hdr = lines[0][0].lstrip('#'), *lines[0][1:]
markers = [dict(zip(hdr, l)) for l in lines[1:]]
S = ['SC_GMJOL5309804', 'SC_GMJOL5309851', 'SC_GMJOL5309829']
CODE = {'DERIVED_SUPPORTED': 'D', 'ANCESTRAL_SUPPORTED': 'A'}
tally = {s: {'concordant': 0, 'discordant': 0, 'no_call': 0} for s in S}
print('| Marker | GRCh38 | Anc>Der | YFull node | ' + ' | '.join(f'{s[-7:]} exp / obs' for s in S) + ' |')
print('|' + '---|' * (4 + len(S)))
for m in markers:
    cells = []
    for s in S:
        r = ev.get((s, m['id']))
        exp = m[f'expect_{s[-7:]}']
        if r is None:
            cells.append(f'{exp} / n/a'); continue
        obs = CODE.get(r['call_classification'])
        detail = f"{r['filtered_ref_reads']}/{r['filtered_alt_reads']}/{r['filtered_other_reads']}"
        if obs is None:
            tally[s]['no_call'] += 1
            short = {'INSUFFICIENT_COVERAGE': 'no-call (cov)', 'UNRELIABLE_MAPPING': 'no-call (map)',
                     'MIXED_AMBIGUOUS': 'MIXED'}.get(r['call_classification'], r['call_classification'])
            cells.append(f'{exp} / {short} {detail}')
        elif obs == exp:
            tally[s]['concordant'] += 1; cells.append(f'{exp} / **{obs}** ✔ {detail}')
        else:
            tally[s]['discordant'] += 1; cells.append(f'{exp} / **{obs}** ✘ {detail}')
    print(f"| {m['id']} | {m['pos_grch38']} | {m['ref']}>{m['alt']} | {m['yfull_node']} | " + ' | '.join(cells) + ' |')
print('\nCell: expected / observed, then filtered ref/alt/other reads. D = derived, A = ancestral.\n')
print('| Sample | Concordant | Discordant | No-call |\n|---|---|---|---|')
for s in S:
    t = tally[s]; print(f"| {s} | {t['concordant']} | {t['discordant']} | {t['no_call']} |")
