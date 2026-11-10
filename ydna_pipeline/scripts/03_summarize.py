#!/usr/bin/env python3
"""03_summarize.py: Markdown evidence matrix from jola_site_evidence.csv + site_annotation.csv."""
import csv, sys
from collections import Counter, defaultdict

ev = list(csv.DictReader(open(sys.argv[1])))
ann = {r['site_id']: r for r in csv.DictReader(open(sys.argv[2]))}
SHORT = {'DERIVED_SUPPORTED': 'DER', 'ANCESTRAL_SUPPORTED': 'ANC', 'MIXED_AMBIGUOUS': 'MIX',
         'INSUFFICIENT_COVERAGE': 'LOWCOV', 'UNRELIABLE_MAPPING': 'BADMAP', 'NO_DATA_ERROR': 'ERR'}
samples = list(dict.fromkeys(r['sample_id'] for r in ev))
cell = {(r['sample_id'], r['site_id']): r for r in ev}

for title, prefix in (('Private variants of YF133256 (19 sites x 3 samples = 57 cells)', 'PV'),
                      ('E-Z15174-block branch markers (positive control)', 'Z')):
    sites = list(dict.fromkeys(r['site_id'] for r in ev if r['site_id'].startswith(prefix)))
    print(f'### {title}\n')
    print('Cell format: `CALL ref/alt/other (filtered depth) MAPQmed BQmed`. '
          'DER = derived allele supported; ANC = ancestral/reference supported; MIX = mixed/ambiguous; '
          'LOWCOV = insufficient coverage (<2 filtered reads); BADMAP = unreliable mapping.\n')
    print('| Site | GRCh38 | Ref>Der | Prelim. mappability | ' + ' | '.join(samples) + ' |')
    print('|' + '---|' * (4 + len(samples)))
    for s in sites:
        a = ann[s]; r0 = cell[(samples[0], s)]
        cells = []
        for smp in samples:
            r = cell.get((smp, s))
            if not r:
                cells.append('n/a'); continue
            c = SHORT[r['call_classification']]
            conf = f"({r['confidence']})" if r['confidence'] not in ('none', '') else ''
            cells.append(f"**{c}**{conf} {r['filtered_ref_reads'] or 0}/{r['filtered_alt_reads'] or 0}/"
                         f"{r['filtered_other_reads'] or 0} ({r['filtered_depth'] or 0}x) "
                         f"MQ{r['mapq_median_filtered'] or '-'} BQ{r['bq_median_filtered'] or '-'}")
        print(f"| {s} | {r0['grch38_pos']} | {r0['my_ref_allele']}>{r0['my_derived_allele']} | "
              f"{a['mappability_flag']} | " + ' | '.join(cells) + ' |')
    print()
    tally = defaultdict(Counter)
    for r in ev:
        if r['site_id'].startswith(prefix):
            tally[r['sample_id']][SHORT[r['call_classification']]] += 1
    print('| Sample | DER | ANC | MIX | LOWCOV | BADMAP | ERR |\n|---|---|---|---|---|---|---|')
    for smp in samples:
        t = tally[smp]
        print(f"| {smp} | " + ' | '.join(str(t.get(k, 0)) for k in ('DER', 'ANC', 'MIX', 'LOWCOV', 'BADMAP', 'ERR')) + ' |')
    print()
