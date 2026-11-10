#!/usr/bin/env python3
"""
12_branch_pilot.py: hierarchical upstream/sister/intermediate branch classification for a
small, fixed set of men (pilot). Stops after the listed samples; never expands automatically.

Markers come from config/tree_markers.tsv and are verified beforehand by 11_verify_markers.py.
Tier-A private SNPs (config/tierA_private_variants.tsv) are CANDIDATE markers: they are reported
as shared/not shared but never used to define a branch.

Read evaluation: 07_screen_site.evaluate (strict filters; depth = DNA fragments; >=2 fragments and
>=90% agreement for a call). Each sample x marker query is cached as JSON.

Node rules (missing data stays unknown):
  node DERIVED      >=1 defining SNP D, 0 A           (confidence 'multi-SNP' if >=2 D)
  node NOT DERIVED  >=2 defining SNPs A, 0 D           (single-SNP nodes: 1 A => 'not derived, 1 SNP only')
  node CONFLICT     D and A both present among defining SNPs
  node UNKNOWN      otherwise (no-calls)
  Nodes and defining SNPs: E-M132 {M132}; E-Z958 {Z958}; E-CTS246 {CTS246, Z906};
  E-Z15174 {Z15174, Z15197, Z15231, Z15232, Z15234}; E-Z5991 {Z5991}; E-FT402211 {FT402211};
  E-Z15196 {Z15196}; E-FTF75935 {FTF75935}.
Gating: step 1 (M132, Z958, CTS246, Z906) for every male; step 2 (E-Z15174 block) if E-CTS246 is
DERIVED/UNKNOWN/CONFLICT; step 3 (Z5991, FT402211, Z15196, FTF75935 + Tier-A candidates) if E-Z15174 is
DERIVED/UNKNOWN/CONFLICT. Controls (--force-full) run all steps.
"""
import argparse, csv, json, os, sys, time, types
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(os.path.join(HERE, '..'))
import importlib.util
spec = importlib.util.spec_from_file_location('scr', os.path.join(HERE, '07_screen_site.py'))
scr = importlib.util.module_from_spec(spec); spec.loader.exec_module(scr)

NODES = [('E-M132', ['M132']), ('E-Z958', ['Z958']), ('E-CTS246', ['CTS246', 'Z906']),
         ('E-Z15174', ['Z15174', 'Z15197', 'Z15231', 'Z15232', 'Z15234']), ('E-Z5991', ['Z5991']),
         ('E-FT402211', ['FT402211']), ('E-Z15196', ['Z15196']), ('E-FTF75935', ['FTF75935'])]
STEP1 = ['M132', 'Z958', 'CTS246', 'Z906']
STEP2 = ['Z15174', 'Z15197', 'Z15231', 'Z15232', 'Z15234']
STEP3 = ['Z5991', 'FT402211', 'Z15196', 'FTF75935']
CODE = {'DERIVED_SUPPORTED': 'D', 'ANCESTRAL_SUPPORTED': 'A', 'MIXED_AMBIGUOUS': 'M'}


def rows(p):
    L = [l.rstrip('\n').split('\t') for l in open(p) if l.strip()]
    h = L[0][0].lstrip('#'), *L[0][1:]
    return [dict(zip(h, x)) for x in L[1:] if not x[0].startswith('#')]


def no_call_reason(e):
    if e.get('call_classification') == 'NO_DATA_ERROR':
        return 'query error: ' + e.get('error', '')
    d = int(e.get('usable_depth') or 0)
    if e['call_classification'] == 'UNRELIABLE_MAPPING':
        return f"majority of reads MAPQ<20 ({e.get('frac_lowmapq')})"
    if d == 0:
        raw = int(e.get('raw_primary_reads') or 0)
        return 'zero coverage' if raw == 0 else (f"{raw} primary read(s), none usable (lowMQ {e.get('low_mapq_reads')}, "
                                                 f"clipped {e.get('clipped_reads')}, XS>AS {e.get('suboptimal_reads')}, lowBQ {e.get('low_bq_reads')})")
    if d == 1:
        return 'single usable fragment'
    if e['call_classification'] == 'MIXED_AMBIGUOUS':
        return f"mixed alleles (ref {e['ref_reads']}, alt {e['alt_reads']}, other {e['other_reads']})"
    return ''


def query(sample, url, mk, a):
    path = os.path.join(a.cache, f"{sample}__{mk['id']}.json")
    if os.path.exists(path):
        return json.load(open(path))
    ns = types.SimpleNamespace(ref_allele=mk['ref'], alt_allele=mk['alt'], pos=int(mk['pos_grch38']), min_mapq=20, min_bq=20)
    o, err = scr.run(['samtools', 'view', '-T', a.ref, url, f"chrY:{mk['pos_grch38']}-{mk['pos_grch38']}"], a.crai_dir, a.timeout)
    e = {'call_classification': 'NO_DATA_ERROR', 'error': err} if o is None else scr.evaluate(o.splitlines(), ns)
    e.update(sample_id=sample, marker=mk['id'], grch38_pos=mk['pos_grch38'], anc=mk['ref'], der=mk['alt'])
    e['call'] = CODE.get(e['call_classification'], 'N'); e['no_call_reason'] = no_call_reason(e) if e['call'] in ('N', 'M') else ''
    json.dump(e, open(path, 'w')); time.sleep(a.delay)
    return e


def node_status(calls, snps):
    c = [calls.get(s, '.') for s in snps]; d, an = c.count('D'), c.count('A')
    if d and an:
        return 'CONFLICT', f'{d}D/{an}A'
    if d:
        return 'DERIVED', f"{d}/{len(snps)} SNP D" + (' (multi-SNP)' if d >= 2 else ' (single SNP)')
    if an >= 2 or (len(snps) == 1 and an == 1):
        return 'NOT_DERIVED', f"{an}/{len(snps)} SNP A" + (' (single SNP only)' if an == 1 else '')
    tested = sum(1 for x in c if x != '.')
    return ('UNKNOWN', f'{tested} tested, {an} A, rest no-call') if tested else ('NOT_TESTED', '')


def placement(st):
    order = ['E-M132', 'E-Z958', 'E-CTS246', 'E-Z15174', 'E-Z5991']
    deepest = 'unplaced (E-M132 not derived)' if st['E-M132'][0] == 'NOT_DERIVED' else 'unresolved'
    for n in order:
        if st[n][0] == 'DERIVED':
            deepest = n
        elif st[n][0] in ('UNKNOWN', 'CONFLICT') and deepest.startswith('E-'):
            return deepest + f' ({n} {st[n][0].lower()})'
        elif st[n][0] == 'NOT_DERIVED' and deepest.startswith('E-'):
            label = {'E-Z958': 'other E-M132 (not E-Z958)', 'E-CTS246': 'E-Z958, not E-CTS246',
                     'E-Z15174': 'E-CTS246, not E-Z15174 (sister lineage)',
                     'E-Z5991': 'E-Z15174, not E-Z5991'}.get(n, deepest)
            if n == 'E-Z5991':
                f = st['E-FTF75935'][0]
                return label + (' and not E-FTF75935: candidate NEW E-Z15174 lineage' if f == 'NOT_DERIVED' else
                                ' and E-FTF75935 DERIVED' if f == 'DERIVED' else f' (E-FTF75935 {f.lower()})')
            return label
    if deepest == 'E-Z5991':
        for sub in ('E-FT402211', 'E-Z15196'):
            if st[sub][0] == 'DERIVED':
                return sub
    return deepest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', required=True, help='TSV: sample, population, dataset, role, cram_url, selection_note')
    ap.add_argument('--ref', default='ref/chrY_analysisset.fa'); ap.add_argument('--crai-dir', default='ref/crai_cache')
    ap.add_argument('--cache', default='results/branch_pilot/cache'); ap.add_argument('--outdir', default='results/branch_pilot')
    ap.add_argument('--timeout', type=int, default=240); ap.add_argument('--delay', type=float, default=0.5)
    a = ap.parse_args(); a.ref = os.path.abspath(a.ref)
    os.makedirs(a.cache, exist_ok=True); os.makedirs(a.crai_dir, exist_ok=True)
    mk = {m['id']: m for m in rows('config/tree_markers.tsv')}
    tier = {m['id']: m for m in rows('config/tierA_private_variants.tsv')}
    allm = {**mk, **tier}
    S = rows(a.samples)
    long, wide, place = [], [], []
    for s in S:
        sid, url = s['sample'], s['cram_url']; force = s['role'].startswith('control')
        sexp = os.path.join(a.cache, f'{sid}__SEX.json')
        if os.path.exists(sexp):
            sex = json.load(open(sexp))
        else:
            o, err = scr.run(['samtools', 'view', '-c', '-q', '20', '-F', '0xF04', '-T', a.ref, url, scr.SEX_WIN], a.crai_dir, a.timeout)
            n = int(o) if o else -1
            sex = dict(sex_window_reads=n, inferred_sex='male' if n >= scr.MALE_MIN else 'female' if 0 <= n <= scr.FEMALE_MAX else 'ambiguous' if n > 0 else 'error')
            json.dump(sex, open(sexp, 'w'))
        print(f"== {sid} ({s['role']}) sexY={sex['sex_window_reads']} {sex['inferred_sex']}", file=sys.stderr, flush=True)
        if sex['inferred_sex'] != 'male':
            place.append(dict(sample_id=sid, population_label=s['population'], dataset=s['dataset'], role=s['role'],
                              selection_note=s['selection_note'], sex_window_reads=sex['sex_window_reads'],
                              inferred_sex=sex['inferred_sex'], placement='not classified (not male)'))
            continue
        calls, ev = {}, {}
        todo = list(STEP1)
        for step in (1, 2, 3):
            for m in todo:
                e = query(sid, url, allm[m], a); calls[m] = e['call']; ev[m] = e
                print(f"   {m:9s} {e['call']} ref={e.get('ref_reads','')} alt={e.get('alt_reads','')} "
                      f"oth={e.get('other_reads','')} d={e.get('usable_depth','')} {e.get('no_call_reason','')}", file=sys.stderr, flush=True)
            st = {n: node_status(calls, snps) for n, snps in NODES}
            if step == 1:
                todo = STEP2 if (force or st['E-CTS246'][0] in ('DERIVED', 'UNKNOWN', 'CONFLICT')) else []
            elif step == 2:
                todo = (STEP3 + list(tier)) if (force or st['E-Z15174'][0] in ('DERIVED', 'UNKNOWN', 'CONFLICT')) else []
            if not todo:
                break
        for m in STEP1 + STEP2 + STEP3 + list(tier):
            e = ev.get(m)
            if e is None:
                continue
            long.append(dict(sample_id=sid, population_label=s['population'], dataset=s['dataset'], role=s['role'],
                             marker=m, marker_kind='tierA_candidate' if m in tier else 'tree_marker', grch38_pos=e['grch38_pos'],
                             anc_der=f"{e['anc']}>{e['der']}", call=e['call'], usable_fragments=e.get('usable_depth', ''),
                             ref_count=e.get('ref_reads', ''), alt_count=e.get('alt_reads', ''), other_count=e.get('other_reads', ''),
                             ref_strands_fwd_rev=e.get('ref_strands', ''), alt_strands_fwd_rev=e.get('alt_strands', ''),
                             mapq_median=e.get('mapq_median', ''), bq_median=e.get('bq_median', ''), confidence=e.get('confidence', ''),
                             alt_in_filtered_out=e.get('alt_reads_in_filtered_out', ''), no_call_reason=e.get('no_call_reason', '')))
        st = {n: node_status(calls, snps) for n, snps in NODES}
        tier_c = Counter(calls.get(t, '.') for t in tier)
        place.append(dict(sample_id=sid, population_label=s['population'], dataset=s['dataset'], role=s['role'],
                          selection_note=s['selection_note'], sex_window_reads=sex['sex_window_reads'], inferred_sex='male',
                          **{n: f"{st[n][0]} [{st[n][1]}]" if st[n][1] else st[n][0] for n, _ in NODES},
                          tierA_candidates=(f"tested {sum(tier_c.values()) - tier_c['.']}: D={tier_c['D']} A={tier_c['A']} "
                                            f"no-call={tier_c['N'] + tier_c['M']}" if any(t in calls for t in tier) else 'not tested'),
                          placement=placement(st)))
        wide.append(dict(sample_id=sid, population_label=s['population'], dataset=s['dataset'], role=s['role'],
                         **{m: calls.get(m, '.') for m in STEP1 + STEP2 + STEP3 + list(tier)}, placement=place[-1]['placement']))
    os.makedirs(a.outdir, exist_ok=True)
    for name, data in (('pilot_marker_evidence_long.csv', long), ('pilot_sample_by_marker.csv', wide), ('pilot_placements.csv', place)):
        if data:
            keys = list(dict.fromkeys(k for d in data for k in d))
            with open(os.path.join(a.outdir, name), 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(data)
    for p in place:
        print(f"{p['sample_id']:18s} {p['role']:22s} {p['placement']}", file=sys.stderr)


if __name__ == '__main__':
    main()
