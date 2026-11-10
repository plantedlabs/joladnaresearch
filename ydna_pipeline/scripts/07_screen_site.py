#!/usr/bin/env python3
"""
07_screen_site.py: cached, rate-limited population screen of ONE chrY site
(default FTF75935, chrY:2895562 A>C) across many remote CRAMs.

Per sample:
  1. Sex / Y-presence: count primary, non-dup, MAPQ>=20 reads in a 100 kb
     window that is 100% Umap-k100-unique (chrY:7,000,001-7,100,000).
     No base decoding is needed. Male >= MALE_MIN reads; female <= FEMALE_MAX; else ambiguous.
  2. Site evidence, using the same strict read filters and haploid classification
     as 02_cram_site_evidence.py (imported, not duplicated).

Results are cached per sample as JSON (re-runs skip finished samples). CRAM
indexes (.crai) are cached in --crai-dir so each index is downloaded once.
A delay between samples and a small worker pool limit load on the EBI FTP.

Usage:
  python3 07_screen_site.py --manifest config/ggvp_manifest.tsv --ref ref/chrY_analysisset.fa \
      --cache results/ftf75935_screen/cache [--samples S1,S2] [--limit N] [--workers 2]
"""
import argparse, json, os, subprocess, sys, time, statistics as st
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
_spec = importlib.util.spec_from_file_location('ev', os.path.join(os.path.dirname(os.path.abspath(__file__)), '02_cram_site_evidence.py'))
ev = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(ev)

SEX_WIN = 'chrY:7000001-7100000'
MALE_MIN, FEMALE_MAX = 300, 60


def run(cmd, cwd, timeout, retries=2):
    err = ''
    for i in range(retries + 1):
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd)
            if p.returncode == 0:
                return p.stdout, None
            err = p.stderr.strip()[-300:]
        except subprocess.TimeoutExpired:
            err = f'timeout {timeout}s'
        time.sleep(2 ** (i + 2))
    return None, err


def screen(sample, pop, url, a):
    out = os.path.join(a.cache, f'{sample}.json')
    if os.path.exists(out):
        return json.load(open(out)), True
    t0 = time.time()
    rec = dict(sample_id=sample, population=pop, source_dataset=a.dataset, cram_url=url,
               site_id=a.site_id, grch38_pos=a.pos, ref=a.ref_allele, alt=a.alt_allele, error='')
    o, err = run(['samtools', 'view', '-c', '-q', '20', '-F', '0xF04', '-T', a.ref, url, SEX_WIN], a.crai_dir, a.timeout)
    if o is None:
        rec.update(sex_window_reads='', inferred_sex='error', error=f'sex query: {err}')
    else:
        n = int(o.strip()); rec['sex_window_reads'] = n
        rec['inferred_sex'] = 'male' if n >= MALE_MIN else 'female' if n <= FEMALE_MAX else 'ambiguous'
    if rec['inferred_sex'] in ('male', 'ambiguous') or a.all_sexes:
        o, err = run(['samtools', 'view', '-T', a.ref, url, f'chrY:{a.pos}-{a.pos}'], a.crai_dir, a.timeout)
        if o is None:
            rec.update(call_classification='NO_DATA_ERROR', callable='no', error=f'site query: {err}')
        else:
            rec.update(evaluate(o.splitlines(), a))
    else:
        rec.update(call_classification='NOT_SCREENED_FEMALE', callable='n/a')
    rec['seconds'] = round(time.time() - t0, 1)
    json.dump(rec, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
    return rec, False


def evaluate(lines, a):
    ref_b, alt_b, pos = a.ref_allele, a.alt_allele, a.pos
    raw = excl = lowmq = lowbq = clipped = subopt = 0; pm = []; good = []; alt_drop = 0
    for ln in lines:
        f = ln.split('\t'); flag, mapq, start, cigar, seq, qual = int(f[1]), int(f[4]), int(f[3]), f[5], f[9], f[10]
        b, bq, dist, soft = ev.base_at(pos, start, cigar, seq, qual)
        if flag & ev.EXCL:
            excl += 1; continue
        raw += 1; pm.append(mapq)
        if b is None:
            continue
        aligned = sum(int(n) for n, o in ev.CIGAR_RE.findall(cigar) if o in 'M=X')
        axs = ev._as_xs(f)
        if mapq < a.min_mapq:
            lowmq += 1; st_ = 'lowmq'
        elif aligned < 0.5 * len(seq):
            clipped += 1; st_ = 'clip'
        elif axs is not None and axs < 0:
            subopt += 1; st_ = 'subopt'
        elif bq < a.min_bq:
            lowbq += 1; st_ = 'lowbq'
        else:
            good.append((b, mapq, '-' if flag & 16 else '+', bq)); st_ = 'used'
        if st_ != 'used' and b == alt_b:
            alt_drop += 1
    c = Counter(g[0] for g in good); r, al = c.get(ref_b, 0), c.get(alt_b, 0)
    oth = sum(v for k, v in c.items() if k not in (ref_b, alt_b)); d = len(good)
    flm = sum(1 for m in pm if m < a.min_mapq) / len(pm) if pm else 0.0
    cls = ev.classify(r, al, oth, flm, 'pass', d, lowmq)
    sup = al if cls == 'DERIVED_SUPPORTED' else r
    strands = len({g[2] for g in good if g[0] == (alt_b if cls == 'DERIVED_SUPPORTED' else ref_b)})
    conf = ev.confidence(cls, sup, strands, 'pass', flm)
    return dict(raw_primary_reads=raw, excluded_dup_secondary=excl, low_mapq_reads=lowmq, clipped_reads=clipped,
                suboptimal_reads=subopt, low_bq_reads=lowbq, usable_depth=d, ref_reads=r, alt_reads=al, other_reads=oth,
                alt_reads_in_filtered_out=alt_drop,
                alt_strands=f"{sum(1 for g in good if g[0]==alt_b and g[2]=='+')}/{sum(1 for g in good if g[0]==alt_b and g[2]=='-')}",
                mapq_median=(st.median(g[1] for g in good) if good else ''), bq_median=(st.median(g[3] for g in good) if good else ''),
                frac_lowmapq=round(flm, 3), call_classification=cls, confidence=conf,
                callable='yes' if cls in ('DERIVED_SUPPORTED', 'ANCESTRAL_SUPPORTED', 'MIXED_AMBIGUOUS') else 'no')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True); ap.add_argument('--ref', required=True)
    ap.add_argument('--cache', required=True); ap.add_argument('--crai-dir', default='ref/crai_cache')
    ap.add_argument('--dataset', default='GGVP (EBI 1000G FTP, GRCh38DH 2015 low-coverage CRAM)')
    ap.add_argument('--site-id', default='FTF75935'); ap.add_argument('--pos', type=int, default=2895562)
    ap.add_argument('--ref-allele', default='A'); ap.add_argument('--alt-allele', default='C')
    ap.add_argument('--samples', default=''); ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--workers', type=int, default=2); ap.add_argument('--delay', type=float, default=1.0)
    ap.add_argument('--timeout', type=int, default=180); ap.add_argument('--all-sexes', action='store_true')
    ap.add_argument('--min-mapq', type=int, default=20); ap.add_argument('--min-bq', type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.cache, exist_ok=True); os.makedirs(a.crai_dir, exist_ok=True)
    a.ref = os.path.abspath(a.ref)
    rows = [l.rstrip('\n').split('\t') for l in open(a.manifest) if l.strip() and not l.startswith('#')]
    if a.samples:
        want = a.samples.split(','); rows = [r for r in rows if r[0] in want]
    if a.limit:
        rows = rows[:a.limit]
    t0 = time.time(); done = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {}
        for r in rows:
            futs[ex.submit(screen, r[0], r[1], r[2], a)] = r[0]
            time.sleep(a.delay / max(1, a.workers))
        for fu in as_completed(futs):
            rec, cached = fu.result(); done += 1
            print(f"[{done}/{len(rows)}] {rec['sample_id']} {rec['population']} sexY={rec.get('sex_window_reads')} "
                  f"{rec['inferred_sex']} -> {rec.get('call_classification')} ref={rec.get('ref_reads','')} alt={rec.get('alt_reads','')} "
                  f"{'(cached)' if cached else str(rec.get('seconds'))+'s'} {rec.get('error','')}", file=sys.stderr, flush=True)
    print(f'elapsed {time.time()-t0:.0f}s', file=sys.stderr)


if __name__ == '__main__':
    main()
