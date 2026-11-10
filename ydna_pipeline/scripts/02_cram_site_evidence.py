#!/usr/bin/env python3
"""
02_cram_site_evidence.py: per-read evidence at target chrY sites in remote CRAMs.

Requires only python3 (stdlib) + samtools >= 1.10 on PATH, and the MD5-verified
local chrY from 00_setup_reference.sh (passed via -T, so htslib never contacts
the slow ENA reference server).

For each sample x site it streams `samtools view` over HTTPS (indexed byte-range
access, ~1 CRAM container per query), parses CIGAR to find the base at the target
position, and applies read-level filters:

  excluded:  unmapped(0x4), secondary(0x100), QC-fail(0x200), duplicate(0x400),
             supplementary(0x800), deletion/ref-skip at the site
  low MAPQ:  MAPQ < --min-mapq (default 20)  -> counted, not used for calling
  low BQ:    base quality < --min-bq (default 20) -> counted, not used
  near-end:  base within 5 bp of a read end (incl. soft clip) -> counted, still used

Classification (haploid chrY; zero/low coverage is NEVER called reference):
  UNRELIABLE_MAPPING     >50% (and >=2) of primary non-dup reads have MAPQ < min-mapq,
                         or preliminary mappability flag = 'fail'
  INSUFFICIENT_COVERAGE  < 2 filtered reads
  DERIVED_SUPPORTED      >= 2 filtered alt reads and alt fraction >= 0.9
  ANCESTRAL_SUPPORTED    >= 2 filtered ref reads and ref fraction >= 0.9
  MIXED_AMBIGUOUS        anything else (incl. third alleles)
Confidence: high = >=4 supporting reads on both strands, no mappability
caveat and <=25% low-MAPQ reads; moderate = >=3 supporting reads; low = 2 supporting reads or caveat.
"""
import argparse, csv, re, subprocess, sys, time, statistics as st
from collections import Counter

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
EXCL = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def read_tsv(path):
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    hdr = lines[0].lstrip('#').split('\t')
    return [dict(zip(hdr, l.split('\t'))) for l in lines[1:] if not l.startswith('#')]


def base_at(pos1, rstart1, cigar, seq, qual):
    """Return (base, bq, dist_to_read_end, n_softclip) or (None,..) if deleted/not covered."""
    rpos, qpos = rstart1, 0
    ops = [(int(n), o) for n, o in CIGAR_RE.findall(cigar)]
    qlen_total = sum(n for n, o in ops if o in 'MIS=X')
    nsoft = sum(n for n, o in ops if o == 'S')
    for n, o in ops:
        if o in 'M=X':
            if rpos <= pos1 < rpos + n:
                q = qpos + (pos1 - rpos)
                bq = ord(qual[q]) - 33 if qual != '*' else -1
                return seq[q].upper(), bq, min(q, qlen_total - 1 - q), nsoft
            rpos += n; qpos += n
        elif o in 'IS':
            qpos += n
        elif o in 'DN':
            if rpos <= pos1 < rpos + n:
                return None, None, None, nsoft
            rpos += n
    return None, None, None, nsoft


def fetch(cram, ref_fa, region, timeout, retries):
    cmd = ['samtools', 'view', '-T', ref_fa, cram, region]
    for attempt in range(retries + 1):
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            if p.returncode == 0:
                return p.stdout.splitlines(), None
            err = p.stderr.strip()[-300:]
        except subprocess.TimeoutExpired:
            err = f'timeout after {timeout}s'
        time.sleep(2 ** (attempt + 1))
    return None, err


def _as_xs(fields):
    t = {x[:2]: x[5:] for x in fields[11:]}
    try:
        return int(t['AS']) - int(t['XS'])
    except (KeyError, ValueError):
        return None


def classify(ref_n, alt_n, oth_n, frac_lowmq, mapflag, depth, n_lowmq):
    if mapflag == 'fail' or (frac_lowmq > 0.5 and n_lowmq >= 2):
        return 'UNRELIABLE_MAPPING'
    if depth < 2:
        return 'INSUFFICIENT_COVERAGE'
    if alt_n >= 2 and alt_n / depth >= 0.9:
        return 'DERIVED_SUPPORTED'
    if ref_n >= 2 and ref_n / depth >= 0.9:
        return 'ANCESTRAL_SUPPORTED'
    return 'MIXED_AMBIGUOUS'


def confidence(cls, sup_n, sup_strands, mapflag, frac_lowmq=0.0):
    if cls not in ('DERIVED_SUPPORTED', 'ANCESTRAL_SUPPORTED'):
        return 'none'
    if sup_n >= 4 and sup_strands == 2 and mapflag == 'pass' and frac_lowmq <= 0.25:
        return 'high'
    if sup_n >= 3 and mapflag != 'caution':
        return 'moderate'
    return 'low'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', required=True)
    ap.add_argument('--sites', required=True, nargs='+', help='one or more site TSVs')
    ap.add_argument('--ref', required=True, help='MD5-verified chrY fasta')
    ap.add_argument('--mappability', help='CSV from 01_annotate_sites.py (site_id,mappability_flag)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--reads-out', required=True, help='per-read evidence TSV')
    ap.add_argument('--min-mapq', type=int, default=20)
    ap.add_argument('--min-bq', type=int, default=20)
    ap.add_argument('--timeout', type=int, default=180)
    ap.add_argument('--retries', type=int, default=2)
    ap.add_argument('--strict', action='store_true',
                    help='also exclude reads with <50%% of bases aligned (heavy clipping) or XS > AS (better alignment elsewhere)')
    ap.add_argument('--roles', default='', help='comma list of sample roles to include (default all)')
    a = ap.parse_args()

    samples = read_tsv(a.samples)
    if a.roles:
        keep = set(a.roles.split(','))
        samples = [s for s in samples if s['role'] in keep]
    sites = []
    for path in a.sites:
        for s in read_tsv(path):
            s['site_set'] = 'private' if s['id'].startswith('PV') else 'tree_marker'
            sites.append(s)
    mapflags = {}
    if a.mappability:
        with open(a.mappability) as f:
            for r in csv.DictReader(f):
                mapflags[r['site_id']] = (r['mappability_flag'], r.get('mappability_notes', ''))

    # reference base check
    refseq = {}
    for s in sites:
        r = subprocess.run(['samtools', 'faidx', a.ref, f"{s['chrom']}:{s['pos_grch38']}-{s['pos_grch38']}"],
                           capture_output=True, text=True, check=True).stdout.split('\n')[1].upper()
        refseq[s['id']] = r

    cols = ['sample_id', 'population', 'ydna_branch_yfull', 'filter_mode', 'site_set', 'site_id', 'grch38_pos',
            'my_ref_allele', 'my_derived_allele', 'reference_base_check',
            'raw_primary_reads', 'dup_or_secondary_excluded', 'low_mapq_reads', 'low_bq_reads',
            'filtered_depth', 'filtered_ref_reads', 'filtered_alt_reads', 'filtered_other_reads',
            'observed_alleles', 'alt_strands_fwd_rev', 'ref_strands_fwd_rev',
            'mapq_median_filtered', 'mapq_min_filtered', 'frac_primary_lowmapq',
            'bq_median_filtered', 'alt_reads_in_filtered_out', 'near_read_end_reads', 'softclipped_reads',
            'mappability_prelim', 'callable', 'call_classification', 'confidence', 'interpretation', 'error']
    out = open(a.out, 'w', newline=''); w = csv.DictWriter(out, fieldnames=cols); w.writeheader()
    rout = open(a.reads_out, 'w')
    rout.write('sample\tsite\tread\tflag\tmapq\tstrand\tbase\tbq\tdist_end\tsoftclip\tstatus\n')

    for smp in samples:
        for s in sites:
            pos = int(s['pos_grch38']); region = f"{s['chrom']}:{pos}-{pos}"
            t0 = time.time()
            lines, err = fetch(smp['cram_url'], a.ref, region, a.timeout, a.retries)
            mf, mnote = mapflags.get(s['id'], ('not_assessed', ''))
            row = dict(sample_id=smp['sample'], population=smp['population'],
                       ydna_branch_yfull=smp['yfull_branch'], filter_mode=('strict' if a.strict else 'standard'), site_set=s['site_set'], site_id=s['id'],
                       grch38_pos=pos, my_ref_allele=s['ref'], my_derived_allele=s['alt'],
                       reference_base_check=('match' if refseq[s['id']] == s['ref'].upper()
                                             else f"MISMATCH(ref={refseq[s['id']]})"),
                       mappability_prelim=mf, error='')
            if lines is None:
                row.update(call_classification='NO_DATA_ERROR', callable='no', confidence='none',
                           interpretation='CRAM query failed; no evidence', error=err)
                w.writerow(row); out.flush(); print(f"{smp['sample']} {s['id']} ERROR {err}", file=sys.stderr)
                continue
            raw = excl = lowmq = lowbq = 0
            primary_mapqs = []
            good = []  # (base, mapq, strand, bq, dist, soft)
            alt_dropped = 0
            for ln in lines:
                f = ln.split('\t')
                flag, mapq, start, cigar, seq, qual = int(f[1]), int(f[4]), int(f[3]), f[5], f[9], f[10]
                b, bq, dist, soft = base_at(pos, start, cigar, seq, qual)
                strand = '-' if flag & 0x10 else '+'
                if flag & EXCL:
                    excl += 1; status = 'excluded_flag'
                elif b is None:
                    raw += 1; primary_mapqs.append(mapq); status = 'deletion_or_skip'
                else:
                    raw += 1; primary_mapqs.append(mapq)
                    if mapq < a.min_mapq:
                        lowmq += 1; status = 'low_mapq'
                    elif a.strict and (sum(int(n) for n, o in CIGAR_RE.findall(cigar) if o in 'M=X') < 0.5 * len(seq)):
                        status = 'clipped'
                    elif a.strict and _as_xs(f) is not None and _as_xs(f) < 0:
                        status = 'suboptimal'
                    elif bq < a.min_bq:
                        lowbq += 1; status = 'low_bq'
                    else:
                        good.append((b, mapq, strand, bq, dist, soft)); status = 'used'
                if status in ('low_mapq', 'low_bq', 'clipped', 'suboptimal') and b == s['alt'].upper():
                    alt_dropped += 1
                rout.write(f"{smp['sample']}\t{s['id']}\t{f[0]}\t{flag}\t{mapq}\t{strand}\t{b}\t{bq}\t{dist}\t{soft}\t{status}\n")
            ref_b, alt_b = s['ref'].upper(), s['alt'].upper()
            cnt = Counter(g[0] for g in good)
            ref_n, alt_n = cnt.get(ref_b, 0), cnt.get(alt_b, 0)
            oth_n = sum(v for k, v in cnt.items() if k not in (ref_b, alt_b))
            depth = len(good)
            frac_lowmq = (sum(1 for m in primary_mapqs if m < a.min_mapq) / len(primary_mapqs)) if primary_mapqs else 0.0
            strands = lambda base: (sum(1 for g in good if g[0] == base and g[2] == '+'),
                                    sum(1 for g in good if g[0] == base and g[2] == '-'))
            af, ar = strands(alt_b); rf, rr = strands(ref_b)
            cls = classify(ref_n, alt_n, oth_n, frac_lowmq, mf, depth, lowmq)
            if cls == 'DERIVED_SUPPORTED':
                conf = confidence(cls, alt_n, (af > 0) + (ar > 0), mf, frac_lowmq)
            elif cls == 'ANCESTRAL_SUPPORTED':
                conf = confidence(cls, ref_n, (rf > 0) + (rr > 0), mf, frac_lowmq)
            else:
                conf = 'none'
            interp = {
                'DERIVED_SUPPORTED': f"carries {alt_b} (same allele as YF133256)",
                'ANCESTRAL_SUPPORTED': f"carries reference {ref_b}; does NOT share this variant",
                'MIXED_AMBIGUOUS': f"conflicting reads (ref={ref_n}, alt={alt_n}, other={oth_n}); no call",
                'INSUFFICIENT_COVERAGE': f"only {depth} usable read(s); no call (not evidence of reference)",
                'UNRELIABLE_MAPPING': f"{frac_lowmq:.0%} of reads MAPQ<{a.min_mapq} / prelim mappability={mf}; no call",
            }[cls]
            if alt_dropped:
                interp += f" [NOTE: {alt_dropped} low-MAPQ/low-BQ read(s) carry {alt_b}; excluded]"
            if cls == 'INSUFFICIENT_COVERAGE' and depth == 1:
                interp += f" [single read shows {good[0][0]}]"
            fm = [g[1] for g in good]; fb = [g[3] for g in good]
            row.update(raw_primary_reads=raw, dup_or_secondary_excluded=excl, low_mapq_reads=lowmq,
                       low_bq_reads=lowbq, filtered_depth=depth, filtered_ref_reads=ref_n,
                       filtered_alt_reads=alt_n, filtered_other_reads=oth_n,
                       observed_alleles=';'.join(f"{k}:{v}" for k, v in sorted(cnt.items())) or '-',
                       alt_strands_fwd_rev=f"{af}/{ar}", ref_strands_fwd_rev=f"{rf}/{rr}",
                       mapq_median_filtered=(st.median(fm) if fm else ''), mapq_min_filtered=(min(fm) if fm else ''),
                       frac_primary_lowmapq=round(frac_lowmq, 3), bq_median_filtered=(st.median(fb) if fb else ''), alt_reads_in_filtered_out=alt_dropped,
                       near_read_end_reads=sum(1 for g in good if g[4] < 5),
                       softclipped_reads=sum(1 for g in good if g[5] > 0),
                       callable=('yes' if cls in ('DERIVED_SUPPORTED', 'ANCESTRAL_SUPPORTED', 'MIXED_AMBIGUOUS') else 'no'),
                       call_classification=cls, confidence=conf, interpretation=interp)
            w.writerow(row); out.flush()
            print(f"{smp['sample']} {s['id']} {pos} depth={depth} ref={ref_n} alt={alt_n} oth={oth_n} "
                  f"lowMQ={lowmq} -> {cls}/{conf} ({time.time()-t0:.1f}s)", file=sys.stderr)
    out.close(); rout.close()


if __name__ == '__main__':
    main()
