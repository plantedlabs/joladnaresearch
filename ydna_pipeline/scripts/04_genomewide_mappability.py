#!/usr/bin/env python3
"""
04_genomewide_mappability.py: query PUBLISHED genome-wide mappability and
problematic-region tracks at each site (remote byte-range access; nothing
large is downloaded).

Tracks (UCSC hg38 gbdb):
  * Umap multi-read mappability, k = 24/36/50/100 (Karimzadeh, Ernst, Kundaje
    & Hoffman 2018, NAR 46:e120), hoffmanMappability/k{K}.Umap.MultiTrackMappability.bw
    Value = fraction of k-mers overlapping the position that are uniquely
    mappable in the whole GRCh38 genome (incl. chrY). 1.0 = every overlapping
    k-mer is unique genome-wide; 0 = none are.
  * Umap single-read unique k100 regions, k100.Unique.Mappability.bb
  * ENCODE blacklist v2 (Amemiya, Kundaje & Boyle 2019), problematic/encBlacklist.bb
  * GRC exclusions, problematic/grcExclusions.bb

Output: per site: value at the site, and min/mean over +/-50 bp (a 100 bp read
covering the site spans this window).
"""
import argparse, csv, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bbi_remote import BBI

G = 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/'


def read_tsv(path):
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    hdr = lines[0].lstrip('#').split('\t')
    return [dict(zip(hdr, l.split('\t'))) for l in lines[1:] if not l.startswith('#')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sites', nargs='+', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    umap = {k: BBI(f'{G}hoffmanMappability/k{k}.Umap.MultiTrackMappability.bw') for k in (24, 36, 50, 100)}
    uniq100 = BBI(f'{G}hoffmanMappability/k100.Unique.Mappability.bb')
    bl = BBI(f'{G}problematic/encBlacklist.bb')
    grc = BBI(f'{G}problematic/grcExclusions.bb')
    cols = ['site_id', 'grch38_pos'] + [f'umap_k{k}_{m}' for k in (24, 36, 50, 100) for m in ('site', 'min_pm50', 'mean_pm50')] + \
           ['k100_unique_region', 'encode_blacklist', 'grc_exclusion']
    with open(a.out, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=cols); w.writeheader()
        for s in [x for p in a.sites for x in read_tsv(p)]:
            p1 = int(s['pos_grch38']); p0 = p1 - 1
            row = {'site_id': s['id'], 'grch38_pos': p1}
            for k, r in umap.items():
                iv = r.query('chrY', p0 - 50, p1 + 50)
                # expand to per-base values (missing = 0: Umap omits unmappable bases)
                vals = {}
                for st, en, v in iv:
                    for b in range(max(st, p0 - 50), min(en, p1 + 50)):
                        vals[b] = v
                win = [vals.get(b, 0.0) for b in range(p0 - 50, p1 + 50)]
                row[f'umap_k{k}_site'] = round(vals.get(p0, 0.0), 3)
                row[f'umap_k{k}_min_pm50'] = round(min(win), 3)
                row[f'umap_k{k}_mean_pm50'] = round(sum(win) / len(win), 3)
            row['k100_unique_region'] = bool(uniq100.query('chrY', p0, p1))
            row['encode_blacklist'] = ';'.join(f'{x[2]}:{x[0]+1}-{x[1]}' for x in bl.query('chrY', p0, p1)) or ''
            row['grc_exclusion'] = ';'.join(f'{x[2]}:{x[0]+1}-{x[1]}' for x in grc.query('chrY', p0, p1)) or ''
            w.writerow(row); fo.flush()
            print(s['id'], {k: row[f'umap_k{k}_site'] for k in (24, 36, 50, 100)}, row['k100_unique_region'],
                  row['encode_blacklist'], row['grc_exclusion'], file=sys.stderr)


if __name__ == '__main__':
    main()
