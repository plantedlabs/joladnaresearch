#!/usr/bin/env python3
"""
01_annotate_sites.py: sequence-context annotation for target chrY sites.

  * reference base check against the MD5-verified analysis-set chrY
  * UCSC hg38 soft-masking (lowercase = RepeatMasker/TRF) at site and +/-50 bp
  * UCSC genomicSuperDups overlap (segmental duplication, best fracMatch)
  * PAR membership
  * PRELIMINARY mappability: tiled 100-mers (10 tiles, site at read offsets
    5..95) aligned with `bwa mem -a` to chrX+chrY ONLY. Reports primary MAPQ,
    whether the primary hit is the true origin, and alternative hits with
    edit distance <= 2 over >= 90 aligned bp.

    *** This does NOT establish genome-wide uniqueness: autosomes, alt contigs
    and decoys are not in the index. A 'pass' means only "no near-identical
    copy on chrX/chrY for 100 bp reads". ***

  Flag:  fail    = >=50% of tiles MAPQ 0 or primary not at origin
         caution = any tile MAPQ < 20, or any alt hit NM <= 2, or segdup
                   fracMatch >= 0.98
         pass    = otherwise
"""
import argparse, csv, re, subprocess, statistics as st, tempfile, os

PAR = [(10001, 2781479), (56887903, 57217415)]


def faidx(fa, ctg, s, e):
    out = subprocess.run(['samtools', 'faidx', fa, f'{ctg}:{s}-{e}'], capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:])


def read_tsv(path):
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    hdr = lines[0].lstrip('#').split('\t')
    return [dict(zip(hdr, l.split('\t'))) for l in lines[1:] if not l.startswith('#')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sites', nargs='+', required=True)
    ap.add_argument('--refdir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    R = a.refdir
    ay = f'{R}/chrY_analysisset.fa'; uy = f'{R}/ucsc_hg38_chrY.fa'; bidx = f'{R}/chrXY_analysisset.fa'

    segd = []
    with open(f'{R}/genomicSuperDups.chrY.tsv') as f:
        for l in f:
            c = l.rstrip('\n').split('\t')
            segd.append((int(c[2]), int(c[3]), c[7], int(c[8]), float(c[26])))  # fracMatch = UCSC genomicSuperDups col 27

    sites = [s for p in a.sites for s in read_tsv(p)]
    # build tiles
    fq = tempfile.NamedTemporaryFile('w', suffix='.fq', delete=False)
    for s in sites:
        p = int(s['pos_grch38'])
        for k, off in enumerate(range(5, 100, 10)):
            start = p - off  # site at 0-based read offset `off`
            seq = faidx(ay, 'chrY', start, start + 99).upper()
            fq.write(f"@{s['id']}|{k}|{start}\n{seq}\n+\n{'I'*len(seq)}\n")
    fq.close()
    sam = subprocess.run(['bwa', 'mem', '-a', '-t', '2', bidx, fq.name], capture_output=True, text=True, check=True).stdout
    os.unlink(fq.name)
    tiles = {}
    for l in sam.splitlines():
        if l.startswith('@'):
            continue
        f = l.split('\t'); sid, k, start = f[0].split('|'); flag = int(f[1])
        nm = next((int(x[5:]) for x in f[11:] if x.startswith('NM:i:')), None)
        t = tiles.setdefault((sid, k), {'mapq': None, 'origin': False, 'alts': []})
        if flag & 0x900 == 0:
            t['mapq'] = int(f[4]); t['origin'] = (f[2] == 'chrY' and abs(int(f[3]) - int(start)) <= 5)
        elif nm is not None and nm <= 2 and sum(int(n) for n, o in re.findall(r'(\d+)([MIDNSHP=X])', f[5]) if o in 'M=X') >= 90:
            t['alts'].append(f"{f[2]}:{f[3]}(NM{nm})")

    cols = ['site_id', 'grch38_pos', 'config_ref', 'config_alt', 'analysisset_ref', 'ucsc_ref',
            'ref_check', 'in_PAR', 'repeatmasked_site', 'repeatmasked_frac_pm50', 'segdup_overlap',
            'segdup_best_fracMatch', 'tile_mapq_min', 'tile_mapq_median', 'tiles_primary_off_origin',
            'tiles_with_alt_hits_nm_le2', 'alt_hit_examples', 'mappability_flag', 'mappability_notes']
    with open(a.out, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=cols); w.writeheader()
        for s in sites:
            p = int(s['pos_grch38'])
            ab = faidx(ay, 'chrY', p, p).upper(); ub_raw = faidx(uy, 'chrY', p, p)
            win = faidx(uy, 'chrY', p - 50, p + 50)
            sd = [x for x in segd if x[0] < p <= x[1]]
            best = max((x[4] for x in sd), default=None)
            ts = [tiles.get((s['id'], str(k)), {'mapq': 0, 'origin': False, 'alts': []}) for k in range(10)]
            mq = [t['mapq'] if t['mapq'] is not None else 0 for t in ts]
            off = sum(1 for t in ts if not t['origin'])
            withalt = sum(1 for t in ts if t['alts'])
            ex = sorted({h for t in ts for h in t['alts']})[:4]
            notes = []
            if sum(1 for m in mq if m == 0) >= 5 or off >= 5:
                flag = 'fail'; notes.append('majority of 100-mer tiles ambiguous on chrX/chrY')
            elif min(mq) < 20 or withalt or (best is not None and best >= 0.98):
                flag = 'caution'
                if min(mq) < 20: notes.append(f'some tiles MAPQ {min(mq)}')
                if withalt: notes.append(f'{withalt}/10 tiles have near-identical alt hit')
                if best is not None and best >= 0.98: notes.append(f'segdup fracMatch {best:.3f}')
            else:
                flag = 'pass'
            if ub_raw.islower(): notes.append('site repeat-masked')
            notes.append('PRELIMINARY chrX/chrY-only test; not genome-wide')
            w.writerow(dict(site_id=s['id'], grch38_pos=p, config_ref=s['ref'], config_alt=s['alt'],
                            analysisset_ref=ab, ucsc_ref=ub_raw.upper(),
                            ref_check='match' if ab == s['ref'].upper() == ub_raw.upper() else 'MISMATCH',
                            in_PAR=any(a_ <= p <= b_ for a_, b_ in PAR),
                            repeatmasked_site=ub_raw.islower(),
                            repeatmasked_frac_pm50=round(sum(c.islower() for c in win) / len(win), 2),
                            segdup_overlap=bool(sd), segdup_best_fracMatch=(round(best, 4) if best else ''),
                            tile_mapq_min=min(mq), tile_mapq_median=st.median(mq),
                            tiles_primary_off_origin=off, tiles_with_alt_hits_nm_le2=withalt,
                            alt_hit_examples=';'.join(ex), mappability_flag=flag, mappability_notes='; '.join(notes)))


if __name__ == '__main__':
    main()
