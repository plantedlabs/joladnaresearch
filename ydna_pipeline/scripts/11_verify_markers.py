#!/usr/bin/env python3
"""
11_verify_markers.py: pre-query verification of every marker used in the branch pilot.

For each marker (config/tree_markers.tsv + config/tierA_private_variants.tsv):
  1. Build/reference: the base at the GRCh38 position in the MD5-verified analysis-set chrY
     (ce3e31103314a704255f3cd90369ecce) equals the configured REF.
  2. Coordinate/name: YBrowse snps_hg38.vcf.gz (tabix, remote) has a record at that exact position whose
     ID contains the claimed name (tree markers) and whose REF/ALT equal the configured alleles.
  3. Orientation: YBrowse INFO AA (ancestral allele) == REF, so derived == ALT.
  4. Equivalence (tree markers): the SNP name is listed in the block of the claimed node on the live YFull
     tree page (visible names + the '+N SNPs' tooltip list).
Tier-A private SNPs are checked for 1-3 only. They are candidate markers, not validated branch definitions.

Output: results/audit/marker_verification.csv. Exit code 1 if any tree marker fails.
"""
import csv, html, os, re, subprocess, sys
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
YB = 'https://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz'
REF = 'ref/chrY_analysisset.fa'
YFULL_PAGES = {'E-M132': 'E-M132', 'E-Z958': 'E-Z958', 'E-CTS246': 'E-CTS246', 'E-CTS246 block': 'E-CTS246',
               'E-Z15174': 'E-Z15174', 'E-Z5991': 'E-Z15174', 'E-FT402211': 'E-Z15174', 'E-Z15196': 'E-Z15174'}


def rows(p):
    L = [l.rstrip('\n').split('\t') for l in open(p) if l.strip()]
    h = L[0][0].lstrip('#'), *L[0][1:]
    return [dict(zip(h, x)) for x in L[1:] if not x[0].startswith('#')]


def yfull_blocks(page):
    t = subprocess.run(['curl', '-sS', '--max-time', '60', f'https://www.yfull.com/tree/{page}/'],
                       capture_output=True, text=True).stdout
    blocks = {}
    # each node: <a ...>E-XXX</a> ... visible SNPs ... title="A * B * ..." (+N SNPs)
    for m in re.finditer(r'>(E-[A-Za-z0-9]+)</a>(.*?)(?=>E-[A-Za-z0-9]+</a>|$)', t, flags=re.S):
        node, body = m.group(1), m.group(2)[:6000]
        txt = html.unescape(re.sub(r'<[^>]+>', ' ', body))
        names = set(re.findall(r'[A-Z]{1,4}[0-9]+', txt.split('formed')[0]))
        for tt in re.findall(r'title="([^"]*\*[^"]*)"', body):
            names |= set(re.findall(r'[A-Z]{1,4}[0-9]+', html.unescape(tt)))
        blocks.setdefault(node, set()).update(names)
    return blocks


def main():
    tree = rows('config/tree_markers.tsv'); tier = rows('config/tierA_private_variants.tsv')
    cache = {}
    out = []; fail = 0
    for kind, mk in [('tree_marker', m) for m in tree] + [('tierA_candidate', m) for m in tier]:
        pos = int(mk['pos_grch38']); ref, alt = mk['ref'], mk['alt']; name = mk['id']
        rb = subprocess.run(['samtools', 'faidx', REF, f'chrY:{pos}-{pos}'], capture_output=True, text=True,
                            check=True).stdout.split('\n')[1].upper()
        yb = subprocess.run(['tabix', YB, f'chrY:{pos}-{pos}'], capture_output=True, text=True, timeout=120).stdout
        recs = [l.split('\t') for l in yb.splitlines() if l.split('\t')[1] == str(pos)]
        match = [r for r in recs if r[3] == ref and r[4] == alt]
        ids = ';'.join(r[2] for r in recs)
        aa = ''
        if match:
            m = re.search(r'AA=([ACGT])', match[0][7]); aa = m.group(1) if m else ''
        name_ok = (kind == 'tierA_candidate') or any(name in r[2].split(',') for r in match)
        equiv = 'n/a (candidate marker)'
        if kind == 'tree_marker' and mk['yfull_node'] not in YFULL_PAGES:
            equiv = 'FTDNA-only node (not on YFull; YF133256 = E-Z15174* there): equivalence not checkable'
        elif kind == 'tree_marker':
            node = mk['yfull_node']; page = YFULL_PAGES[node]
            if page not in cache:
                cache[page] = yfull_blocks(page)
            target = node.replace(' block', '')
            equiv = 'listed in YFull block' if name in cache[page].get(target, set()) else 'NOT FOUND in YFull block'
        ok = rb == ref and bool(match) and aa == ref and name_ok and (kind != 'tree_marker' or equiv.startswith(('listed', 'FTDNA-only')))
        if not ok and kind == 'tree_marker':
            fail += 1
        out.append(dict(marker=name, kind=kind, grch38_pos=pos, config_ref=ref, config_alt=alt, analysisset_ref=rb,
                        ybrowse_ids_at_pos=ids, ybrowse_ref_alt_match=bool(match), ybrowse_AA=aa,
                        orientation='ancestral=REF, derived=ALT' if aa == ref else 'CHECK', name_at_coordinate=name_ok,
                        yfull_node=mk.get('yfull_node', ''), yfull_equivalence=equiv, verdict='PASS' if ok else 'FAIL'))
        print(f"{name:9s} {kind:16s} {pos:>9} ref {rb}=={ref} yb[{ids[:30]}] AA={aa} {equiv} -> {'PASS' if ok else 'FAIL'}", file=sys.stderr)
    with open('results/audit/marker_verification.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
