#!/usr/bin/env python3
"""
bbi_remote.py: minimal, dependency-free reader for remote bigWig / bigBed files
(UCSC BBI format, Kent et al. 2010) using HTTP byte-range requests via curl.

Only the chromosome B+ tree, the R-tree nodes overlapping the query, and the
relevant zlib-compressed data blocks are downloaded, typically a few kB-MB,
so whole-genome tracks (e.g. 825 MB Umap k100) can be queried at single sites.

Usage as a library:
    r = BBI(url); r.query('chrY', start0, end0) -> list of (start0, end0, value_or_rest)
CLI:
    python3 bbi_remote.py URL chrY 10973239 10973240
"""
import struct, subprocess, sys, zlib

BW_MAGIC, BB_MAGIC = 0x888FFC26, 0x8789F2EB


class BBI:
    def __init__(self, url, timeout=120):
        self.url, self.timeout = url, timeout
        self._cache = {}
        h = self._get(0, 64)
        (self.magic, self.version, self.zoom, self.chrom_tree, self.full_data, self.full_index,
         self.field_count, self.defined_fc, self.autosql, self.total_summary,
         self.uncomp_buf) = struct.unpack('<IHHQQQHHQQI', h[:56])
        if self.magic not in (BW_MAGIC, BB_MAGIC):
            raise ValueError(f'not a little-endian bigWig/bigBed: {self.magic:#x}')
        self.is_bw = self.magic == BW_MAGIC
        self.chroms = self._read_chrom_tree()

    def _get(self, off, n):
        key = (off, n)
        if key not in self._cache:
            p = subprocess.run(['curl', '-sS', '--fail', '--retry', '3', '--max-time', str(self.timeout),
                                '-r', f'{off}-{off + n - 1}', self.url], capture_output=True, check=True)
            if len(p.stdout) != n:
                raise IOError(f'short read {len(p.stdout)}/{n} at {off}')
            self._cache[key] = p.stdout
        return self._cache[key]

    def _read_chrom_tree(self):
        magic, bs, ks, vs, n, _ = struct.unpack('<IIIIQQ', self._get(self.chrom_tree, 32))
        chroms = {}

        def node(off):
            is_leaf, _, cnt = struct.unpack('<BBH', self._get(off, 4))
            item = ks + 8
            data = self._get(off + 4, cnt * item)
            for i in range(cnt):
                rec = data[i * item:(i + 1) * item]
                key = rec[:ks].rstrip(b'\0').decode()
                if is_leaf:
                    cid, csize = struct.unpack('<II', rec[ks:])
                    chroms[key] = (cid, csize)
                else:
                    node(struct.unpack('<Q', rec[ks:])[0])
        node(self.chrom_tree + 32)
        return chroms

    def _blocks(self, cid, start, end):
        out = []

        def ov(sc, sb, ec, eb):
            return (sc, sb) < (cid, end) and (ec, eb) > (cid, start)

        def node(off):
            is_leaf, _, cnt = struct.unpack('<BBH', self._get(off, 4))
            item = 32 if is_leaf else 24
            data = self._get(off + 4, cnt * item)
            for i in range(cnt):
                rec = data[i * item:(i + 1) * item]
                sc, sb, ec, eb = struct.unpack('<IIII', rec[:16])
                if not ov(sc, sb, ec, eb):
                    continue
                if is_leaf:
                    out.append(struct.unpack('<QQ', rec[16:]))
                else:
                    node(struct.unpack('<Q', rec[16:])[0])
        node(self.full_index + 48)
        return out

    def query(self, chrom, start, end):
        if chrom not in self.chroms:
            return []
        cid = self.chroms[chrom][0]
        res = []
        for off, size in self._blocks(cid, start, end):
            raw = self._get(off, size)
            buf = zlib.decompress(raw) if self.uncomp_buf else raw
            if self.is_bw:
                c, s0, e0, step, span, typ, _, cnt = struct.unpack('<IIIIIBBH', buf[:24])
                p = 24
                for i in range(cnt):
                    if typ == 1:
                        s, e, v = struct.unpack('<IIf', buf[p:p + 12]); p += 12
                    elif typ == 2:
                        s, v = struct.unpack('<If', buf[p:p + 8]); e = s + span; p += 8
                    else:
                        v, = struct.unpack('<f', buf[p:p + 4]); s = s0 + i * step; e = s + span; p += 4
                    if c == cid and s < end and e > start:
                        res.append((s, e, v))
            else:
                p = 0
                while p < len(buf):
                    c, s, e = struct.unpack('<III', buf[p:p + 12]); p += 12
                    z = buf.index(b'\0', p); rest = buf[p:z].decode(); p = z + 1
                    if c == cid and s < end and e > start:
                        res.append((s, e, rest))
        return sorted(res)


if __name__ == '__main__':
    url, chrom, s, e = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    for r in BBI(url).query(chrom, s, e):
        print(chrom, *r, sep='\t')
