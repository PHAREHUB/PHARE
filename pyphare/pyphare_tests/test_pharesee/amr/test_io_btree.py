"""Tests of pharesee.amr.io.btree, the HDF5 v1 chunk B-tree walker.

The walker's results are compared with HDF5's own chunk listing (chunk_iter,
fine on these small files; the reader itself never calls it).
"""

import os
import struct
import unittest

import numpy as np

from pyphare.pharesee.amr.io import btree, rawio
from pyphare.pharesee.amr.io.vtkhdf import VtkhdfBackend
from pyphare_tests.test_pharesee.amr._common import h5py, tmp_dir
from pyphare_tests.test_pharesee.amr.vtkhdf_files import tiling, write_vtkhdf


def chunk_map(dset):
    """{(row, col): byte address} from HDF5 itself."""
    out = {}

    def cb(info):
        out[tuple(int(x) for x in info.chunk_offset)] = int(info.byte_offset)

    dset.id.chunk_iter(cb)
    return out


def as_map(rows, cols, addrs):
    """{(row, col): byte address} from the walker's three arrays."""
    return {(r, c): a for r, c, a in zip(rows.tolist(), cols.tolist(), addrs.tolist())}


class BTreeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # small chunks: a deep tree (root level >= 2) in a small file
        cls.deep = str(tmp_dir(cls) / "deep.vtkhdf")
        steps = [[tiling([240, 160], [20, 16])] for _ in range(3)]
        write_vtkhdf(cls.deep, steps, chunk_rows=64)

    def open_tree(self, path, level=0):
        f = h5py().File(path, "r")
        self.addCleanup(f.close)
        dset = f[f"VTKHDF/Level{level}/PointData/data"]
        fd = rawio.open_readonly(path)
        self.addCleanup(os.close, fd)
        bt = btree.ChunkBTree(
            fd, VtkhdfBackend._header_addr(f, level), dset.dtype, dset.chunks, dset.shape[1]
        )
        return f, dset, fd, bt

    def test_all_leaves_give_hdf5_chunk_map(self):
        _, dset, fd, bt = self.open_tree(self.deep)
        leaves = bt.leaf_addresses(fd, 0, dset.shape[0])
        rows, cols, addrs = bt.read_leaves(fd, leaves)
        mine = as_map(rows, cols, addrs)
        self.assertEqual(len(mine), len(rows))  # no chunk twice
        self.assertEqual(mine, chunk_map(dset))

    def test_sub_range_leaves_cover_the_range(self):
        _, dset, fd, bt = self.open_tree(self.deep)
        ref = chunk_map(dset)
        nrows = dset.shape[0]
        for lo, hi in [(0, 1), (5000, 5100), (nrows // 3, 2 * nrows // 3), (nrows - 7, nrows)]:
            leaves = bt.leaf_addresses(fd, lo, hi)
            got = as_map(*bt.read_leaves(fd, leaves))
            need = {k: v for k, v in ref.items() if k[0] + 64 > lo and k[0] < hi}
            self.assertTrue(need.items() <= got.items(), (lo, hi))
            # the walk does not read the whole index
            self.assertLess(len(leaves), len(bt.leaf_addresses(fd, 0, nrows)) // 2 + 2)

    def test_parse_leaves_from_memory_equals_single_reads(self):
        _, dset, fd, bt = self.open_tree(self.deep)
        span = bt.level_span(fd, 0, dset.shape[0])
        lo, hi = span["span_lo"], span["span_hi"]
        buf = rawio.read_ranges(fd, [lo], [hi - lo])
        a = bt.parse_leaves(buf, span["main"] - lo)
        b = bt.read_leaves(fd, span["main"])
        for x, y in zip(a, b):
            np.testing.assert_array_equal(x, y)

    def test_corrupt_leaf_is_refused(self):
        _, _, fd, bt = self.open_tree(self.deep)
        leaf = int(bt.leaf_addresses(fd, 0, 10)[0])
        buf = rawio.read_ranges(fd, [leaf], [bt.nodesize]).copy()
        # (byte, new value): signature, node type (not a leaf), filter mask of chunk 0
        for byte, value in ((0, ord("X")), (5, 1), (24 + 4, 1)):
            bad = buf.copy()
            bad[byte] = value
            with self.subTest(byte=byte), self.assertRaises(btree.Unsupported):
                bt.parse_leaves(bad, [0])
        with self.assertRaises(btree.Unsupported):
            bt.parse_leaves(buf, [10])  # node beyond the buffer


class RefusalTest(unittest.TestCase):
    """Everything the walker does not know raises Unsupported."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tmp_dir(cls)
        cls.steps = [[tiling([60, 40], [20, 20])]]

    def tree(self, path):
        with h5py().File(path, "r") as f:
            dset = f["VTKHDF/Level0/PointData/data"]
            fd = rawio.open_readonly(path)
            try:
                return btree.ChunkBTree(
                    fd, VtkhdfBackend._header_addr(f, 0), dset.dtype, dset.chunks, dset.shape[1]
                )
            finally:
                os.close(fd)

    def write(self, name, **kw):
        path = str(self.tmp / name)
        write_vtkhdf(path, self.steps, **kw)
        return path

    def test_big_endian(self):
        with self.assertRaisesRegex(btree.Unsupported, "endian"):
            self.tree(self.write("be.vtkhdf", dtype=">f8"))

    def test_superblock_parsing(self):
        path = self.tmp / "sb.bin"

        def check(version, k=None, base=0, sizes=(8, 8), sig=btree._SIGNATURE):
            head = bytearray(96)
            head[:8] = sig
            head[8] = version
            head[13], head[14] = sizes
            if version == 1:
                struct.pack_into("<H", head, 24, k)
                struct.pack_into("<Q", head, 28, base)
            else:
                struct.pack_into("<Q", head, 24, base)
            path.write_bytes(head)
            fd = rawio.open_readonly(path)
            try:
                return btree.read_superblock(fd)
            finally:
                os.close(fd)

        self.assertEqual(check(0)["chunk_k"], 32)
        self.assertEqual(check(1, k=16)["chunk_k"], 16)
        for kw in (dict(version=2), dict(version=3), dict(version=0, base=512),
                   dict(version=0, sizes=(4, 8)), dict(version=0, sig=b"notHDF5!")):
            with self.assertRaises(btree.Unsupported, msg=str(kw)):
                check(**kw)


if __name__ == "__main__":
    unittest.main()
