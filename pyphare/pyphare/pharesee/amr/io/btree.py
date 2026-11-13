"""Walker of the HDF5 version-1 chunk B-tree, for raw reads of chunked datasets.

Why: h5py reads a large chunked dataset at ~200 MB/s cold on network file
systems, and every HDF5 call that lists chunk addresses (``get_chunk_info``,
``chunk_iter``, ``h5o.get_info``) scans the whole chunk index (minutes on a
100 GB file). The chunks of one dump of one level sit in one contiguous byte
region of the file, interleaved with the B-tree leaf nodes that index them.
So we walk the B-tree down to the level just above the leaves (a few dozen
nodes), which gives the leaf addresses, hence the byte span of the dump; the
caller reads that span with large parallel reads and parses the leaves from
memory (``parse_leaves``) to locate each chunk inside the buffer.

Supported, and checked (anything else raises ``Unsupported`` and the caller
falls back to h5py):
- superblock version 0 or 1, base address 0, offsets and lengths of 8 bytes;
  the B-tree order K is read from the superblock (32 when not stored, v0);
- object header version 1 for the dataset;
- data layout message version 3, chunked class (v1 B-tree chunk index);
  version 4 (``libver="latest"``: other chunk indexes) is refused;
- no filter pipeline message (no compression, shuffle, checksum);
- chunks of exactly ``(chunk_rows, 1)`` elements of a little-endian dtype,
  each stored uncompressed (on-disk size == chunk bytes, filter mask 0).

Format reference: HDF5 File Format Specification version 2/3, sections II.A
(superblock), III.A.1 (v1 B-trees), IV.A.1.a (v1 object headers), IV.A.2.i
(data layout message).
"""

from __future__ import annotations

import struct
import sys

import numpy as np

from . import rawio

_SIGNATURE = b"\x89HDF\r\n\x1a\n"
_MSG_NIL = 0x00
_MSG_LAYOUT = 0x08
_MSG_FILTERS = 0x0B
_MSG_CONTINUATION = 0x10
_UNDEF = 0xFFFFFFFFFFFFFFFF


class Unsupported(Exception):
    """The file uses an HDF5 feature this walker does not handle.

    Private to the io layer: the reader catches it and falls back to h5py.
    The message says which check failed.
    """


def read_superblock(fd):
    """Parse the superblock at offset 0. Returns dict(version, chunk_k).

    Raises Unsupported unless: signature at 0 (no user block), version 0/1,
    8-byte offsets and lengths, base address 0.
    """
    head = rawio.read_ranges(fd, [0], [96], nthreads=1).tobytes()
    if head[:8] != _SIGNATURE:
        raise Unsupported("no HDF5 signature at offset 0 (user block?)")
    version = head[8]
    if version not in (0, 1):
        raise Unsupported(f"superblock version {version} (only 0 and 1 are handled)")
    size_offsets, size_lengths = head[13], head[14]
    if (size_offsets, size_lengths) != (8, 8):
        raise Unsupported(f"sizes of offsets/lengths {(size_offsets, size_lengths)} != (8, 8)")
    if version == 0:
        chunk_k = 32  # not stored in v0: HDF5 default for the chunk B-tree
        base_at = 24
    else:
        chunk_k = struct.unpack_from("<H", head, 24)[0]
        base_at = 28
    base = struct.unpack_from("<Q", head, base_at)[0]
    if base != 0:
        raise Unsupported(f"superblock base address {base} != 0")
    if not 1 <= chunk_k <= 32768:
        raise Unsupported(f"chunk B-tree K = {chunk_k} out of range")
    return dict(version=version, chunk_k=int(chunk_k))


def read_layout(fd, header_addr, file_size):
    """Parse the v1 object header at `header_addr`; return the chunk layout.

    Returns dict(btree=address of the chunk B-tree root, dims=tuple of chunk
    dimensions including the trailing element size in bytes).
    Raises Unsupported for an object header other than v1, a layout message
    other than version 3 chunked, or the presence of a filter pipeline.
    """
    if not 0 <= header_addr < file_size:
        raise Unsupported(f"object header address {header_addr} outside the file")
    head = rawio.read_ranges(fd, [header_addr], [16], nthreads=1).tobytes()
    if head[:4] == b"OHDR":
        raise Unsupported("object header version 2 (file written with a recent libver)")
    version, _, _nmsg, _refcount, size = struct.unpack_from("<BBHII", head, 0)
    if version != 1:
        raise Unsupported(f"object header version {version}")
    blocks = [(header_addr + 16, size)]
    layout = None
    seen = 0
    while blocks:
        addr, length = blocks.pop(0)
        seen += 1
        if seen > 64 or length > (1 << 20) or addr + length > file_size:
            raise Unsupported("object header continuation chain looks corrupt")
        body = rawio.read_ranges(fd, [addr], [length], nthreads=1).tobytes()
        p = 0
        while p + 8 <= length:
            mtype, msize, _flags = struct.unpack_from("<HHB", body, p)
            msg = body[p + 8 : p + 8 + msize]
            if mtype == _MSG_CONTINUATION:
                blocks.append(struct.unpack_from("<QQ", msg, 0))
            elif mtype == _MSG_FILTERS:
                raise Unsupported("dataset has a filter pipeline (compression/shuffle/checksum)")
            elif mtype == _MSG_LAYOUT:
                layout = _parse_layout(msg)
            p += 8 + msize
    if layout is None:
        raise Unsupported("no data layout message found")
    return layout


def _parse_layout(msg):
    version = msg[0]
    if version != 3:
        raise Unsupported(f"data layout message version {version} (not a v1 B-tree index)")
    layout_class = msg[1]
    if layout_class != 2:
        raise Unsupported(f"data layout class {layout_class} is not 'chunked'")
    ndims = msg[2]
    btree_addr = struct.unpack_from("<Q", msg, 3)[0]
    dims = struct.unpack_from(f"<{ndims}I", msg, 11)
    return dict(btree=int(btree_addr), dims=tuple(int(d) for d in dims))


class ChunkBTree:
    """Chunk index of one 2-D chunked dataset ``(nrows, ncols)``, chunks ``(rows, 1)``.

    Construction reads the superblock and the dataset object header (3-4 small
    reads) and validates everything listed in the module docstring.

    fd           open file descriptor of the HDF5 file
    header_addr  address of the dataset object header (``links.get_info(name).u``)
    dtype        numpy dtype of the dataset elements
    chunk_shape  dataset chunk shape as h5py reports it, must be ``(rows, 1)``
    ncols        number of columns of the dataset

    Attributes: ``root`` (B-tree root address), ``k`` (B-tree order: 2K
    children per node at most), ``nodesize`` (bytes of one node), ``chunk_rows``,
    ``chunk_bytes``, ``ncols``, ``dtype``.
    Raises Unsupported.
    """

    def __init__(self, fd, header_addr, dtype, chunk_shape, ncols):
        file_size = rawio.stat_signature(fd)[0]
        dtype = np.dtype(dtype)
        if sys.byteorder != "little" or dtype.byteorder not in ("<", "=", "|"):
            raise Unsupported(f"dtype {dtype.str} is not little-endian on a little-endian host")
        if chunk_shape is None or len(chunk_shape) != 2 or chunk_shape[1] != 1:
            raise Unsupported(f"chunk shape {chunk_shape} is not (rows, 1)")
        superblock = read_superblock(fd)
        layout = read_layout(fd, header_addr, file_size)
        expect = (int(chunk_shape[0]), 1, dtype.itemsize)
        if layout["dims"] != expect:
            raise Unsupported(f"layout chunk dims {layout['dims']} != {expect}")
        self.file_size = file_size
        self.root = layout["btree"]
        if self.root == _UNDEF or not 0 <= self.root < file_size:
            raise Unsupported("chunk B-tree root undefined or outside the file")
        self.k = superblock["chunk_k"]
        self.ndims = len(layout["dims"])  # rank + 1 (element size)
        self.chunk_rows = int(chunk_shape[0])
        self.ncols = int(ncols)
        self.dtype = dtype
        self.chunk_bytes = self.chunk_rows * dtype.itemsize
        keysize = 8 + 8 * self.ndims
        self.entry = np.dtype(
            [("size", "<u4"), ("mask", "<u4"), ("off", "<u8", (self.ndims,)), ("child", "<u8")]
        )
        assert self.entry.itemsize == keysize + 8
        self.key = np.dtype([("size", "<u4"), ("mask", "<u4"), ("off", "<u8", (self.ndims,))])
        self.nodesize = 24 + (2 * self.k + 1) * keysize + 2 * self.k * 8

    # ------------------------------------------------------------ nodes

    def _check_header(self, raw, level=None):
        if raw[:4] != b"TREE" or raw[4] != 1:
            raise Unsupported("bad B-tree node signature or node type (not raw data)")
        node_level = raw[5]
        n = struct.unpack_from("<H", raw, 6)[0]
        if n > 2 * self.k:
            raise Unsupported(f"B-tree node with {n} entries > 2K = {2 * self.k}")
        if level is not None and node_level != level:
            raise Unsupported(f"B-tree node level {node_level}, expected {level}")
        return node_level, n

    def _parse_internal(self, raw):
        """-> (level, child row lower bounds (n+1,), child addresses (n,))."""
        level, n = self._check_header(raw)
        entries = np.frombuffer(raw, dtype=self.entry, count=n, offset=24)
        last = np.frombuffer(raw, dtype=self.key, count=1, offset=24 + n * self.entry.itemsize)
        rows = np.concatenate([entries["off"][:, 0], last["off"][:, 0]]).astype(np.int64)
        return level, rows, entries["child"].astype(np.int64)

    def _read_nodes(self, fd, addrs, executor=None):
        addrs = np.asarray(addrs, dtype=np.int64)
        if ((addrs < 0) | (addrs + self.nodesize > self.file_size)).any():
            raise Unsupported("B-tree node address outside the file")
        buf = rawio.read_ranges(
            fd, addrs, np.full(len(addrs), self.nodesize), executor=executor,
            nthreads=None if executor is not None else min(8, rawio.io_threads()),
        )
        return buf.reshape(len(addrs), self.nodesize)

    def leaf_addresses(self, fd, row_lo, row_hi, executor=None):
        """Addresses (sorted, unique, int64) of the leaves indexing rows [row_lo, row_hi).

        Reads only the internal nodes (the levels above the leaves) whose key
        range overlaps the rows: ~0.01 s on a 100 GB file. A conservative
        overlap test may keep a neighbouring leaf; that only costs a few bytes.
        """
        lo = (int(row_lo) // self.chunk_rows) * self.chunk_rows
        hi = int(row_hi)
        frontier = np.array([self.root], dtype=np.int64)
        root_raw = self._read_nodes(fd, frontier, executor)[0].tobytes()
        root_level, _ = self._check_header(root_raw)
        if root_level == 0:
            return frontier
        level = root_level
        leaves = []
        depth = 0
        while len(frontier):
            depth += 1
            if depth > 32:
                raise Unsupported("B-tree deeper than 32 levels: corrupt?")
            raws = self._read_nodes(fd, frontier, executor)
            nxt = []
            for raw in raws:
                raw = raw.tobytes()
                node_level, rows, children = self._parse_internal(raw)
                if node_level != level:
                    raise Unsupported("B-tree nodes of one depth have different levels")
                # child i holds keys in [rows[i], rows[i+1]]: keep it if it overlaps [lo, hi)
                keep = (rows[:-1] < hi) & (rows[1:] >= lo)
                nxt.append(children[keep])
            nxt = np.concatenate(nxt) if nxt else np.empty(0, np.int64)
            if level == 1:
                leaves = nxt
                break
            frontier = nxt
            level -= 1
        return np.unique(np.asarray(leaves, dtype=np.int64))

    def level_span(self, fd, row_lo, row_hi, margin=1 << 20, executor=None):
        """Byte span holding the chunks of rows [row_lo, row_hi).

        Returns dict:
        leaves    all leaf addresses (sorted) indexing these rows
        main      the leaves inside [span_lo, span_hi) (largest dense cluster)
        others    leaves outside it (to be read one by one)
        span_lo, span_hi  byte span = main cluster +- `margin` bytes, clipped
                  to the file. Chunks are interleaved with their leaves, so this
                  span holds almost all chunks; the few outside (e.g. the first
                  chunk row, allocated with the previous dump) are read singly.
        """
        leaves = self.leaf_addresses(fd, row_lo, row_hi, executor)
        if len(leaves) == 0:
            raise Unsupported("no B-tree leaf found for the requested rows")
        gap = max(1 << 20, 4 * 2 * self.k * self.chunk_bytes * max(1, self.ncols))
        cuts = np.nonzero(np.diff(leaves) > gap)[0] + 1
        clusters = np.split(leaves, cuts)
        main = max(clusters, key=len)
        others = np.setdiff1d(leaves, main)
        span_lo = max(0, int(main[0]) - int(margin))
        span_hi = min(self.file_size, int(main[-1]) + self.nodesize + int(margin))
        return dict(leaves=leaves, main=main, others=others, span_lo=span_lo, span_hi=span_hi)

    # ------------------------------------------------------------ leaves

    def parse_leaves(self, buf, rel):
        """Parse leaf nodes found in memory.

        buf  uint8 ndarray holding the nodes (e.g. a window of the file)
        rel  int64 offsets of the leaf nodes inside `buf`

        Returns (rows, cols, addrs) int64 arrays, one entry per chunk: first
        row and column of the chunk in the dataset, byte address in the file.
        Vectorised (one gather for all leaves). Raises Unsupported when a
        node is not a raw-data leaf, a chunk is filtered or not of the
        expected size, or an address is outside the file.
        """
        rel = np.asarray(rel, dtype=np.int64)
        if len(rel) == 0:
            e = np.empty(0, np.int64)
            return e, e, e
        if (rel < 0).any() or (rel + self.nodesize > len(buf)).any():
            raise Unsupported("leaf node outside the buffer")
        nodes = buf[rel[:, None] + np.arange(self.nodesize)]  # (nleaves, nodesize) uint8
        if (nodes[:, :4] != np.frombuffer(b"TREE", np.uint8)).any():
            raise Unsupported("bad B-tree leaf signature")
        if (nodes[:, 4] != 1).any() or (nodes[:, 5] != 0).any():
            raise Unsupported("B-tree node is not a raw-data leaf")
        counts = nodes[:, 6].astype(np.int64) | (nodes[:, 7].astype(np.int64) << 8)
        if (counts > 2 * self.k).any():
            raise Unsupported("B-tree leaf with more than 2K entries")
        width = 2 * self.k * self.entry.itemsize
        table = np.ascontiguousarray(nodes[:, 24 : 24 + width]).view(self.entry)
        table = table.reshape(len(rel), 2 * self.k)
        used = np.arange(2 * self.k)[None, :] < counts[:, None]
        ent = table[used]
        if (ent["mask"] != 0).any() or (ent["size"] != self.chunk_bytes).any():
            raise Unsupported("filtered or odd-sized chunk")
        off = ent["off"].astype(np.int64)
        rows, cols = off[:, 0], off[:, 1]
        if (off[:, 2:] != 0).any() or (rows % self.chunk_rows != 0).any():
            raise Unsupported("chunk key not aligned on the chunk grid")
        if ((cols < 0) | (cols >= self.ncols)).any():
            raise Unsupported("chunk column outside the dataset")
        addrs = ent["child"].astype(np.int64)
        if ((addrs < 0) | (addrs + self.chunk_bytes > self.file_size)).any():
            raise Unsupported("chunk address outside the file")
        return rows, cols, addrs

    def read_leaves(self, fd, addrs, executor=None):
        """Read leaf nodes at `addrs` one by one and parse them (see parse_leaves)."""
        addrs = np.asarray(addrs, dtype=np.int64)
        if len(addrs) == 0:
            return self.parse_leaves(np.empty(0, np.uint8), addrs)
        raws = self._read_nodes(fd, addrs, executor).reshape(-1)
        return self.parse_leaves(raws, np.arange(len(addrs), dtype=np.int64) * self.nodesize)
