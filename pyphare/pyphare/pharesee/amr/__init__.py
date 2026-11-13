"""Flat, per-level storage of AMR fields for pharesee.

One array per (time, level, component), node centred, no ghost nodes, the same
for phareh5 and pharevtkhdf output. Sub-modules are imported on demand; nothing
is imported here so that importing ``pyphare.pharesee`` stays cheap.
"""
