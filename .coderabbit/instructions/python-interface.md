---
id: python-interface
name: Python user interface and post-processing
paths:
  - "pyphare/pyphare/**"
---

- `pharein` is the user's input language: a new or renamed keyword must be
  validated (unknown keywords rejected) and documented; a renamed keyword breaks
  existing input scripts, so flag it.
- `pharesee` reads runs written by older PHARE versions; flag a change that makes
  old HDF5 or VTKHDF files unreadable.
- Prefer numpy operations over Python loops on field or particle arrays.
