---
id: python-interface
name: Python user interface and post-processing
paths:
  - "pyphare/pyphare/**"
---

- `pharein` is the user's input language. The gate flags removed keywords. You
  check: a new keyword is validated (unknown keywords rejected) and documented; a
  changed default value or a changed meaning of an existing keyword (same name,
  different physics) is flagged, because input scripts keep running and give
  different results.
- `pharesee` reads runs written by older PHARE versions; flag a change that makes
  old HDF5 or VTKHDF files unreadable.
- Prefer numpy operations over Python loops on field or particle arrays.
