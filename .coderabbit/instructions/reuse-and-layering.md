---
id: reuse-and-layering
name: Reuse existing PHARE code, respect the layers
paths:
  - "src/**"
  - "pyphare/**"
---

For every new helper, type, free function or loop, search the repository for an
existing equivalent before accepting it. Look first in:
- `src/core/def.hpp`: `core::isUsable(...)`, `core::isSettable(...)`, `NO_DISCARD`;
- `src/core/utilities/`: `box/`, `point/`, `index/`, `range/`, `meta/`,
  `algorithm.hpp`, `types.hpp`, `span.hpp`, `indexer.hpp`, `cellmap.hpp`;
- `src/amr/utilities/box/`, and the `core/data/` types (`GridLayout`, `NdArray`,
  fields, `ParticleArray`);
- Python: `pyphare/pyphare/core/` (`box.py`, `gridlayout.py`, `operators.py`,
  `phare_utilities.py`).
Flag a duplicate with both locations (new file:line and existing file:line) and
say whether the existing one is a drop-in or needs a small extension. Flag a
second abstraction next to an existing one (another box type, another usability
check beside `isUsable`).

Layers, bottom-up: `core` (physics, no SAMRAI) < `amr` (only layer that uses
SAMRAI) < `simulator` < `python3` bindings < `pyphare`. Flag SAMRAI types reaching
`core` indirectly, through template parameters or forwarded headers.

Runtime choices use factories and strategies (`MessengerFactory`,
`TaggerFactory`, load balancers). A new variant is a new strategy class; flag a
PR that edits or removes existing strategies to add one.

New (dimension, interpolation order, refined particle count) combinations must be
listed in `possibleSimulators()` (`src/simulator/simulator_def.hpp`); each one
adds a compiled module and compile time.
