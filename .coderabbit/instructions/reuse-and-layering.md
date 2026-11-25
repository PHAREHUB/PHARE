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

Layers, bottom-up: `core` (physics, never SAMRAI) < `amr` (wraps SAMRAI) <
`simulator` < `python3` bindings < `pyphare`. A few files above `amr` use SAMRAI
directly (`simulator.hpp`, `python3/cpp_etc.cpp`, `diagnostic_model_view.hpp`):
not a finding by itself, but prefer the `amr` wrappers in new code. Flag SAMRAI
types reaching `core`, including indirectly, through template parameters or
forwarded headers.

Runtime choices use factories and strategies (`MessengerFactory`,
`TaggerFactory`, load balancers). A new variant is a new strategy class; flag a
PR that edits or removes existing strategies to add one.

New (dimension, interpolation order, refined particle count) combinations must be
listed in `possibleSimulators()` (`src/simulator/simulator_def.hpp`); each one
adds a compiled module and compile time.
