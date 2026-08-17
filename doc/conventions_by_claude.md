# PHARE Conventions

### Reference document for the code base

Process conventions (issues, pull requests, review) are in [CONTRIBUTING.md](../CONTRIBUTING.md) and
are not repeated here.

# Sections

1. C++
1. Python
1. CMake
1. Tests
1. Etc

<br/>

# 1. C++

## 1.1 General

The language standard, warning set and compiler flags are set by CMake, not per file. `devMode=ON`
builds with `-Wall -Wextra -pedantic -Werror` (`-Wno-unused-variable -Wno-unused-parameter`
excepted, see `res/cmake/def.cmake`). CI builds with `devMode=ON`, so code must compile warning free
under it.

Formatting is entirely delegated to [`.clang-format`](../.clang-format).

### 1.1.1 Formatting enforcement — open point

`clang-format` is committed but not enforced: 20 of 241 headers under `src/` currently fail
`clang-format --dry-run -Werror`. There is no CI job checking it.

Options:

- **A** — add a CI job running `clang-format --dry-run -Werror` on changed files only. Cheap, no
  repo-wide diff, blocks new drift.
- **B** — one-off repo-wide `clang-format -i`, then CI on all files. Clean state, but a large
  mechanical commit that touches history/blame broadly.
- **C** — add a `pre-commit` hook and leave CI as is. Catches drift earlier, but only for
  contributors who install the hook.

## 1.2 Files

- Headers use the `.hpp` extension.
- Implementation files use the `.cpp` extension.
- File names are `snake_case`.
- Directory layout under `src/` follows the library split: `core/`, `amr/`, `initializer/`,
  `diagnostic/`, `restarts/`, `simulator/`, `python3/`, `phare/`. A header's namespace matches its
  directory (see 1.4), with two established exceptions: `simulator/` lives directly in `PHARE`, and
  `python3/` in `PHARE::pydata`.

### 1.2.1 Include guards — open point

| Form                          | Example                                          | Count |
| ----------------------------- | ------------------------------------------------ | ----- |
| `PHARE_<FILE>_HPP`            | `PHARE_TAGGER_HPP`                               | 202   |
| `<FILE>_HPP`, no prefix       | `TAGGER_STRATEGY_HPP`                            | 32    |
| `PHARE_SRC_<PATH>_<FILE>_HPP` | `PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_DATA_HPP` | 12    |

Options:

- **A** — `PHARE_<PATH_FROM_SRC>_<FILE>_HPP` everywhere. Collision-proof, verbose; a scripted rename
  touches ~234 files but is mechanical and reviewable.
- **B** — `PHARE_<FILE>_HPP` everywhere. Matches the existing majority, cheapest (44 files to
  change). Risk: two headers sharing a basename in different directories would collide silently.

## 1.3 Includes

- Project headers use quotes and a path from `src/`: `#include "core/data/grid/gridlayout.hpp"` (861
  of 981 quoted includes). Relative paths (`"../"`, `"./"`) are not used (0 occurrences).
- Standard library and third-party headers use angle brackets (590 occurrences).
- `SortIncludes` is off, so include order is manual. The prevailing order is project headers first,
  then third-party, then standard library — but this is not enforced and not uniform enough to be a
  rule.

## 1.4 Namespaces

All code lives under the `PHARE` namespace, with a sub-namespace per library: `core`, `amr`,
`solver`, `diagnostic`, `restarts`, `initializer`, `hdf5`, `pydata`. Implementation details that
must be visible in a header go in a nested `detail` namespace.

### 1.4.1 Namespace form — open point

Two forms coexist:

| Form                                                | Count                                                          |
| --------------------------------------------------- | -------------------------------------------------------------- |
| Nested blocks: `namespace PHARE { namespace core {` | 82                                                             |
| Qualified: `namespace PHARE::core`                  | 163 (`PHARE::core` 55, `PHARE::amr` 55, `PHARE::solver` 17, …) |

Note that `NamespaceIndentation: Inner` in `.clang-format` means the nested form costs an
indentation level for the whole file body, while the qualified form does not.

Options:

- **A** — qualified `namespace PHARE::core` everywhere (C++17). Already the majority, removes one
  indent level. Converting a file is a whole-file reindent, so the diff is large per file even
  though the change is trivial.
- **B** — keep nested blocks everywhere. Requires converting the majority; no upside identified.
- **C** — convert opportunistically: any file otherwise touched by a PR moves to the qualified form.
  No dedicated commit, but the split lingers for a long time and PR diffs get noisier.

## 1.5 Naming

Consistent:

- **Types** (classes, structs, enums, aliases-as-types) are `PascalCase`. The only lowercase type
  names are trait / metafunction helpers behaving like standard library traits
  (`is_field_or_tensorfield`, `box_iterator`, `get_value_type`) — this mirrors `std::` and is
  accepted.
- **Private data members and private member functions** end with a trailing underscore: `layout_`,
  `advancePosition_()`. ~1146 sites.
- A **leading** underscore marks a constructor parameter that would otherwise shadow the member it
  initializes: `Box(std::array<Type, dim> _lower, …) : Super{core::Point{_lower}, …}`. Used in 33
  files. It is legal on function parameters (only `_Capital` and `__double` forms are reserved), but
  do not use it on any other kind of identifier.
- **`enum class` only** (26 uses); plain unscoped `enum` is not used.
- **`typedef` is never used** (0 uses) — always `using`.
- `struct` is used for aggregates and trait types whose members are all public; `class` is used as
  soon as there is private state or an invariant to protect. Only 2 structs in the whole tree
  declare `private:` / `protected:`.

### 1.5.1 Function naming — open point

Both cases are widespread:

Out of 448 distinct function names declared in `src/`, 258 are `camelCase` and 169 `snake_case`. By
directory (declaration sites):

| Case         | `core` | `amr` | `diagnostic` | other                                               |
| ------------ | ------ | ----- | ------------ | --------------------------------------------------- |
| `camelCase`  | 104    | 268   | 66           | `python3` 35, `restarts` 4, `hdf5` 1                |
| `snake_case` | 73     | 13    | 10           | `hdf5` 11, `simulator` 9, `python3` 2, `restarts` 1 |

The rough — but not clean — split is: member functions of simulation classes are `camelCase`
(`advanceLevel`, `fillMomentsGhosts`), while free utility functions and standard-library-like
helpers are `snake_case` (`as_unsigned`, `append_to`, `amr_lcl_idx`). `core` violates this in both
directions.

Options:

- **A** — codify the observed split: `camelCase` for member functions, `snake_case` for free
  functions and anything mimicking a `std::` interface. No mass rename; only the outliers move.
  Boundary is fuzzy for static member helpers.
- **B** — `camelCase` for all functions, member or free. One rule, no judgement call; renames ~169
  functions and reads oddly next to `std::` interfaces.
- **C** — `snake_case` for all functions. Matches the standard library and the Python side; renames
  ~258 functions, the largest churn of the three.

### 1.5.2 Template parameter naming — open point

Template parameters mixing several suffix conventions for the same kind of parameter:

- `Field` (67), `FieldT` (25), `Field_t` (16)
- `GridLayout` (53), `GridLayoutT` (21)
- non-type parameters are lowercase: `dim` (81), `direction` (51), `rank` (19), `interp_order`

Options:

- **A** — bare `PascalCase` for type parameters (`Field`, `GridLayout`), the current majority; use a
  `_t` suffix only when the parameter name would shadow a type in scope. Documents existing
  practice; the shadowing exception has to be judged case by case.
- **B** — always suffix type parameters with `_t` (`Field_t`). Never shadows, but renames most
  template headers.
- **C** — always suffix with `T` (`FieldT`). Same benefit as B, smaller visual weight, same churn.

### 1.5.3 Type alias suffixes — open point

`using` aliases split three ways (828 total):

| Form                                       | Count |
| ------------------------------------------ | ----- |
| bare `PascalCase` (`using GridLayout = …`) | 406   |
| `_t` suffix (`using Field_t = …`)          | 240   |
| `_type` suffix (`using value_type = …`)    | 106   |

Part of the `_type` group is required, not a choice: names like `value_type` and `size_type` are
standard-library interface points and must keep that spelling.

Options:

- **A** — bare `PascalCase` for PHARE aliases, `_type` reserved for standard library interface
  points (`value_type`, `size_type`, …). Follows the majority.
- **B** — `_t` for every PHARE alias, `_type` reserved as in A. Makes aliases visually distinct from
  real class names, at the cost of renaming ~406 aliases.

### 1.5.4 `dim` / `interp_order` spelling — open point

Static constexpr members and template parameters for the same two quantities are spelled
inconsistently: `dimension` (60) vs `dim` (5) for the dimensionality, and `interp_order` (9) vs
`interpOrder` (3) for the interpolation order.

Options:

- **A** — `dimension` and `interp_order` for static constexpr members exposed as part of a type's
  interface; `dim` and `interp_order` for template parameters. This is what most of the code does
  today.
- **B** — `dim` and `interp_order` everywhere; shorter, one spelling per concept, renames ~60
  interface members that external code may name.

## 1.6 `NO_DISCARD` and other macros

- Use the `NO_DISCARD` macro from `core/def.hpp`, not `[[nodiscard]]` directly (680 vs 1). Place it
  before the return type.
- `PHARE_DEBUG_DO(...)` wraps code that must only exist in debug builds; it expands to nothing under
  `NDEBUG` unless `PHARE_FORCE_DEBUG_DO` is set.
- New macros go in `core/def.hpp` and are prefixed `PHARE_`; macro helpers not meant for direct use
  are prefixed `_PHARE_`.

## 1.7 Errors and assertions

- **`throw std::runtime_error`** (180 uses) for conditions that can occur at runtime from
  configuration or data: bad input, unmet precondition coming from outside the code, unsupported
  combination.
- **`assert`** (123 uses) for internal invariants that a correct program cannot violate. Asserts are
  compiled out in release builds, so never use one to validate user input.
- Never call `core::Errors::instance()` from code compiled into more than one shared library — the
  singleton is duplicated per DSO. Query global error state through `mpi::any_errors()`.

## 1.8 Return types

Trailing return types (`auto f() -> T`) are effectively unused (2 occurrences). Write the return
type in the leading position, or use `auto` where the type is obvious from the returned expression.

<br/>

# 2. Python

## 2.1 General

- `snake_case` for functions and variables (457 of 538 function definitions), `PascalCase` for
  classes, `snake_case` for module file names (no camelCase module names in `pyphare/`).
- Double quotes for strings (6583 vs 253 single).
- A leading underscore marks module-private names (`_simulator.py`, `_compute_divB`). Note the
  difference with C++, where privateness is marked by a *trailing* underscore (see 1.5).

### 2.1.1 camelCase in the Python API — open point

81 of 538 function definitions are `camelCase`. They are not scattered randomly: they are mostly (a)
the `pharein` configuration DSL and `pyphare.cpp` bindings, which mirror C++ names (`clearDict`,
`changeCentering`, `allocSize`), and (b) AMR helpers named after C++ counterparts (`AMRToLocal`,
`AMRBoxToLocal`).

Options:

- **A** — accept camelCase where a Python name deliberately mirrors a C++ symbol, `snake_case`
  everywhere else. No renames, documents the intent; the boundary needs judgement.
- **B** — `snake_case` throughout, with deprecation aliases for the public `pharein` names. PEP 8
  clean, but `pharein` names appear in every user simulation script and in tests, so the blast
  radius is large.
- **C** — `snake_case` for new code only, leave existing camelCase untouched and undocumented as a
  target. Zero cost, split persists indefinitely.

### 2.1.2 Formatting — open point

The Python code is close to `black` output but not enforced: 8 of 59 files under `pyphare/` would be
reformatted by `black`, and the repo has no `pyproject.toml`, `setup.cfg`, `.pre-commit-config.yaml`
or CI formatting job.

Options:

- **A** — adopt `black` explicitly: add a `pyproject.toml`, run it once over `pyphare/` and
  `tests/`, add a CI check. Matches the de facto style; the one-off commit touches 8 files.
- **B** — adopt `ruff format` plus `ruff check` for lint as well. Same formatting result, also
  catches unused imports and shadowing; one more tool to pin.
- **C** — document the style in prose and leave it unenforced. No new tooling, drift continues.

## 2.2 dependencies and imports

Third party dependencies are stated in the file `requirements.txt` in the project root. Fewer
dependencies is generally better but there should be a cost/benefit assessment for adding new
dependencies.

### 2.2.1 Python file import structure.

Generally, we want to avoid importing any dependency at the top of a python script that may rely on
binary libraries.

Exceptions to this are things like numpy, which are widely used and tested.

Things to expressly avoid importing at the top of a python script are

- h5py
- mpi4py
- scipy.optimize

The first two are noted as they can, and will pull in system libraries such as libmpi.so and
libhdf5.so, which may not be the libraries which were used during PHARE build time, this can cause
issues at runtime.

scipy.optimize relies on system libraries which may not be available at runtime.

The gist here is to only import these libraries at function scope when you actually need them, so
that python files can be imported or scanned for tests and not cause issues during these operations,
until the functions are used at least.

## 2.3 Simulation state

`pyphare.pharein` is stateful: the current simulation lives in `global_vars.sim`. Code that builds a
simulation must reset it (`ph.global_vars.sim = None`) before constructing a new one, and release
C++ resources afterwards via `ph.clearDict()` or `Simulator.reset()`. This applies to tests and to
any script that builds more than one simulation in a process.

<br/>

# 3. CMake

## 3.1 General

- Commands are lowercase: `add_subdirectory(...)`, `set(...)`, `if(...)`. Four legacy uppercase
  calls remain (`SET(`, `ENDIF(`) and should be lowercased when the surrounding code is touched.
- Closing keywords repeat the opening name: `endif(...)`, `endfunction(name)`.
- 2-space indent inside blocks.
- Options are declared in `res/cmake/options.cmake` with `option()` and a default; reusable
  functions in `res/cmake/def.cmake`; test registration in `res/cmake/test.cmake`. Feature-specific
  logic goes in its own `res/cmake/<feature>.cmake` (`coverage.cmake`, `cppcheck.cmake`,
  `bench.cmake`).
- Internal helper functions end with a trailing underscore, as in C++: `add_phare_test_`,
  `set_exe_paths_`, `phare_sanitize_`. Functions without the underscore are the public interface for
  `CMakeLists.txt` files.

## 3.2 Simulation permutations

Supported combinations of time integrator, reconstruction, slope limiter and Riemann solver are
listed in `res/sim/all.txt`. `src/python3/CMakeLists.txt` generates one `pybindlibs/cpp_<id>` module
per line. Adding a permutation means adding a line there, not editing generated code; at runtime
`pyphare/pyphare/cpp/__init__.py` picks the module via `simulator_id(sim)`.

<br/>

# 4. Tests

## 4.1 General

- `tests/` mirrors `src/`: `tests/core/`, `tests/amr/`, `tests/diagnostic/`, `tests/initializer/`,
  `tests/simulator/`, plus `tests/functional/` for end-to-end physics cases.
- C++ tests use GoogleTest. `TEST` (166), `TYPED_TEST` (118) and `TEST_F` (74) are all in normal
  use; `TYPED_TEST_SUITE` (61) is the standard way to cover the `dim` × interpolation-order matrix.
- The C++ test source of a directory is `test_main.cpp`, holding `::testing::InitGoogleTest` and
  `RUN_ALL_TESTS`. Tests needing SAMRAI construct a `PHARE::SamraiLifeCycle` in `main`.
- Python test files are named `test_*.py`, test classes `PascalCase` ending in `Test`
  (`AdvanceTest1D`, `HarrisTest`).
- Embedded-Python C++ tests (`simulator`, `multiphysics`) need `PYTHONPATH` set and must be run
  through `ctest`, not by invoking the binary directly.

## 4.2 Registration

Always register tests with the macros in `res/cmake/def.cmake` — never a raw `add_test()`:

| Macro                                                       | Use                                               |
| ----------------------------------------------------------- | ------------------------------------------------- |
| `add_phare_test(binary directory)`                          | C++ test; runs under `mpirun` when `-DtestMPI=ON` |
| `add_no_mpi_phare_test(binary directory)`                   | C++ test that must stay serial                    |
| `add_python3_test(name file directory)`                     | Python test; `mpirun` when `-DtestMPI=ON`         |
| `add_no_mpi_python3_test(name file directory)`              | Python test that must stay serial                 |
| `add_mpi_python3_test(N name file directory)`               | Python test on `N` ranks                          |
| `phare_exec(level target exe directory)`                    | C++ test gated on execution level                 |
| `phare_python3_exec(level target file directory ...)`       | Python test gated on execution level              |
| `phare_mpi_python3_exec(level N target file directory ...)` | both, gated                                       |

Notes:

- Python tests appear in CTest as `py3_<name>`, and as `py3_<name>_mpi_n_<N>` when `N > 1`.
- The non-levelled macros register at level 1 and are therefore skipped when
  `PHARE_EXEC_LEVEL_MIN > 1`.
- Do not pass CMake variables as arguments to `phare_python3_exec` / `phare_mpi_python3_exec` when
  the file runs Python unit tests — they interfere with the test runner.

### 4.2.1 CMake test project names

Each test directory has its own `CMakeLists.txt` with `project(test-<something>)` — dashes, 70 of 70
occurrences — while the source file is `test_<something>.cpp` with an underscore. Both are
consistent within their own domain; keep the two spellings.

The exceptions are `tests/diagnostic/test-diagnostics_{1,2,3}d.cpp`, the only three C++ test sources
using a dash. Rename them to `test_diagnostics_<n>d.cpp` when that directory is next touched.

### 4.2.2 Execution levels — open point

Levels are compared against `PHARE_EXEC_LEVEL_MIN` (default 1) and `PHARE_EXEC_LEVEL_MAX` (default
10), so a test registered above 10 is excluded from a default build. The values actually in use are
9 (25 tests), 11 (12), 21 (7), 101 (1) and 121 (1). What each value means is not written down
anywhere; the conventional reading is "9 and below = default suite, 11+ = heavy, 101+ = very heavy
3D/MHD".

Options:

- **A** — document the tiers as they are (default ≤ 10, heavy 11–100, very heavy 101+) in
  `res/cmake/def.cmake` and stop introducing new values. No code change.
- **B** — define named CMake constants (e.g. `PHARE_LEVEL_DEFAULT`, `PHARE_LEVEL_HEAVY`,
  `PHARE_LEVEL_VERY_HEAVY`) and use them at every call site. Self-documenting, ~46 call sites to
  update.

## 4.3 Heavy tests

Tests above the default level are excluded unless `-DPHARE_EXEC_LEVEL_MAX` is raised (e.g. `101`).
Many simulator and functional Python tests are additionally guarded by `if(HighFive)` because they
need HDF5 diagnostics.

`lowResourceTests=ON` is a second, independent exclusion switch, but it is only consulted in one
place (`pyphare/pyphare_tests/test_pharesee/CMakeLists.txt`) and CI does not set it — CI configures
with `devMode=ON -Dbench=ON -Dphare_configurator=ON -DPHARE_PYTEST_SIMULATORS=1` and leaves the
execution level at its default. New heavy tests should be gated by execution level rather than by
`lowResourceTests`.

<br/>

# 5. Etc

## 5.1 Static analysis

`cppcheck` is available as a build target:

```bash
cmake -S . -B build -Dcppcheck=ON
cmake --build build --target cppcheck-xml   # or cppcheck-html
```

CodeQL runs in CI (`.github/workflows/codeql.yml`).

## 5.2 Logging and run artifacts

- Rank-specific logs go to `.log/`.
- Runtime metadata goes to `.phare/`.
- Neither directory is committed.

## 5.3 Documentation

Sphinx sources live in `doc/source/`. Build with:

```bash
uv run python -m pip install -r doc/source/requirements.txt
uv run make -C doc html
```

Doxygen is configured through `doc/Doxyfile.in`. There is no rule yet on which C++ entities must
carry Doxygen comments — the code base has very few of them.

## 5.4 Open points summary

| #   | Topic                                                 | Section |
| --- | ----------------------------------------------------- | ------- |
| 1   | C++ formatting not enforced                           | 1.1.1   |
| 2   | Include guard naming                                  | 1.2.1   |
| 3   | Nested vs qualified namespaces                        | 1.4.1   |
| 4   | camelCase vs snake_case functions                     | 1.5.1   |
| 5   | Template parameter suffixes                           | 1.5.2   |
| 6   | Type alias suffixes                                   | 1.5.3   |
| 7   | `dim` vs `dimension`, `interp_order` vs `interpOrder` | 1.5.4   |
| 8   | camelCase in the Python API                           | 2.1.1   |
| 9   | Python formatting not enforced                        | 2.1.2   |
| 10  | Execution level semantics                             | 4.2.2   |

Each open point should become its own issue, decided independently. Until an option is chosen, match
the file you are editing.
