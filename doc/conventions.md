# PHARE Conventions

### Reference document for the code base

# Sections

1. [C++](#c)
1. [Python](#python)
1. [CMake](#cmake)
1. [Tests](#tests)
1. [Etc](#etc)

<br/>

# C++

## Project architecture

C++ sources live in `src/`, one subfolder per component. Each component is a CMake library named
`phare_<component>` (except `hdf5/`, which is headers only) and a namespace (see
[Namespaces](#namespaces)).

Components are layered: a component includes headers of the components below it, never of the
components above it. From top to bottom:

- `python3/`: pybind11 modules, one per simulation permutation, plus `cpp_etc`; `phare/`: the
  standalone executable.
- `simulator/`: ties together the model, the solvers, the diagnostics and the restarts.
- `diagnostic/`, `restarts/`: output and checkpoints, using `core/`, `mpi/` and `hdf5/` (and `amr/`
  for `diagnostic/`).
- `amr/`: the SAMRAI-based layer: hierarchy, refinement, transfer of data between levels
  (messengers), solvers, tagging. Uses `core/`, `mpi/` and SAMRAI.
- `mpi/`, `hdf5/`: thin wrappers around MPI and HDF5, using `core/` utilities.
- `core/`: physics and numerics on a single patch: grids, fields, particles, ion and electron state,
  numerical operators. `core/` never includes `amr/` or SAMRAI.

`initializer/` holds the input dictionary (`PHAREDict`) built from the Python configuration, read by
`core/`, `amr/` and most components above them.

Code that needs patches, levels or SAMRAI goes in `amr/`; code that works on a single patch goes in
`core/`.

**Open point**: two includes go against the layering:

- `core/` and `initializer/` include each other: 14 `core/` headers include
  `initializer/data_provider.hpp`, which includes `core/def.hpp` and `core/utilities/span.hpp`.
- `initializer/dictator.cpp` includes `python3/pybind_def.hpp`.

## Header files

### Self-contained headers

A header includes everything it uses, so that it compiles on its own. A header that needs a
build-time macro fails with `#error` when the macro is missing (`diagnostic/diagnostics.hpp`,
`restarts/restarts.hpp`).

### `#define` guard

Headers use `#ifndef` guards; `#pragma once` is not used.

**Open point**: guard names follow several forms, and 30 no longer match their file
(`amr/messengers/communicator.hpp` uses `PHARE_QUANTITY_REFINER_HPP`). Out of 236 guarded headers:

| Form                               | Count |
| ---------------------------------- | ----- |
| `PHARE_<FILE>_HPP`                 | 76    |
| `PHARE_<PATH_FROM_SRC>_<FILE>_HPP` | 50    |
| `PHARE_<PARTIAL_PATH>_<FILE>_HPP`  | 41    |
| `PHARE_SRC_<PATH>_<FILE>_HPP`      | 12    |
| no `PHARE_` prefix                 | 27    |
| name of another file               | 30    |

- **A** — `PHARE_<PATH_FROM_SRC>_<FILE>_HPP` everywhere. Unique by construction and scriptable;
  renames about 186 guards.
- **B** — `PHARE_<FILE>_HPP` everywhere. Shortest, but not unique: `diagnostic/detail/h5writer.hpp`
  and `restarts/detail/h5writer.hpp` would get the same guard.
- **C** — any unique `PHARE_`-prefixed name. No renaming, but stale names stay.

### Include what you use

A file should include a header if and only if it uses a symbol of the header. This means that:

- if a file does not use any symbols defined in a header, this header should not be included.
- you should not rely on transitive inclusions: include directly the symbol you are using.

`mpi.h` is only included through `phare_mpi.hpp`. The umbrella headers `phare_core.hpp`,
`phare_amr.hpp` and `phare_solver.hpp` gather types for tests and bindings.

### Forward declarations

Include the header that declares a type rather than forward-declaring it.

### Defining functions in header files

The code is header-only by default: templates and most classes are defined in headers, and
`diagnostic/`, `restarts/` and `simulator/` are header-only libraries. A non-template function
defined at namespace scope in a header is `inline` or `constexpr` (the one exception,
`build_config()` in `core/def/phare_config.hpp`, is a latent multiple-definition error).

The few `.cpp` files (13) hold what cannot live in a header: SAMRAI and MPI set-up, singletons,
Python module definitions and `main`.

### Names of includes

`#include` statements for PHARE headers should refer to the full path with respect to the `src/`
folder, like:

```cpp
#include "core/data/grid/grid.hpp"
```

Do not use relative paths. Standard library headers use angle brackets.

**Open point**:

- 74 includes name a header of the same directory without its path
  (`#include "resources_manager.hpp"`), and 3 PHARE headers are included with angle brackets.
  - **A** — always the full path from `src/`, in quotes (77 edits).
  - **B** — allow same-directory names.
- SAMRAI headers are included as `<SAMRAI/...>` 158 times and as `"SAMRAI/..."` 21 times (9 files
  mix both).
  - **A** — angle brackets, the majority.
  - **B** — quotes.

### Order of includes

**Open point**: `.clang-format` does not sort includes (`SortIncludes: false`). In the 183 files
that include more than one kind of header, the blocks come in this order:

| Order                | Files |
| -------------------- | ----- |
| PHARE, std           | 81    |
| PHARE, SAMRAI, std   | 36    |
| std, PHARE           | 27    |
| PHARE, SAMRAI        | 10    |
| SAMRAI before PHARE  | 5     |
| other or interleaved | 24    |

Proposal, which matches the majority:

- start by phare's headers, grouped by lead namespace (`amr`, `core`, etc ...). Each group appear in
  alphabetical order. Inside a group, order by alphabetical order
- then SAMRAI's headers, sorted by alphabetical order.
- finally system headers, ordered by alphabetical order.

Once chosen, the order could be enforced by clang-format (`SortIncludes` with `IncludeCategories`).

## Scoping

### Namespaces

All code is in the `PHARE` namespace. A file's namespace matches the name of its top-level folder in
`src/` (see [Project architecture](#project-architecture)), with these established exceptions:

- files in `simulator/` and `phare/` live directly in `PHARE`;
- files in `python3/` live in `PHARE::pydata`;
- `amr/solvers/`, `amr/physical_models/`, `amr/level_initializer/`, `amr/messenger_registration.hpp`
  and `amr/multiphysics_integrator.hpp` live in `PHARE::solver`;
- a few files live directly in `PHARE`: `amr/samrai.hpp`, `core/logger.hpp`,
  `core/def/phare_config.hpp`, `core/models/options/`.

Implementation details that must be visible in a header go in a nested `detail` namespace.

**Open point**: two syntaxes coexist. Per file: `namespace PHARE::core` 158, nested
`namespace PHARE { namespace core {` 65, both 9. `diagnostic/`, `python3/`, `restarts/`, `mpi/` and
`hdf5/` only use the first. With `NamespaceIndentation: Inner`, the nested form indents the whole
file.

- **A** — `namespace PHARE::core` everywhere. Already the majority; converting a file re-indents all
  of it.
- **B** — either form, one per file. Only the 9 mixed files change.

**Open point**: 7 headers have `using namespace` at namespace scope (`using namespace PHARE::core`
in `amr/data/particles/refine/split_{1,2,3}d.hpp`, `using namespace hdf5::h5` in 4
`diagnostic/detail/` writers), and 2 declare `namespace py = pybind11;` at global scope
(`initializer/python_data_provider.hpp`, `python3/cpp_simulator.hpp`). None has `using namespace` at
global scope.

- **A** — no `using namespace` in headers, as in Google's guide (7 edits).
- **B** — allowed inside `PHARE` namespaces only.

### Internal linkage

**Open point**: an unnamed namespace or a `static` variable at namespace scope gives every
translation unit its own copy. Headers have one unnamed namespace (`simulator/simulator.hpp`) and a
few `static` variables (`to_primitive_converter.hpp`, `hdf5/detail/h5/h5_file.hpp`, and some
`static constexpr` ones).

- **A** — neither in headers: use `inline` or `inline constexpr` variables (about 4 edits).
- **B** — allow `static constexpr`, forbid the rest.

### Static and global variables

There are no mutable variables at namespace scope. State shared across the program is a
function-local static returned by an accessor (`core::Errors::instance()`,
`PHAREDictHandler::INSTANCE()`) or a `static inline` class member. A function-local static defined
in a header exists once per shared library that uses it: see [Exceptions](#exceptions).

## Classes

### Doing work in constructors

Constructors may throw on invalid input or configuration; many take a `PHAREDict` and read it. They
never call virtual functions. Initialization in two steps is reserved for what SAMRAI imposes: views
whose data buffers are set later by the `ResourcesManager` (`setBuffer`, `isUsable()`), and initial
conditions set by `initialize()`.

### Implicit conversions

Implicit conversions are deliberate between PHARE and SAMRAI boxes (`amr/utilities/box/amr_box.hpp`)
and from `std::array` to `Point`.

**Open point**: 18 of 104 single-argument constructors are `explicit`, and counterparts disagree:
`SolverPPC` and `MHDModel` are `explicit`, `SolverMHD` and `HybridModel` are not. 14 `explicit` are
on constructors that need two or more arguments.

- **A** — `explicit` on every single-argument constructor and conversion operator, except the
  deliberate conversions above (Google's rule). About 75 constructors to annotate.
- **B** — `explicit` required on constructors taking a `PHAREDict` or a scalar. About 36
  constructors.

### Copyable and movable types

Most types (410 of 439) declare no copy or move operation and rely on the ones the compiler
generates. A type that must not be copied deletes its copy operations.

**Open point**: about 10 types declare part of the set only. `GridLayout` and `Ions` default their
copy and move constructors, which leaves them without copy assignment; 3 types delete copy without
declaring move, so they cannot be moved either. `LoadBalancerManager` can be copied although its
destructor unregisters a SAMRAI variable.

- **A** — a type that declares any copy or move operation, or a destructor with side effects,
  declares all four copy and move operations (`= default` or `= delete`).
- **B** — no rule beyond relying on the generated operations.

### Structs vs. classes

**Open point**: Google uses `struct` for passive data only. In PHARE, 205 of 209 structs have only
public members, but most of them carry behaviour: 78 are traits or tags, 62 hold only functions, 48
hold data and functions. 41 classes have only public members.

- **A** — `struct` when every member is public, `class` as soon as one is private or protected. This
  is what the code does; 4 structs with private members (`CompositeFill`, `MinimizingVector`,
  `Reconstructor`, `Range`) and the all-public classes change.
- **B** — Google's rule: `struct` for passive data, stateless traits and functors; everything else
  is a `class`. The 48 structs holding data and functions change.

### Structs vs. pairs and tuples

`std::tuple` holds heterogeneous compile-time lists (`std::tuple<Patterns...>`). Functions returning
several values return a pair or a tuple, which callers unpack with structured bindings.

**Open point**: a few pairs and tuples are records accessed by position: `std::get<0>` to
`std::get<4>` on particle arrays from Python (`python3/particles.hpp`),
`core/data/particles/particle_packer.hpp`, and `std::pair<std::string, Attributes>` in the virtual
interface of diagnostic writers (`diagnostic/detail/h5typewriter.hpp`).

- **A** — pairs and tuples for compile-time lists and for returned values unpacked at once; a named
  struct when members are accessed by position or cross an interface. About 3 places change.
- **B** — named structs everywhere (Google's rule). About 55 functions returning
  `std::make_pair`/`std::make_tuple` change.

### Inheritance

- An overriding function is marked `override`, or `final` (as on SAMRAI `PatchData`,
  `PatchDataFactory`, `BoxGeometry` overrides). It is never also marked `virtual` (4 functions are).
- A polymorphic base class has a public virtual destructor.

**Open point**: 14 overrides lack `override` today, 9 of them in
`amr/messengers/field_operate_transaction.hpp`.

- **A** — fix them and add `-Wsuggest-override` to the `devMode` flags, so that the compiler
  enforces the rule.
- **B** — fix them only.

### Operator overloading

Operators keep their usual meaning: `operator()` for element access and functors, `operator[]`,
arithmetic on `Point` and `Box`. There are no user-defined literals and no overloads of `,`, `&&` or
`||`.

**Open point**: a few overloads have an unexpected meaning: unary `Grid::operator&` returns the
underlying field (`core/data/grid/grid.hpp`), `operator>>` copies the border of a `MaskedView`
(`core/data/ndarray/ndarray_vector.hpp`), `Box::operator*` is the intersection and returns an
optional, `Field::operator=(Field*)` rebinds the data buffer.

- **A** — overload operators only with their usual meaning; these become named functions.
- **B** — keep these, and add no new ones.

### Access control

**Open point**: Google makes all class data members private. In PHARE, 32 classes have public data
members (133 members: physics states such as `MHDState`, the `*MessengerInfo` classes, the
`FieldData` and `ParticlesData` payloads, Riemann solver scratch next to private members), and 13
have protected data members (32 members, e.g. `Interpolator`).

- **A** — public data only in types whose members are all public, which then become structs; private
  data elsewhere. 10 to 15 classes change.
- **B** — public data allowed in state and payload types; protected data forbidden.

### Declaration order

Group similar declarations together, placing public parts earlier.

A class definition should usually start with a `public:` section, followed by `protected:`, then
`private:`. Omit sections that would be empty.

A class may open with a private block of type aliases and static constants used by its public
interface (93 classes do). 21 classes put other private or protected members before public ones.

Within each section, prefer grouping similar kinds of declarations together, and prefer the
following order:

1. Types and type aliases (typedef, using, enum, nested structs and classes, and friend types)
1. (Optionally, for structs only) non-static data members
1. Static constants
1. Factory functions
1. Constructors and assignment operators
1. Destructor
1. All other functions (static and non-static member functions, and friend functions)
1. All other data members (static and non-static)

## Functions

### Inputs and outputs

Prefer returning values; return several values as a pair, a tuple or a struct (see
[Structs vs. pairs and tuples](#structs-vs-pairs-and-tuples)). Non-const reference parameters are
for data filled in place (fields) and objects modified in place; a parameter that is not modified is
a `const&`. Output parameters are references, never pointers.

**Open point**: position of output parameters. In a sample of 16 functions with a pure output
parameter, 10 put it after the inputs (`Faraday(B, E, Bnew, dt)`, `Ampere(B, J)`) and 5 first
(`linear_time_interpolate(dst, ...)`, `Interpolator::operator()(Field&, ...)`). Objects modified in
place come first (`updatePopulations(ions, ...)`).

- **A** — inputs, then outputs (Google's rule, and most of `core/` numerics); objects modified in
  place first.
- **B** — destination first, as `memcpy` and SAMRAI do.

### Default arguments

Default arguments are allowed on non-virtual functions. A virtual function has default arguments
only to repeat those of an external base class exactly (the 4 cases override SAMRAI functions).

### Trailing return type syntax

Write the return type before the function name, or `auto`. Trailing return types are used only where
needed: return types that depend on the parameters (`-> decltype(...)`) and deduction guides.

## Other C++ features

### Ownership and smart pointers

- Never `new` or `delete` directly.
- `std::unique_ptr` for objects with a single owner, created with `std::make_unique`.
- `std::shared_ptr` where SAMRAI or pybind11 require it (most of the 273 uses), and for objects
  shared on purpose such as the `ResourcesManager`. Create them with `std::make_shared`.

### Exceptions

Unlike Google's guide, PHARE uses exceptions.

- throw a **`std::runtime_error`** with a meaningful error message for conditions that can occur at
  runtime from configuration or data: bad input, unmet precondition coming from outside the code,
  unsupported combination.
- a failure during a time step that all MPI ranks must learn about (a particle leaving the ghost
  box, a solver step failing) throws a `core::DictionaryException`, whose key/value pairs describe
  the failure.
- catch by `const&`. `catch (...)` is only used where C++ returns to Python
  (`Simulator::initialize`, `Simulator::advance`), which rethrow a `std::runtime_error`.
- use **`assert`** for internal invariants that a correct program cannot violate. Asserts are
  compiled out in release builds, so never use one to validate user input.

An error detected on one MPI rank must stop all ranks. Solvers catch `core::DictionaryException` and
pass it to `PHARE_LOG_ERROR` (`mpi/mpi_utils.hpp`), which calls `MPI_Abort` by default. With
`useExceptionsInsteadOfMPIAbort` (forced by `devMode`), `PHARE_LOG_ERROR` records the error in
`core::Errors` instead; the solver then calls `mpi::any_errors()`, which is true on every rank as
soon as one rank recorded an error, and every rank throws up to the simulator. Only reach
`core::Errors` through `PHARE_LOG_ERROR` and `mpi::any_errors()`: `Errors::instance()` is a
function-local static in a header, so every shared library that calls it gets its own copy.

**Open point**: two assertion macros are used: `assert` (67) and SAMRAI's `TBOX_ASSERT` (35, in
`amr/` only).

- **A** — `assert` everywhere.
- **B** — `TBOX_ASSERT` allowed in `amr/`.

### Run-time type information (RTTI)

`dynamic_cast` downcasts SAMRAI abstractions (`PatchData`, `BoxOverlap`, `BoxGeometry`, ..., about
67 casts), to a reference when the type is known so that a mismatch throws. `typeid` is not used.

**Open point**: about 34 casts downcast PHARE's own interfaces: `IPhysicalModel` to `HybridModel` or
`MHDModel` (22), `IMessengerInfo` (6), `IMessenger` (5).

- **A** — `dynamic_cast` to a reference wherever the type is guaranteed (current practice).
- **B** — `static_cast` for PHARE's own interfaces when the type is guaranteed by construction.
- **C** — remove the need for these downcasts with templates or visitors.

### Casting

Use C++ casts, never C-style casts. `static_cast` is the common one; `reinterpret_cast` only views
an array as another array type (4 uses); `const_cast` mostly implements a non-const getter with its
const version.

### Streams

Output uses iostreams; `printf` and `fmt` are not used. `std::cout` is for information messages: the
simulator redirects it to a per-rank file in `.log/`. `PHARE_LOG_LINE_STR` and `PHARE_LOG_LINE_SS`
(`core/logger.hpp`) print debug traces and are compiled out in release builds. Errors are reported
with exceptions, not printed.

### Preincrement and predecrement

Use `++i` rather than `i++` (228 of 255 `for` loops).

### Use of const

Use `const` wherever possible, to the right of what it qualifies: `auto const&`, `T const&`,
`T const*`, `double const x`. clang-format enforces it (`QualifierAlignment: Right`); `const` on the
left only appears where a SAMRAI signature is repeated.

### Use of constexpr

Values known at compile time are `constexpr`, and `if constexpr` is the usual way to branch on
template parameters (439 uses).

**Open point**: the order of `static`, `constexpr` and the type varies: `static constexpr auto` 115,
`auto static constexpr` 80, `constexpr static auto` 23, `auto constexpr static` 11,
`static auto constexpr` 6; for local variables, `auto constexpr` 41 and `constexpr auto` 37. The
same member `dimension` is declared in 5 ways.

- **A** — `static constexpr T`, the most common order and the standard library's.
- **B** — `T static constexpr`, in the spirit of east `const`. About 230 edits.

Either order can be enforced by clang-format (`QualifierAlignment: Custom` with `QualifierOrder`,
replacing the current `Right`), which would also rewrite the existing code.

### Integer types

- Sizes and counts are `std::size_t`. Fixed-width types are written with their namespace:
  `std::uint32_t` (23 bare `size_t` remain).
- AMR indices are signed `int`, as in SAMRAI (`Point<int, dim>`, `Box<int, dim>`, the particle
  `iCell`); local indices into the arrays of a patch are `std::uint32_t`.
- Unlike Google's guide, which prefers `int`, PHARE uses unsigned types for sizes and local indices.

**Open point**: loop indices are `auto` 114 times (35 with a `0u` initializer), `std::size_t` 97,
`int` 13 (MPI ranks), `std::uint16_t` 10, bare `size_t` 7.

- **A** — `std::size_t`.
- **B** — `auto i = 0u`.
- **C** — the type of the loop bound.

### Floating-point types

Physics is computed in `double`; `long double` is not used. `float` is used for particle splitting
coefficients and for diagnostics output, which is single precision unless PHARE is built with
`PHARE_DIAG_DOUBLES`.

### Preprocessor macros

- Use the `NO_DISCARD` macro from `core/def.hpp`, not `[[nodiscard]]` directly. Place it first in
  the declaration: `NO_DISCARD auto size() const`.
- `PHARE_DEBUG_DO(...)` wraps code that must only exist in debug builds; it expands to nothing under
  `NDEBUG` unless `PHARE_FORCE_DEBUG_DO` is set.
- New macros go in `core/def.hpp` and are prefixed `PHARE_`; macro helpers not meant for direct use
  are prefixed `_PHARE_`.

**Open point**:

- 7 macros lack the `PHARE_` prefix: `NO_DISCARD`, `DISABLE_WARNING`, `ENABLE_WARNING` and four
  `DIAG_*` helpers (`core/def/pragma_disable.hpp`).
- macros are also defined outside `core/def.hpp`: in `core/def/`, `core/logger.hpp`,
  `mpi/mpi_utils.hpp`, `diagnostic/diagnostics.hpp`, `hdf5/phare_hdf5.hpp` and
  `python3/cpp_simulator.hpp`.
- a name starting with an underscore and a capital letter is reserved by the language, so the
  `_PHARE_` prefix for helpers (`_PHARE_TO_STR`, `_PHARE_WITH_HIGHFIVE`) uses reserved names.

Options:

- **A** — every macro prefixed `PHARE_`, helpers `PHARE_DETAIL_`. `NO_DISCARD` becomes
  `PHARE_NO_DISCARD` (717 lines).
- **B** — as A, but `NO_DISCARD` stays the one unprefixed macro.
- **C** — as A, but `NO_DISCARD` is replaced by the standard `[[nodiscard]]` it expands to (717
  lines).

### 0 and nullptr/NULL

Use `nullptr` (2 `NULL` remain, in `mpi/mpi_utils.cpp`).

### Type deduction (including auto)

`auto` is the default for local variables (about 70% of them), `auto const` when the variable is not
modified. Structured bindings unpack pairs, tuples and structs (186 uses).

### Class template argument deduction

Class template argument deduction is used with `std::array` and with PHARE types that provide
deduction guides (`Point`, `Box`, the transformers in `amr/solvers/`). Add a deduction guide before
relying on it for a new type.

### Lambda expressions

Local lambdas capture by reference with `[&]` (267 of 370 lambdas). `[=]` is not used. A lambda
stored or passed beyond its scope lists its captures.

### Template metaprogramming

Compile-time code is central to PHARE: each simulation permutation of `res/sim/all.txt` is compiled
into its own Python module, and `if constexpr`, fold expressions and `std::index_sequence` are
common. Template parameters are declared with `typename`, not `class` (1105 vs 14).

Concepts are named in `PascalCase` (`Spannable`, `CountedEnum`, `HasInterpOrder`).

**Open point**: templates are constrained both with `std::enable_if` (17, mostly older code) and
with C++20 `requires` clauses and concepts (17, all from 2025 on).

- **A** — `requires` and concepts in new code; existing `std::enable_if` stays until the code around
  it changes.
- **B** — replace the 17 `std::enable_if` now.

### Third-party libraries

Each library is confined to the components that need it:

- SAMRAI: `amr/`, plus `simulator/simulator.hpp`, `diagnostic/diagnostic_model_view.hpp` and
  `phare_amr.hpp`. Never in `core/`.
- pybind11: `python3/`, and `initializer/` to read the Python configuration.
- HighFive: `hdf5/`, the HDF5 writers of `diagnostic/` and `restarts/`, and `python3/cpp_etc.cpp`.
- MPI: the C API is called in `mpi/` (and `hdf5/`), and `mpi.h` is only included through
  `phare_mpi.hpp`.
- cppdict: the `PHAREDict` type, used everywhere.

### Aliases

Use `using`, never `typedef`.

Omit `typename` where C++20 makes it optional: alias declarations (`using X = T::type;`), return
types and member declarations in class scope, `static_cast<T::type>`, trailing return types. Keep it
where it is still required: template arguments (`Foo<typename T::type>`), parameters of
namespace-scope function templates, and expressions (`typename T::type{...}`). 82 aliases still
write `using X = typename ...`.

### Switch statements

Switches are rare (11); values of an `enum` are mostly dispatched with `if` or `if constexpr`. Every
switch over an `enum class` has a `default:` that throws.

**Open point**: a `default:` hides the compiler's check that every enumerator is handled
(`-Wswitch`, an error under `devMode`).

- **A** — no `default:` on a switch over an `enum class`, so that adding an enumerator fails to
  compile until every switch handles it.
- **B** — keep a `default:` that throws (Google requires one).

## Naming

Physics quantities keep their usual symbols: see
[Exceptions to naming rules](#exceptions-to-naming-rules).

### File names

- file extension is `.hpp` for headers, `.cpp` for implementation files.
- file names are `snake_case`.

### Type names

Classes, structs, enums, type aliases and concepts are `PascalCase` (326 of 353 class names). Type
traits and helpers that mimic the standard library follow it: `is_field_or_tensorfield`,
`box_iterator`.

### Variable names

**Open point**: local variables and parameters are mostly `camelCase`. Counting names of several
words (74 of 215 files mix both styles):

|                 | `camelCase` | `snake_case` |
| --------------- | ----------- | ------------ |
| local variables | 1149        | 315          |
| parameters      | 1084        | 188          |

- **A** — `camelCase`, the majority.
- **B** — `snake_case`, as in Google's guide and in the Python code. About 80% of the names of
  several words change.

### Class data members

Private and protected data members end with an underscore: `layout_` (435 of 538; `m_` is not used).
Public data members of structs have no underscore.

A constructor parameter that would shadow the public member it initializes may start with an
underscore: `Box(..., _lower, _upper)`. Never start a name with an underscore followed by a capital
letter, which is reserved by the language (`_B`, `_Splitter`, `_ModelView` exist today).

**Open point**: 173 private member functions end with an underscore (`advancePosition_`), 86 do not.

- **A** — private member functions end with an underscore too, the majority.
- **B** — only data members do, as in Google's guide. About 170 functions are renamed.

### Constant names

**Open point**: Google's `kName` is not used. `static constexpr` members (381) are mostly
`snake_case` or one lowercase word (`dimension` alone: 91), but 66 are `camelCase` (`nbRefinedPart`,
`rootLevelNumber`), 33 end with an underscore, 7 are `UPPER_CASE` (`NO_ROTATE`) and 7 `PascalCase`.
Namespace-scope constants (32) are `snake_case` 20, `camelCase` 7 (`refinementRatio`), `UPPER_CASE`
4 (`MAX_LEVEL_IDX`).

- **A** — constants are named like variables (see [Variable names](#variable-names)).
- **B** — `kName`, Google's rule; nothing uses it.
- **C** — `UPPER_CASE`; 11 uses, and the same style as macros.

### Function names

Functions that implement a standard library interface keep its names (`begin`, `size`,
`emplace_back`).

**Open point**: of 953 distinct function names, 535 are `camelCase` and 190 `snake_case` (most
others are one lowercase word). Member functions: 447 `camelCase`, 125 `snake_case`; free functions:
94 and 66. By component (`camelCase` : `snake_case`): `amr` 253:38, `core` 175:95, `diagnostic`
67:12, `python3` 38:9, `restarts` 7:2, `simulator` 10:12, `hdf5` 4:15, `mpi` 3:10. Some accessors
are `UPPER_CASE` (`INSTANCE()`, `FROM()`).

- **A** — `camelCase` everywhere, except standard library interfaces. About 190 functions renamed.
- **B** — `camelCase` for member functions, `snake_case` for free functions. Free functions are
  split 60:40 today, so both sides change.
- **C** — `snake_case` everywhere, like the standard library and the Python code. About 535
  functions renamed.

### Namespace names

Namespace names are lowercase (`core`, `amr`, `solver`, `detail`, ...), except the top-level `PHARE`
and `PHARE::MHDOpts` (`phare_simulator_options.hpp`).

### Enumerator names

Enumerations are `enum class`; unscoped `enum` is not used.

**Open point**: of 151 enumerators, 72 are `PascalCase`, 17 acronyms (`WENO3`, `TVDRK3`), 17 physics
symbols (`X`, `B`), 30 lowercase, `camelCase` or `snake_case` (`rho`, `rhoVx`, `primal`,
`domain_only`), 15 `Pascal_with_underscore` (`VecFlux_x`). Enumerations used with the `CountedEnum`
concept end with a `count` enumerator.

- **A** — `PascalCase`, with acronyms and physics symbols kept as they are, and `count` where
  `CountedEnum` requires it.
- **B** — `kName`, Google's rule; nothing uses it.

### Template parameter names

Non-type template parameters are lowercase (`dim`, `direction`, `rank`).

**Open point**: type parameters are bare `PascalCase` 870 times out of 1137 (`Field` 65,
`GridLayout` 51), a single capital 93 times (`T`), `XxxT` 67 times (`GridLayoutT` 31, `FieldT` 19)
and `Xxx_t` 59 times (`Field_t` 15).

- **A** — bare `PascalCase`, with a `T` suffix only when the name clashes with a type alias in
  scope.
- **B** — always a `T` suffix (`GridLayoutT`). Never clashes, but renames most template parameters.

### Macro names

Macros are `UPPER_CASE` and prefixed `PHARE_` (see [Preprocessor macros](#preprocessor-macros) for
the exceptions).

### Aliases

**Open point**: of 833 type aliases, 410 are bare `PascalCase`, 234 end with `_t` and 111 with
`_type` (partly standard library interfaces such as `value_type`, which must keep that spelling).
The alias for the grid layout alone is spelled `GridLayout` (30), `gridlayout_type` (12),
`GridLayout_t` (10) and `GridLayoutT` (7).

- **A** — bare `PascalCase`; `_type` only for standard library interfaces. Follows the majority.
- **B** — `_t` for every PHARE alias, `_type` as in A. Tells aliases from class names, but renames
  about 410 aliases.

### `dim` and `interp_order`

**Open point**: the dimension and the interpolation order are spelled differently for members and
for template parameters. `static constexpr` members: `dimension` 91, `dim` 11; `interp_order` 12,
`interpOrder` 3. Non-type template parameters: `dim` 65, `dimension` 29; `interpOrder` 5,
`interp_order` 2.

- **A** — `dimension` and `interp_order` for members, `dim` and `interp_order` for template
  parameters, the majority for each.
- **B** — `dim` and `interp_order` everywhere. Renames about 90 members that other code refers to.

### Exceptions to naming rules

Physics quantities keep their usual symbols, whatever the rule for that kind of name: variables `B`,
`Ex`, `Etot`, enumerators `X`, `rhoVx`, functions `BxToEx`. Acronyms keep their case: `MHD`, `HLLD`,
`WENO3`.

## Comments

### Comment style

Comments are written in English.

**Open point**: `//` is the usual style. Doxygen comments are also used (`/** */` 311, `//!` 30,
`///` 10, with both `@brief` 254 and `\brief` 69), mostly in `amr/`, but no documentation is
generated from them: `doc/Doxyfile.in` is not configured by any CMake target.

- **A** — Doxygen for interface comments, written `/** @brief ... */` with `@` commands only, and
  `doc/Doxyfile.in` wired into the build.
- **B** — `//` comments only; Doxygen is dropped.

### File comments

Files have no license or author header; the license is in `LICENCE.TXT`.

### TODO comments

**Open point**: the 21 TODO comments use 4 spellings (`TODO`, `TODO :`, `TODO:`, `@TODO`), and 6 of
them link an issue.

- **A** — `// TODO(#1234): what remains to do`, always with an issue.
- **B** — free form.

### Commented-out code

**Open point**: commented-out code appears in about 38 files, most of it in one 183-line block of
`core/numerics/riemann_solvers/hlld.hpp`.

- **A** — no commented-out code; git keeps the old versions.
- **B** — short snippets tolerated.

## Formatting

Formatting is entirely delegated to [`.clang-format`](../.clang-format): 100 columns, 4 spaces and
no tabs, braces on their own line after functions, classes and control statements, `const` to the
right of the type.

The CI checks that the C++ lines changed by a pull request follow it, with clang-format
`CLANG_FORMAT_VERSION` from [`gate.yml`](../.github/workflows/gate.yml) (other versions format
differently), and fails otherwise unless the pull request has the `format-ok` label. Format your
changed lines before pushing:

```bash
git clang-format $(git merge-base origin/master HEAD) -- src tests tools
```

Right now [.clang-format](../.clang-format) does not sort header includes. It is expected that you
follow [this order](#order-of-includes).

Code and string literals are ASCII, except two `static_assert` messages in the Yee grid layouts;
comments may use other characters (mathematical symbols appear in 15 files).

**Open point**: `.clang-format` declares `Standard: Cpp11` while the code is C++20, and 16 of the
250 files of `src/` are not formatted (the CI only checks changed lines).

- **A** — set `Standard: Latest` and reformat `src/` once, in a commit of its own.
- **B** — leave both as they are.

<br/>

# Python

## General

...

## Dependencies and imports

Third party dependencies are stated in the file `requirements.txt` in the project root. Fewer
dependencies is generally better but there should be a cost/benefit assessment for adding new
dependencies.

### Python file import structure

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

<br/>

# CMake

## General

...

<br/>

# Tests

## General

...

<br/>

# Etc

## General

...
