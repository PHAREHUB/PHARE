# PHARE Conventions

### Reference document for the code base

# Sections

1. C++
1. Python
1. CMake
1. Tests
1. Etc

<br/>

# 1. C++

## Files

### Project organization

PHARE C++ source files are located in `src/`, and are organized in main subfolders: `core/`, `amr/`,
`initializer/`, `diagnostic/`, `restarts/`, `simulator/`, `python3/`, `phare/`. A files's namespace
matches the name of its parent directory among the ones just listed, with two established
exceptions: files in `simulator/` live directly in `PHARE`, and files in `python3/` live in
`PHARE::pydata`.

### File naming

- file extension is `.hpp` for headers, `.cpp` for implementation files.
- file names are `snake_case`.

### `#define` guard

**Open-point**: Options:

- **A** — `PHARE_<PATH_FROM_SRC>_<FILE>_HPP` everywhere. Collision-proof, verbose; a scripted rename
  touches ~234 files but is mechanical and reviewable.
- **B** — `PHARE_<FILE>_HPP` everywhere. Matches the existing majority, cheapest (44 files to
  change). Risk: two headers sharing a basename in different directories would collide silently.

### Include what you use

A file should include a header if and only if it uses a symbol of the header. This means that:

- if a file does not use any symbols defined in a header, this header should not be included.
- you should not rely on transitive inclusions: include directly the symbol you are using.

### Names of includes

`#include` statements for PHARE headers should refer to the full path with respect to the `\src`
folder, like:

```cpp
#include core/data/grid/grid.hpp
```

Do not use relative paths.

### Order of includes

**Open point** No established convention. One should be chosen, and it could be enforced via
clang-format.

Proposal:

- start by phare's headers, grouped by lead namespace (`amr`, `core`, etc ...). Each group appear in
  alphabetical order. Inside a group, order by alphabetical order
- then SAMRAI's headers, sorted by alphabetical order.
- finally system headers, ordered by alphabetical order.

<!-- ## Scoping -->

<!-- ### Namespace -->

## Classes

### Structs vs. classes

`struct` should be used for passive objects that carry data; everything else is a `class`.

### Declaration order

Group similar declarations together, placing public parts earlier.

A class definition should usually start with a `public:` section, followed by `protected:`, then
`private:`. Omit sections that would be empty.

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

## Other C++ features

### Macros

- Use the `NO_DISCARD` macro from `core/def.hpp`, not `[[nodiscard]]` directly.
- `PHARE_DEBUG_DO(...)` wraps code that must only exist in debug builds; it expands to nothing under
  `NDEBUG` unless `PHARE_FORCE_DEBUG_DO` is set.
- New macros go in `core/def.hpp` and are prefixed `PHARE_`; macro helpers not meant for direct use
  are prefixed `_PHARE_`.

### Errors and assertions

- throw a **`std::runtime_error`** with a meaningful error message for conditions that can occur at
  runtime from configuration or data: bad input, unmet precondition coming from outside the code,
  unsupported combination.
- use **`assert`** for internal invariants that a correct program cannot violate. Asserts are
  compiled out in release builds, so never use one to validate user input.

## Naming

## Comments

## Formatting

Formatting is entirely delegated to [`.clang-format`](../.clang-format).

For now formatting is not enforced by the CI, therefore it is expected that developpers properly
format their files before

Right now [.clang-format](../.clang-format) does not sort header includes. It is expected that you
follow [this order](#order-of-includes).

<br/>

# 2. Python

## 2.1 General

...

## 2.2 dependencies and imports

Third party depenencies are stated in the file `requirements.txt` in the project root. Fewer
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

<br/>

# 3. CMake

## 3.1 General

...

<br/>

# 4. Tests

## 4.1 General

...

<br/>

# 5. Etc

## 5.1 General

...
