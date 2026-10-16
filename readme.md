

<div align="center">
<img src="https://github.com/user-attachments/assets/5e330ec4-99c7-4872-a5d4-73c618aaba07" alt="PHARE snapshot">


<p align="center">
<a  href="https://www.gnu.org/licenses/gpl-3.0">
<img src="https://img.shields.io/badge/PHARE-GPL%20v3-blue.svg" alt="License: GPL v3" />
</a>
<a  href="">
<img src="https://img.shields.io/badge/Language-C++20-blue.svg" alt="CPP20" />
</a>
<a  href="https://phare.readthedocs.io/en/latest/?badge=latest">
<img src="https://readthedocs.org/projects/phare/badge/?version=latest" alt="Documentation Status" />
</a>
<a href="https://doi.org/10.1016/j.cpc.2023.108966">
<img src="https://img.shields.io/badge/DOI%20(paper)-10.1016/j.cpc.2023.108966-blue.svg"
</a>
<a href="https://arxiv.org/abs/2210.14580">
<img src="https://img.shields.io/badge/DOI%20(arXiv)-https://arxiv.org/abs/2210.14580-blue.svg"
</a>
</p>
</div>

# PHARE: Parallel Hybrid PIC code with Adaptive mesh REfinement


PHARE is a Hybrid Particle-In-Cell (PIC) code. It solves the evolution of the Vlasov equation
of an arbitrary number of ion populations in a Lagrangian way. Electrons are modeled as a single fluid.
Their momentum equation is used to compute the electric field, assuming quasineutrality.

Using Adaptive Mesh Refinement, provided by the library [SAMRAI](https://github.com/llnl/samrai),
PHARE aims at filling the gap between sub-ion scales and large "MHD" scales by increasing the mesh 
resolution wherever the solution needs it.

**WARNING** - PHARE is under development ;-)




## Software Licence

PHARE is an open-source projet licenced under the GPLv3. Please refer to [LICENCE.TXT](LICENCE.TXT)


## Install PHARE

Basic tools and library requirements:

```
- git
- cmake
- make/ninja
- C++ and Fortran compiler
- MPI
- Parallel HDF5
- Python 3.11+ devel package
```

From a source checkout, `pip install .` builds PHARE with CMake and installs the `pyphare` python package together with the compiled simulator modules. Nothing needs to be added to `PYTHONPATH` or `LD_LIBRARY_PATH` afterwards:

```
  python3 -m venv phare_venv
  source phare_venv/bin/activate
  pip install .
```

PHARE depends on [SAMRAI](https://github.com/llnl/samrai) to manage the adaptive mesh refinement. By default the latest version is built from `subprojects/`. To use a pre-installed SAMRAI instead:

```
  pip install . -C cmake.define.SAMRAI_ROOT=/path/to/SAMRAI/install
```

By default every simulation permutation listed in `res/sim/all.txt` is compiled, which takes a long time. To build only the permutations listed in your own file:

```
  pip install . -C cmake.define.PHARE_PERMUTATIONS=/path/to/permutations.txt
```

The number of simulator modules compiled concurrently is limited according to the available memory. Override it with `-C cmake.define.PHARE_MODULE_COMPILE_JOBS=<n>` (`0` disables the limit) or `-C cmake.define.PHARE_MODULE_COMPILE_MEMORY_MB=<mb>` (default 2048, memory assumed per compile job).

The build directory `_skbuild/` is kept between installs, so rebuilds are incremental. `SAMRAI_ROOT` and `PHARE_PERMUTATIONS` are not remembered from one install to the next: pass them again each time. Delete `_skbuild/` after changing compilers, flags or loaded modules.

With [uv](https://docs.astral.sh/uv/), `uv sync` creates a `.venv`, installs the python dependencies and builds PHARE the same way. Since `-C` options are not available, pass CMake definitions through the `SKBUILD_CMAKE_DEFINE` environment variable (semicolon-separated `KEY=VALUE`):

```
  SKBUILD_CMAKE_DEFINE="SAMRAI_ROOT=/path/to/SAMRAI/install;PHARE_PERMUTATIONS=/path/to/permutations.txt" uv sync
```

uv does not detect changes to these options, nor to a permutation file outside the repository. After changing them, run `uv sync --reinstall-package pyphare` (pip reinstalls on every `pip install .`).

Do not run `ninja` or `cmake --build` by hand inside `_skbuild/`: after a CMake file changes, the automatic reconfigure does not have your `-C` options and falls back to the defaults. Reinstall with pip or uv instead.

To install on a cluster, or to choose compilers and flags, see [Run PHARE from python](https://phare.readthedocs.io/en/latest/usage/run_from_python.html).

## Documentation

- [https://phare.readthedocs.io/en/latest/](https://phare.readthedocs.io/en/latest/)


## Developers


### Contributing

All contributions are welcome. If you are interested in participating to the project for an internship, PhD, PostDoc, [contact us](mailto:phare@lpp.polytechnique.fr).

### Development environnment

For system library requirements see the following [Docker File](https://github.com/PHARCHIVE/phare-teamcity-agent/blob/master/Dockerfile)

Developers who keep their own CMake build:

```
  mkdir build; cd build; cmake .. -DSAMRAI_ROOT=/path/to/SAMRAI/install; make
```

should install only the python dependencies, since a plain `uv sync` builds PHARE with every permutation:

```
  uv sync --no-install-project
  export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/PHARE/build:$PYTHONPATH
```

Do not mix this build with an environment where `uv sync` or `uv run` installed PHARE: that install is editable, and its import finder takes precedence over `PYTHONPATH`, so the installed `pybindlibs` is used instead of your build (a non-editable `pip install .` does not do this). `uv sync` and `uv run` also manage `.venv` as the project environment, and may recreate it (for example with a uv-managed python), dropping packages installed by hand. Use a separate environment for your own build, or keep uv from touching it: `uv sync --no-install-project`, `uv run --no-sync <cmd>`, or `UV_NO_SYNC=1` for scripts calling `uv run`. A `.venv` shared between checkouts (for example symlinked into git worktrees) is synced by whichever checkout runs uv.

The clang-format and ruff versions checked by CI are pinned in the `lint` dependency group of `pyproject.toml`. `uv sync` installs them, or with pip >= 25.1:

```
  python3 -m pip install --group lint
```

### Active core team

- [Nicolas Aunai](https://github.com/nicolasaunai)
- [Philip Deegan](https://github.com/PhilipDeegan)
- [Roch Smets](https://github.com/rochsmets)
- [Andrea Ciardi](https://sites.google.com/site/andreaciardihomepage/home)
- [Ulysse Caromel](https://github.com/UCaromel)
- [Ivan Girault](https://github.com/igirault)


### Former collaborators

- [Thibault Payet](https://github.com/monwarez)



## Publications

To cite PHARE : 

```bib
@article{AUNAI2024108966,
title = {PHARE: Parallel hybrid particle-in-cell code with patch-based adaptive mesh refinement},
journal = {Computer Physics Communications},
volume = {295},
pages = {108966},
year = {2024},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2023.108966},
url = {https://www.sciencedirect.com/science/article/pii/S0010465523003119},
author = {Nicolas Aunai and Roch Smets and Andrea Ciardi and Philip Deegan and Alexis Jeandet and Thibault Payet and Nathan Guyot and Loic Darrieumerlou},
keywords = {Particle in cell, Adaptive mesh refinement, Collisionless plasmas},
}
```
- [ArXiv free version of PHARE's CPC paper](https://arxiv.org/abs/2210.14580)
- [Particle Splitting](https://doi.org/10.1016/j.cpc.2020.107666): Roch Smets et al. *A new method to dispatch split particles in Particle-In-Cell codes.* Computer Physics Communications, 2021, 261, pp.107666. ⟨10.1016/j.cpc.2020.107666⟩. ([Free Version here](https://universite-paris-saclay.hal.science/hal-03203893/)) 


## More

- follow us on [Youtube](https://www.youtube.com/@PHARE-v1u)


## Acknowledgments


We acknowledge support from:

- [Plas@Par](http://www.plasapar.sorbonne-universite.fr/)
- [Laboratory of Plasma Physics (LPP)](https://www.lpp.polytechnique.fr)
- [LERMA](https://lerma.obspm.fr)
- [CNES](www.cnes.fr)
- [CIEDS](https://www.ip-paris.fr/cieds) and [AID](https://www.defense.gouv.fr/aid)

