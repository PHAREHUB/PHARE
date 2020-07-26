# PHARE


[![License: GPL v3](https://img.shields.io/badge/PHARE-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CPP17](https://img.shields.io/badge/Language-C++17-blue.svg)]()
[![LGTM Total alerts](https://img.shields.io/lgtm/alerts/g/PHAREHUB/PHARE.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/PHAREHUB/PHARE/alerts/)
[![codefactor](https://www.codefactor.io/repository/github/PHAREHUB/PHARE/badge?style=plastic)](https://www.codefactor.io/repository/github/PHAREHUB/PHARE/badge?style=plastic)

## Introduction


PHARE stands for Parallel Hybrid particle in cell code with Adaptive Mesh REfinement.
Hybrid particle-in-cell (PIC) codes solve the Vlasov equation for ions in a Lagrangian way
and assume an electron fluid.

PHARE is in development, and aims at filling the gap between sub-ion scales and large "MHD" scales
by using Adaptive Mesh Refinement, enabling to dynamically increasing the resolution in critical regions while
keeping a lower resolution elsewhere.


## Software Licence

PHARE is an open-source projet licenced under the GPLv3. Please refer to [LICENCE.TXT](LICENCE.TXT)


## Install

For system library requirements see the file:
  https://github.com/PHARCHIVE/phare-teamcity-fc31/blob/master/Dockerfile

For Python API, install dependencies:
  python3 -m pip install -r requirements.txt

Build with CMake and latest SAMRAI
  mkdir build; cd build; cmake ..; make

Build with CMake and specific SAMRAI
  mkdir build; cd build; cmake .. -DSAMRAI_ROOT=/path/to/SAMRAI/install; make


## User Documentation

User documentation can be found here...


## Developers


### Contributing

- See [how to contribute]()
- read our [developer documentation]()

### Active core team

- Nicolas Aunai ([LPP](https://www.lpp.polytechnique.fr))
- Roch Smets ([LPP](https://www.lpp.polytechnique.fr))
- Andrea Ciardi ([LERMA](https://lerma.obspm.fr))
- Philip Deegan (https://github.com/Dekken)

### Former collaborators

- Thibault Payet ([LPP](https://github.com/monwarez))
- Mathieu Drouin



Acknowledgments
---------------

The project PHARE is supported by:

- [Plas@Par](http://www.plasapar.com)
- [Laboratory of Plasma Physics (LPP)](https://www.lpp.polytechnique.fr)
- [LERMA](https://lerma.obspm.fr)
