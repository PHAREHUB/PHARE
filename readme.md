
[![License: GPL v3](https://img.shields.io/badge/PHARE-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CPP17](https://img.shields.io/badge/Language-C++17-blue.svg)]()
[![LGTM Total alerts](https://img.shields.io/lgtm/alerts/g/PHAREHUB/PHARE.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/PHAREHUB/PHARE/alerts/)
[![codefactor](https://www.codefactor.io/repository/github/PHAREHUB/PHARE/badge?style=plastic)](https://www.codefactor.io/repository/github/PHAREHUB/PHARE/badge?style=plastic)

![](https://user-images.githubusercontent.com/3200931/95620089-f13ebb80-0a6f-11eb-9af3-a1db08004bcc.png)


# Parallel Hybrid PIC code with Adaptive mesh REfinement



PHARE is a Hybrid Particle-In-Cell (PIC) code. It solves the evolution of the Vlasov
equation of an arbitrary number of ion populations in a Lagrangian way. Electrons are
modeled as a single fluid. Their momentum equation is used to compute the electric field, 
assuming quasineutrality.

Using Adaptive Mesh Refinement, provided by the library [SAMRAI](https://github.com/llnl/samrai),
PHARE aims at filling the gap between sub-ion scales and large "MHD" scales by increasing the mesh
resolution wherever the solution needs it.

PHARE is still under development.


## Software Licence

PHARE is an open-source projet licenced under the GPLv3. Please refer to [LICENCE.TXT](LICENCE.TXT)


## Install

For system library requirements see the following [Docker File](https://github.com/PHARCHIVE/phare-teamcity-agent/blob/master/Dockerfile)


For Python API, install dependencies:

```
  python3 -m pip install -r requirements.txt
```

Build with CMake and latest [SAMRAI](https://github.com/llnl/samrai):

```
  mkdir build; cd build; cmake ..; make
```

Build with CMake and specific [SAMRAI](https://github.com/llnl/samrai):

```
  mkdir build; cd build; cmake .. -DSAMRAI_ROOT=/path/to/SAMRAI/install; make
```



## Developers


### Contributing

All contributions are welcome, although the code is still at an early phase of
its development. If you are interested in participating to the project for an internship,
PhD, PostDoc, [contact us](mailto:phare@lpp.polytechnique.fr).


### Active core team

- [Nicolas Aunai](https://github.com/nicolasaunai)
- [Philip Deegan](https://github.com/PhilipDeegan)
- [Roch Smets](https://github.com/rochsmets)
- [Andrea Ciardi](https://sites.google.com/site/andreaciardihomepage/home)


### Former collaborators

- [Thibault Payet](https://github.com/monwarez)



## Publications

- Smets R., Aunai, N., Ciardi, A., Drouin M., Campos-Pinto M., Deegan P., A new method to dispatch split particles in Particle-In-Cell codes, Computer Physics Communications, https://doi.org/10.1016/j.cpc.2020.107666

## Acknowledgments


The project PHARE is supported by:

- [Plas@Par](http://www.plasapar.com)
- [Laboratory of Plasma Physics (LPP)](https://www.lpp.polytechnique.fr)
- [LERMA](https://lerma.obspm.fr)
