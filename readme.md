

<div align="center">
<img src="https://user-images.githubusercontent.com/3200931/213649230-da966bb5-9384-4c58-8e77-24ebe46cf2bf.png" alt="PHARE snapshot">

<p align="center">
<a  href="https://www.gnu.org/licenses/gpl-3.0">
<img src="https://img.shields.io/badge/PHARE-GPL%20v3-blue.svg" alt="License: GPL v3" />
</a>
<a  href="">
<img src="https://img.shields.io/badge/Language-C++17-blue.svg" alt="CPP17" />
</a>
<a  href="https://phare.readthedocs.io/en/latest/?badge=latest">
<img src="https://readthedocs.org/projects/phare/badge/?version=latest" alt="Documentation Status" />
</a>
</p>
</div>

# Parallel Hybrid PIC code with Adaptive mesh REfinement


PHARE is a Hybrid Particle-In-Cell (PIC) code. It solves the evolution of the Vlasov equation
of an arbitrary number of ion populations in a Lagrangian way. Electrons are modeled as a single fluid.
Their momentum equation is used to compute the electric field, assuming quasineutrality.

Using Adaptive Mesh Refinement, provided by the library [SAMRAI](https://github.com/llnl/samrai),
PHARE aims at filling the gap between sub-ion scales and large "MHD" scales by increasing the mesh 
resolution wherever the solution needs it.

**WARNING** - PHARE is under development ;-)




## Software Licence

PHARE is an open-source projet licenced under the GPLv3. Please refer to [LICENCE.TXT](LICENCE.TXT)


## Building the code

Basic tools and library requirements:

```
- git
- cmake
- make/ninja
- C++ and Fortran compiler
- MPI
- Parallel HDF5
- Python 3.x devel package
```



PHARE input and post-processing scripts are in python. Install dependencies with:

```
  python3 -m pip install -r requirements.txt
```

PHARE depends on [SAMRAI](https://github.com/llnl/samrai) to manager the adaptive mesh refinement. You can either

- build PHARE with the latest version of [SAMRAI](https://github.com/llnl/samrai):

```
  mkdir build; cd build; cmake ..; make
```

- build PHARE with a pre-installed version of [SAMRAI](https://github.com/llnl/samrai):

```
  mkdir build; cd build; cmake .. -DSAMRAI_ROOT=/path/to/SAMRAI/install; make
```

## Documentation

- [https://phare.readthedocs.io/en/latest/](https://phare.readthedocs.io/en/latest/)


## Developers


### Contributing

All contributions are welcome. If you are interested in participating to the project for an internship, PhD, PostDoc, [contact us](mailto:phare@lpp.polytechnique.fr).

### Development environnment

For system library requirements see the following [Docker File](https://github.com/PHARCHIVE/phare-teamcity-agent/blob/master/Dockerfile)

### Active core team

- [Nicolas Aunai](https://github.com/nicolasaunai)
- [Philip Deegan](https://github.com/PhilipDeegan)
- [Roch Smets](https://github.com/rochsmets)
- [Andrea Ciardi](https://sites.google.com/site/andreaciardihomepage/home)


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




## Acknowledgments


We acknowledge support from:

- [Plas@Par](http://www.plasapar.sorbonne-universite.fr/)
- [Laboratory of Plasma Physics (LPP)](https://www.lpp.polytechnique.fr)
- [LERMA](https://lerma.obspm.fr)
- [CNES](www.cnes.fr)

