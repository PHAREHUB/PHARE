# PHARE Conventions

### Reference document for the code base


# Sections

1. C++
2. Python
3. CMake
4. Tests
5. Etc

<br/>

# 1. C++

## 1.1 General

...


<br/>

# 2. Python

## 2.1 General

...

## 2.2 dependencies and imports

Third party depenencies are stated in the file `requirements.txt` in the project root.
Fewer dependencies is generally better but there should be a cost/benefit assessment for adding new dependencies.

### 2.2.1 Python file import structure.

Generally, we want to avoid importing any dependency at the top of a python script that may rely on binary libraries.

Exceptions to this are things like numpy, which are widely used and tested.

Things to expressly avoid importing at the top of a python script are

- h5py
- mpi4py
- scipy.optimize

The first two are noted as they can, and will pull in system libraries such as libmpi.so and libhdf5.so, which may not be the libraries which were used during PHARE build time, this can cause issues at runtime.

scipy.optimize relies on system libraries which may not be available at runtime.

The gist here is to only import these libraries at function scope when you actually need them, so that python files can be imported
or scanned for tests and not cause issues during these operations, until the functions are used at least.

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

