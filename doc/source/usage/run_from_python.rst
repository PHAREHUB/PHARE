
=====================
Run PHARE from python
=====================

PHARE can run as a python script.

Python dependencies
-------------------

PHARE requires a minimmum version of python 3.8 to run properly.
Make sure `python3` shows the version is at least 3.8. 
Python package dependencies are listed in `requirements.txt` file.
Install dependencies for the user:

.. code-block:: bash

     pip3 install --user -r requirements.txt

Install dependencies in a virtual environment:

.. code-block:: bash

     python3 -m venv phare_venv
     source phare_venv/bin/activate
     pip3 install -r requirements.txt


Running PHARE
-------------
First, make sure it is accessible to python.
Assuming PHARE source directory is in  `/path/to/PHARE`, and the build directory is
 `/path/to/build/`, then use the following to let python know where to find PHARE:


.. code-block:: bash

     export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH


Write a [simulation input script](../simulation_inputs.md) and run the following command:


.. code-block:: bash

      python3 /path/to/my_script.py

