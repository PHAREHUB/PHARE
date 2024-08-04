
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

Write a :doc:`simulation_inputs` and run the following command:


.. code-block:: bash

      python3 /path/to/my_script.py


Executable PHARE initialization script
--------------------------------------

Running PHARE from python basically consists in declaring desired blocks
detailed in a :doc:`simulation_inputs` file, and adding a `main` function to run
the simulation.

If running the simulation from python, you need to define a `main` function that
basically runs the simulation. Here is a small example:


.. code-block:: python

    from pyphare.simulator.simulator import Simulator
    from pyphare.pharein import Simulation
    from pyphare.pharein import MaxwellianFluidModel


    # define your simulation parameters
    # and initial condition...

    sim = Simulation(
    # ...
    )

    MaxwellianFluidModel(
    # ...
    )


    def main():
        simulator = Simulator(sim)
        simulator.run()



The Simulator
-------------

The simulator is used to initialize and run the simulation.
It is given the `Simulation` object as an argument, as in the 
simple example above. The simplest way to start the simulation is
then to call the `run()` method.


.. autoclass:: pyphare.simulator.simulator.Simulator
   :members:

