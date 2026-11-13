"""Where the amr tests find the sample runs written by ``gen_sample_data.py``.

``PHARESEE_SAMPLE_DATA`` (set by ctest to ``<build dir>/sample_data``), else
``<repo>/build/sample_data``. The 3-rank run lives next to it, in ``<root>_mpi3``.
Layout: ``<root>/{1d,2d_static,2d_tagged,3d}/{h5,vtk}`` and ``<root>_mpi3/2d_tagged/{h5,vtk}``.
"""

import os
import unittest

ENV = "PHARESEE_SAMPLE_DATA"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), *[os.pardir] * 4))


def sample_root():
    """Directory holding the single-process sample runs."""
    return os.path.abspath(os.environ.get(ENV) or os.path.join(REPO, "build", "sample_data"))


def mpi3_root(root=None):
    """Directory holding the 3-rank sample run (sibling of ``root``)."""
    return (root or sample_root()).rstrip(os.sep) + "_mpi3"


SAMPLES = sample_root()
SAMPLES_MPI3 = mpi3_root(SAMPLES)


def run_dirs(case, root=SAMPLES):
    """(h5 dir, vtk dir) of one sample run, e.g. run_dirs("2d_tagged")."""
    return os.path.join(root, case, "h5"), os.path.join(root, case, "vtk")


def have(root=SAMPLES, cases=("2d_tagged",)):
    """True when every case of ``root`` has both its h5 and vtk folders."""
    return all(os.path.isdir(d) for c in cases for d in run_dirs(c, root))


def needs_samples(root=SAMPLES, cases=("2d_tagged",)):
    """Class/test decorator: skip unless ``have(root, cases)``."""
    return unittest.skipUnless(have(root, cases), f"needs {root} ({', '.join(cases)})")
