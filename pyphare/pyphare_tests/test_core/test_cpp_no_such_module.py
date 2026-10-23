import sys
import types
import unittest
from unittest.mock import patch

from pyphare import cpp


class CppNoSuchModuleTest(unittest.TestCase):
    def setUp(self):
        self.sim = types.SimpleNamespace(
            ndim=1, interp_order=1, refined_particle_nbr=3, mhd_timestepper=None
        )

    def test_message_names_location_and_fixes(self):
        with patch("pyphare.cpp.built_simulator_ids", return_value=["1_1_2"]), patch(
            "pyphare.cpp._pybindlibs_location", return_value="/some/where/pybindlibs"
        ):
            message = cpp._no_such_module_message(self.sim, "pybindlibs.cpp_1_1_3")
        for expected in [
            "1,1,3",
            "1_1_2",
            "/some/where/pybindlibs",
            "res/sim",
            "pip install . -C cmake.define.PHARE_PERMUTATIONS=",
            "every permutation",
            "a copy of res/sim/all.txt plus this line",
            "installs exactly",
        ]:
            self.assertIn(expected, message)

    def test_location_when_not_importable(self):
        with patch.dict(sys.modules, {"pybindlibs": None}):
            self.assertEqual(cpp._pybindlibs_location(), "pybindlibs not importable")


if __name__ == "__main__":
    unittest.main()
