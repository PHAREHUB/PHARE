import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from pyphare import cpp


def make_pybindlibs(parent):
    package = Path(parent) / "pybindlibs"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    return package


def fake_module(package):
    module = types.ModuleType("pybindlibs")
    module.__file__ = str(package / "__init__.py")
    module.__path__ = [str(package)]
    return module


class PybindlibsShadowingTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.developer = make_pybindlibs(self.tmp / "build")
        self.installed = make_pybindlibs(self.tmp / "site-packages")
        self.sim = types.SimpleNamespace(
            ndim=1, interp_order=1, refined_particle_nbr=2, mhd_timestepper=None
        )
        checked = patch.object(cpp, "_pybindlibs_checked", False)
        checked.start()
        self.addCleanup(checked.stop)

    def environment(self, loaded, *pythonpath):
        modules = patch.dict(sys.modules, {"pybindlibs": fake_module(loaded)})
        env = patch.dict(
            os.environ, {"PYTHONPATH": os.pathsep.join(str(p) for p in pythonpath)}
        )
        modules.start()
        env.start()
        self.addCleanup(modules.stop)
        self.addCleanup(env.stop)

    def test_shadowed_developer_build_raises(self):
        self.environment(self.installed, self.tmp / "empty", self.tmp / "build")
        with self.assertRaises(ImportError) as raised:
            cpp.built_simulator_ids()
        self.assertIs(type(raised.exception), ImportError)
        message = str(raised.exception)
        for expected in [
            str(self.installed.resolve()),
            str(self.developer.resolve()),
            "PYTHONPATH",
            "uv sync --no-install-project",
            "pip uninstall pyphare",
        ]:
            self.assertIn(expected, message)

    def test_every_entry_point_checks(self):
        self.environment(self.installed, self.tmp / "build")
        for entry in [
            lambda: cpp.cpp_lib(self.sim),
            cpp.cpp_etc_lib,
            cpp.built_simulator_ids,
            cpp._pybindlibs_location,
        ]:
            with self.subTest(entry=entry), self.assertRaises(ImportError) as raised:
                entry()
            self.assertIs(type(raised.exception), ImportError)

    def test_developer_build_in_use(self):
        self.environment(self.developer, self.tmp / "build")
        self.assertEqual(cpp.built_simulator_ids(), [])

    def test_symlinked_pythonpath_entry(self):
        (self.tmp / "link").symlink_to(self.tmp / "build")
        self.environment(self.developer, self.tmp / "link")
        self.assertEqual(cpp.built_simulator_ids(), [])

    def test_no_pybindlibs_on_pythonpath(self):
        self.environment(self.installed, self.tmp / "empty")
        self.assertEqual(cpp.built_simulator_ids(), [])

    def test_checked_once(self):
        self.environment(self.developer, self.tmp / "build")
        with patch.object(
            cpp,
            "_check_pybindlibs_not_shadowed",
            wraps=cpp._check_pybindlibs_not_shadowed,
        ) as check:
            cpp.built_simulator_ids()
            cpp._pybindlibs_location()
            cpp.built_simulator_ids()
        self.assertEqual(check.call_count, 1)


if __name__ == "__main__":
    unittest.main()
