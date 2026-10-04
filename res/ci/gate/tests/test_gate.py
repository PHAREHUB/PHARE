"""Each diff check against small fixture repositories: what it must catch, and what it must leave alone.

run: python3 -m unittest discover -s res/ci/gate/tests
"""

import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import added_lines  # noqa: E402
import assert_scan  # noqa: E402
import gitdiff  # noqa: E402
import pharein_keywords  # noqa: E402
import pr_body  # noqa: E402
import tests_removed  # noqa: E402
import threads  # noqa: E402


def sh(cwd, *args):
    subprocess.run(args, cwd=cwd, check=True, capture_output=True)


class Repo:
    """A throwaway git repo: write files at 'base', commit, change them on 'head', commit."""

    def __init__(self, files):
        self.dir = tempfile.TemporaryDirectory()
        self.path = self.dir.name
        sh(self.path, "git", "init", "-q", "-b", "base")
        sh(self.path, "git", "config", "user.email", "t@t")
        sh(self.path, "git", "config", "user.name", "t")
        self.write(files)
        sh(self.path, "git", "add", "-A")
        sh(self.path, "git", "commit", "-qm", "base")
        sh(self.path, "git", "checkout", "-qb", "head")

    def write(self, files):
        for name, content in files.items():
            full = os.path.join(self.path, name)
            os.makedirs(os.path.dirname(full), exist_ok=True)
            with open(full, "w") as f:
                f.write(textwrap.dedent(content))

    def change(self, files=None, rename=None):
        if rename:
            sh(self.path, "git", "mv", *rename)
        self.write(files or {})
        sh(self.path, "git", "add", "-A")
        sh(self.path, "git", "commit", "-qm", "head")
        return self

    def run(self, module):
        try:
            return module.check("base", "head", cwd=self.path)
        finally:
            self.dir.cleanup()


PY_TEST = """\
    import unittest
    import numpy as np

    atol = 1e-8

    class T(unittest.TestCase):
        def test_a(self):
            np.random.seed(42)
            np.testing.assert_allclose(1.0, 1.0, atol=1e-6)
            self.assertAlmostEqual(1.0, 1.0, 7)
            self.assertTrue(True)

        def test_b(self):
            assert abs(1.0 - 1.0) < 1e-10
"""

CPP_TEST = """\
    #include <gtest/gtest.h>
    TEST(Suite, a)
    {
        EXPECT_NEAR(x, y, 1e-12);
        EXPECT_DOUBLE_EQ(u, v);
        EXPECT_EQ(n, 3);
    }
"""


def py(**replacements):
    text = PY_TEST
    for old, new in replacements.values():
        assert old in text, old
        text = text.replace(old, new)
    return text


def cpp(*pairs):
    text = CPP_TEST
    for old, new in pairs:
        assert old in text, old
        text = text.replace(old, new)
    return text


class AddedLines(unittest.TestCase):
    def flagged(self, files, change=None, rename=None):
        return Repo(files).change(change, rename).run(added_lines)

    def test_python_skip(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": py(a=("        def test_b", "        @unittest.skip('slow')\n        def test_b"))})
        self.assertEqual(len(f), 1)

    def test_self_skiptest(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": py(a=("np.random.seed(42)", "self.skipTest('later')"))})
        self.assertEqual(len(f), 1)

    def test_gtest_skip_and_disabled(self):
        f = self.flagged({"tests/test_x.cpp": CPP_TEST}, {"tests/test_x.cpp": cpp(("TEST(Suite, a)", "TEST(Suite, DISABLED_a)"), ("EXPECT_EQ(n, 3);", "GTEST_SKIP();"))})
        self.assertEqual(len(f), 2)

    def test_commented_registration(self):
        cm = "phare_python3_exec(9 x test_x.py ${DIR})\n"
        f = self.flagged({"tests/CMakeLists.txt": cm}, {"tests/CMakeLists.txt": "#" + cm})
        self.assertEqual(len(f), 1)

    def test_commented_python_test(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": py(a=("        def test_b(self):\n            assert", "        # def test_b(self):\n        #     assert"))})
        self.assertGreaterEqual(len(f), 1)

    def test_rename_off(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, rename=("tests/test_x.py", "tests/test_x.py.off"))
        self.assertEqual(len(f), 1)

    def test_new_test_is_fine(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": PY_TEST + "\n    def test_c(self):\n        self.assertTrue(True)\n"})
        self.assertEqual(f, [])

    def test_non_test_file_ignored(self):
        f = self.flagged({"src/x.py": "x = 1\n"}, {"src/x.py": "@unittest.skip\ndef f(): pass\n"})
        self.assertEqual(f, [])


    def test_core_includes_samrai(self):
        f = self.flagged({"src/core/x.hpp": "#include <vector>\n"}, {"src/core/x.hpp": '#include <vector>\n#include "SAMRAI/hier/Box.h"\n#include <amr/y.hpp>\n'})
        self.assertEqual(len(f), 2)
        self.assertTrue(all(x.label is None for x in f))  # layering cannot be overridden

    def test_amr_may_include_samrai(self):
        self.assertEqual(self.flagged({"src/amr/x.hpp": "\n"}, {"src/amr/x.hpp": "#include <SAMRAI/hier/Box.h>\n"}), [])

    def test_using_namespace_in_header(self):
        f = self.flagged({"src/amr/x.hpp": "\n"}, {"src/amr/x.hpp": "using namespace PHARE::core;\nvoid f() {\n    using namespace std;\n}\n"})
        self.assertEqual(len(f), 1)  # namespace scope only (NamespaceIndentation: Inner)
        self.assertEqual(f[0].label, "convention-ok")

    def test_new_cpp_under_src(self):
        f = self.flagged({"src/x.hpp": "\n"}, {"src/core/y.cpp": "int y;\n"})
        self.assertEqual(len(f), 1)
        self.assertEqual(f[0].label, "convention-ok")


class AssertScan(unittest.TestCase):
    def flagged(self, files, change):
        return Repo(files).change(change).run(assert_scan)

    def py_flags(self, **replacements):
        return self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": py(**replacements)})

    def test_keyword_tolerance_loosened(self):
        f = self.py_flags(a=("atol=1e-6", "atol=1e-4"))
        self.assertEqual(len(f), 1)
        self.assertIn("loosened", f[0].message)

    def test_positional_places_loosened(self):
        f = self.py_flags(a=("1.0, 7)", "1.0, 3)"))
        self.assertEqual(len(f), 1)

    def test_assert_bound_loosened(self):
        f = self.py_flags(a=("< 1e-10", "< 1e-3"))
        self.assertEqual(len(f), 1)

    def test_module_tolerance_loosened(self):
        f = self.py_flags(a=("atol = 1e-8", "atol = 1e-6"))
        self.assertEqual(len(f), 1)

    def test_assertion_removed(self):
        f = self.py_flags(a=("            self.assertTrue(True)\n", ""))
        self.assertEqual(len(f), 1)
        self.assertIn("removed", f[0].message)

    def test_seed_changed(self):
        f = self.py_flags(a=("seed(42)", "seed(7)"))
        self.assertEqual(len(f), 1)
        self.assertIn("seed", f[0].message)

    def test_tightened_is_fine(self):
        self.assertEqual(self.py_flags(a=("atol=1e-6", "atol=1e-9"), b=("1.0, 7)", "1.0, 9)")), [])

    def test_unrelated_edit_is_fine(self):
        self.assertEqual(self.py_flags(a=("import numpy as np", "import numpy as np\n    import os")), [])

    def cpp_flags(self, *pairs):
        return self.flagged({"tests/test_x.cpp": CPP_TEST}, {"tests/test_x.cpp": cpp(*pairs)})

    def test_cpp_near_loosened(self):
        f = self.cpp_flags(("1e-12", "1e-8"))
        self.assertEqual(len(f), 1)

    def test_cpp_exact_replaced(self):
        f = self.cpp_flags(("EXPECT_DOUBLE_EQ(u, v);", "EXPECT_NEAR(u, v, 1e-3);"))
        self.assertEqual(len(f), 1)

    def test_cpp_check_removed(self):
        f = self.cpp_flags(("        EXPECT_EQ(n, 3);\n", ""))
        self.assertEqual(len(f), 1)

    def test_cpp_new_check_is_fine(self):
        self.assertEqual(self.cpp_flags(("EXPECT_EQ(n, 3);", "EXPECT_EQ(n, 3);\n    EXPECT_EQ(m, 4);")), [])


class TestsRemoved(unittest.TestCase):
    def flagged(self, files, change=None, rename=None):
        return Repo(files).change(change, rename).run(tests_removed)

    def test_python_test_removed(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": py(a=("        def test_b(self):\n            assert abs(1.0 - 1.0) < 1e-10\n", ""))})
        self.assertEqual([x.message for x in f], ["test T.test_b removed or renamed"])
        self.assertEqual(f[0].label, "tests-removed-ok")

    def test_python_test_hidden_by_underscore(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": py(a=("        def test_b(self):", "        def _test_b(self):"))})
        self.assertEqual(len(f), 1)

    def test_gtest_removed_and_commented(self):
        f = self.flagged({"tests/test_x.cpp": CPP_TEST}, {"tests/test_x.cpp": cpp(("TEST(Suite, a)", "// TEST(Suite, a)"))})
        self.assertEqual(len(f), 1)
        self.assertIn("Suite.a", f[0].message)

    def test_disabled_is_left_to_added_lines(self):
        f = self.flagged({"tests/test_x.cpp": CPP_TEST}, {"tests/test_x.cpp": cpp(("TEST(Suite, a)", "TEST(Suite, DISABLED_a)"))})
        self.assertEqual(f, [])

    def test_deleted_test_file(self):
        repo = Repo({"tests/test_x.cpp": CPP_TEST, "tests/keep.txt": "x\n"})
        os.remove(os.path.join(repo.path, "tests/test_x.cpp"))
        f = repo.change().run(tests_removed)
        self.assertEqual(len(f), 1)

    def test_test_moved_to_other_file(self):
        repo = Repo({"tests/test_x.cpp": CPP_TEST, "tests/test_y.cpp": "\n"})
        os.remove(os.path.join(repo.path, "tests/test_x.cpp"))
        self.assertEqual(repo.change({"tests/test_y.cpp": CPP_TEST}).run(tests_removed), [])

    def test_registration_removed_or_level_changed(self):
        cm = "phare_python3_exec(9 x test_x.py ${DIR})\nadd_no_mpi_phare_test(${PROJECT_NAME}\n    ${CMAKE_CURRENT_BINARY_DIR})\n"
        f = self.flagged({"tests/CMakeLists.txt": cm}, {"tests/CMakeLists.txt": cm.replace("(9 ", "(11 ").replace("add_no_mpi", "# add_no_mpi")})
        self.assertEqual(len(f), 2)

    def test_registration_reformatted_is_fine(self):
        cm = "add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})\n"
        self.assertEqual(self.flagged({"tests/CMakeLists.txt": cm}, {"tests/CMakeLists.txt": cm.replace(" $", "\n  $")}), [])

    def test_new_unregistered_test_file_warns(self):
        f = self.flagged({"tests/CMakeLists.txt": "\n"}, {"tests/test_new.py": PY_TEST})
        self.assertEqual([x.level for x in f], ["warning"])

    def test_new_registered_test_file_is_fine(self):
        f = self.flagged({"tests/CMakeLists.txt": "\n"}, {"tests/test_new.py": PY_TEST, "tests/CMakeLists.txt": "add_python3_test(t test_new.py ${D})\n"})
        self.assertEqual(f, [])

    def test_added_test_is_fine(self):
        f = self.flagged({"tests/test_x.py": PY_TEST}, {"tests/test_x.py": PY_TEST + "\n    def test_c(self):\n        self.assertTrue(True)\n"})
        self.assertEqual(f, [])


SIM = """\
    def checker(func):
        def wrapper(sim, **kwargs):
            accepted_keywords = [
                "cells",
                "dl",
            ]
            accepted_keywords += ["strict"]
"""


class PhareinKeywords(unittest.TestCase):
    def flagged(self, new):
        path = "pyphare/pyphare/pharein/simulation.py"
        return Repo({path: SIM}).change({path: new}).run(pharein_keywords)

    def test_keyword_removed(self):
        f = self.flagged(SIM.replace('            "dl",\n', ""))
        self.assertEqual(len(f), 1)
        self.assertEqual((f[0].level, f[0].label), ("error", "api-break-ok"))
        self.assertIn("'dl'", f[0].message)

    def test_augmented_keyword_removed(self):
        self.assertEqual(len(self.flagged(SIM.replace('["strict"]', "[]"))), 1)

    def test_mandatory_made_optional_is_fine(self):
        sim = SIM.replace('["strict"]', '["strict"]\n            mandatory_keywords = ["path"]')
        path = "pyphare/pyphare/pharein/simulation.py"
        f = Repo({path: sim}).change({path: sim.replace('["strict"]', '["strict", "path"]').replace('mandatory_keywords = ["path"]', "mandatory_keywords = []")}).run(pharein_keywords)
        self.assertEqual(f, [])

    def test_keyword_added_is_notice(self):
        f = self.flagged(SIM.replace('"dl",', '"dl",\n            "dt",'))
        self.assertEqual([x.level for x in f], ["notice"])


TEMPLATE_BODY = """## Issue

<!-- Fixes #... -->
{issue}
## What this implements

<!-- ... -->

## Tests

<!-- What tests?
-->
{tests}
"""


class PrBody(unittest.TestCase):
    def test_filled(self):
        self.assertEqual(pr_body.check(TEMPLATE_BODY.format(issue="Fixes #12", tests="test_x.py")), [])

    def test_written_justification_is_fine(self):
        self.assertEqual(pr_body.check(TEMPLATE_BODY.format(issue="Typo, no issue needed.", tests="none: docs only")), [])

    def test_empty_sections(self):
        f = pr_body.check(TEMPLATE_BODY.format(issue="", tests=""))
        self.assertEqual([x.level for x in f], ["error", "warning"])

    def test_missing_tests_is_warning(self):
        self.assertEqual([x.level for x in pr_body.check("## Issue\n#3\n")], ["warning"])

    def test_closes_anywhere_counts_as_issue(self):
        for body in ("closes #1139", "Fixes #1196.\n\nlong text", "relates to #1110 but maybe doesn't close", "resolves PHAREHUB/PHARE#12"):
            self.assertEqual([x.level for x in pr_body.check(body)], ["warning"], body)

    def test_plain_text_fails(self):
        self.assertIn("error", [x.level for x in pr_body.check("enables dry run of harris")])

    def test_closes_in_template_comment_does_not_count(self):
        self.assertIn("error", [x.level for x in pr_body.check(TEMPLATE_BODY.format(issue="", tests="t"))])

    def test_crlf(self):
        self.assertEqual(pr_body.check(TEMPLATE_BODY.format(issue="#1", tests="t").replace("\n", "\r\n")), [])


def thread(resolved, starter, typename="User"):
    return {
        "isResolved": resolved,
        "path": "src/x.hpp",
        "line": 3,
        "comments": {"nodes": [{"url": "u", "author": {"login": starter, "__typename": typename}}]},
    }


class Threads(unittest.TestCase):
    def test_unresolved_human_fails(self):
        self.assertEqual([f.level for f in threads.evaluate([thread(False, "rev")])], ["error"])

    def test_bot_threads_ignored(self):
        bots = [thread(False, "coderabbitai", typename="Bot"), thread(False, "copilot[bot]")]
        self.assertEqual(threads.evaluate(bots), [])

    def test_resolved_is_fine(self):
        self.assertEqual(threads.evaluate([thread(True, "rev")]), [])


class Report(unittest.TestCase):
    def run_report(self, findings, labels):
        os.environ["GATE_LABELS"] = json.dumps(labels)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                return gitdiff.report(findings, "t")
        finally:
            del os.environ["GATE_LABELS"]

    def test_label_overrides_its_errors_only(self):
        f = [gitdiff.Finding("a", 1, "c", "m", label="ok")]
        self.assertEqual(self.run_report(f, ["ok"]), 0)
        self.assertEqual(f[0].level, "warning")
        self.assertEqual(self.run_report([gitdiff.Finding("a", 1, "c", "m", label="ok")], ["other"]), 1)
        self.assertEqual(self.run_report([gitdiff.Finding("a", 1, "c", "m")], ["ok"]), 1)

    def test_warnings_do_not_fail(self):
        self.assertEqual(self.run_report([gitdiff.Finding("a", 1, "c", "m", level="warning")], []), 0)

    def test_annotation_escaping(self):
        a = gitdiff.Finding("a,b:c", 2, "chk", "50%\nnext").annotation()
        self.assertEqual(a, "::error file=a%2Cb%3Ac,line=2,title=chk::50%25%0Anext")


if __name__ == "__main__":
    unittest.main()
