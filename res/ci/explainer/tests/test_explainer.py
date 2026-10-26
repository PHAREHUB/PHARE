"""Sorting, excerpts, report parsing and comment cleaning, on snippets of real TeamCity logs.

run: python3 -m unittest discover -s res/ci/explainer/tests
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import agent
import ci
import comment
import excerpt
import triage


def build(status_text, status="FAILURE"):
    return {
        "id": 1,
        "config": "Phare_Phare_BuildGithubPrClang",
        "status": status,
        "status_text": status_text,
        "pr": 1,
        "sha": "0" * 40,
    }


def log(*lines):
    return triage.clean("\n".join(lines))


# a harmless warning: TeamCity then creates a container per command and the build goes on
CONTAINER_WARNING = "[20:08:37]W:\t [Step 1/3] Failed to create a reusable container from the image registry/phare:43. A new container will be created per each command"
PULL_FAILED = "[20:09:34]W:\t\t [pull 129.104.6.172:32219/phare/teamcity-fedora_dep:43] Process exited with code 1"
CLANG_ERROR = (
    "[14:30:14] :\t [Step 3/5] FAILED: [code=1] src/python3/CMakeFiles/cpp_1_1_2.dir/cpp_simulator.cpp.o ",
    "[14:30:14] :\t [Step 3/5] In file included from /opt/buildagent/work/3bd323ee5a73a841/src/python3/cpp_simulator.cpp:9:",
    "[14:30:14] :\t [Step 3/5] /opt/buildagent/work/3bd323ee5a73a841/src/simulator/simulator.hpp:395:20: error: variable 'already_setup' set but not used [-Werror,-Wunused-but-set-variable]",
    "[14:30:14] :\t [Step 3/5]   395 |     if (bool const already_setup = multiphysInteg_ != nullptr)",
    "[14:30:14] :\t [Step 3/5] [96/277] ccache /usr/sbin/clang++ -DH5_HAVE_PARALLEL -I/opt/buildagent/work/3bd323ee5a73a841/src",
    "[14:30:14] :\t [Step 3/5] FAILED: [code=1] src/python3/CMakeFiles/cpp_1_1_3.dir/cpp_simulator.cpp.o ",
)
CTEST_FAILURE = (
    "[18:13:53] :\t [Step 5/5] 27/49 Test #29: py3_init-densities_mpi_n_2 ...............***Failed   32.41 sec",
    "[18:13:53] :\t [Step 5/5] validating dim=1",
    "[18:13:53] :\t [Step 5/5] AssertionError: densities differ",
    "[18:13:54] :\t [Step 5/5] 28/49 Test #30: py3_init-2d-fields_mpi_n_2 ...............   Passed   12.00 sec",
    "[18:20:00] :\t [Step 5/5] 98% tests passed, 1 tests failed out of 49",
    "[18:20:00] :\t [Step 5/5] The following tests FAILED:",
    "[18:20:00] :\t [Step 5/5] \t 29 - py3_init-densities_mpi_n_2 (Failed)",
    "[18:20:00] :\t [Step 5/5] Process exited with code 8",
)


class Triage(unittest.TestCase):
    def test_python_test_failure_hidden_by_teamcity_count_is_the_prs(self):
        verdict = triage.classify(
            build(
                "Tests passed: 1621; exit code 8 (Step: run mpi test (Command Line))"
            ),
            log(*CTEST_FAILURE),
        )
        self.assertEqual(
            (verdict.kind, verdict.owner, verdict.tests),
            ("test", "pr", ("py3_init-densities_mpi_n_2",)),
        )

    def test_container_warning_does_not_hide_a_compile_error(self):
        verdict = triage.classify(
            build("Exit code 1 (Step: build (Command Line))"),
            log(CONTAINER_WARNING, *CLANG_ERROR),
        )
        self.assertEqual((verdict.kind, verdict.owner), ("compile", "pr"))

    def test_failed_image_pull_is_the_cis(self):
        verdict = triage.classify(
            build("Exit code 1 (Step: configure + build (Command Line))"),
            log(CONTAINER_WARNING, PULL_FAILED),
        )
        self.assertEqual((verdict.kind, verdict.owner), ("infra", "ci"))

    def test_error_only_in_third_party_code_is_not_blamed_on_the_pr(self):
        line = "[16:03:42] :\t [Step 3/4] /opt/buildagent/work/3bd/subprojects/samrai/source/SAMRAI/hier/Patch.h:14:10: fatal error: 'SAMRAI/SAMRAI_config.h' file not found"
        verdict = triage.classify(
            build("Exit code 1 (Step: build (Command Line))"), log(line)
        )
        self.assertEqual((verdict.kind, verdict.owner), ("compile", "unknown"))

    def test_timeout(self):
        verdict = triage.classify(
            build("Execution timeout (new); tests passed: 102"),
            log("[21:57:36]E:\t [Step 1/3] Execution timeout"),
        )
        self.assertEqual((verdict.kind, verdict.owner), ("timeout", "ci"))

    def test_report_parse_failure_with_no_failed_test_is_a_false_red(self):
        verdict = triage.classify(
            build("Tests passed: 1610; error message: Failed to parse 1 report (new)"),
            log("[10:00:00] :\t [Step 6/14] 100% tests passed"),
        )
        self.assertEqual((verdict.kind, verdict.owner), ("false_red", "ci"))

    def test_phlop_failure(self):
        line = "[03:25:43] :\t [Step 1/1]   /usr/lib64/openmpi/bin/mpirun -n 4 --bind-to none python3 -um test_divb_refinement DivBRefinementTest.test_x FAILED: exited 1"
        verdict = triage.classify(
            build("Exit code 1 (Step: run (Command Line))"), log(line)
        )
        self.assertEqual(
            verdict.tests, ("test_divb_refinement DivBRefinementTest.test_x",)
        )


class Excerpt(unittest.TestCase):
    def test_compile_excerpt_is_repo_relative_without_commands_and_stops_at_next_file(
        self,
    ):
        text = excerpt.excerpt(triage.Verdict("compile", "pr", ""), log(*CLANG_ERROR))
        self.assertIn("src/simulator/simulator.hpp:395:20: error:", text)
        self.assertNotIn("/opt/buildagent", text)
        self.assertNotIn("ccache", text)
        self.assertNotIn("cpp_1_1_3", text.split("First failure:")[1])

    def test_test_excerpt_is_the_failed_tests_output_only(self):
        lines = log(*CTEST_FAILURE)
        text = excerpt.excerpt(
            triage.classify(build("Tests passed: 1; exit code 8"), lines), lines
        )
        self.assertIn("AssertionError: densities differ", text)
        self.assertNotIn("py3_init-2d-fields", text)

    def test_cap_bounds_length(self):
        self.assertLessEqual(
            len(excerpt.cap(["x" * 200] * 1000)), excerpt.MAX_CHARS + 50
        )


REPORT = """<<<REPORT
Owner: pr
What: src/simulator/simulator.hpp:395, `already_setup` is set but never read.
Why: clang's -Werror turns this warning into an error; gcc doesn't warn.
Fix prompt:
> Remove the unused variable at src/simulator/simulator.hpp:395, keep behaviour.
REPORT>>>"""


class Report(unittest.TestCase):
    def test_last_complete_block_wins(self):
        echoed = (
            "the prompt shows <<<REPORT\nOwner: pr | ci | unsure\nREPORT>>> as format\n"
        )
        self.assertTrue(agent.parse(echoed + REPORT).startswith("Owner: pr"))

    def test_missing_field_is_rejected(self):
        self.assertIsNone(agent.parse(REPORT.replace("Why:", "Because:")))
        self.assertIsNone(agent.parse("no markers at all"))

    def test_context_lists_other_builds_of_the_commit_with_their_flags(self):
        red = {**build("x"), "id": 1, "config": "Phare_Phare_BuildGithubPrGccSs"}
        siblings = [
            red,
            {
                "id": 2,
                "config": "Phare_Phare_BuildGithubPr",
                "status": "SUCCESS",
                "state": "finished",
            },
            {
                "id": 3,
                "config": "Phare_Next_PrChecksV1",
                "status": "FAILURE",
                "state": "finished",
            },
            {
                "id": 4,
                "config": "Phare_Phare_BuildGithubPrClang",
                "status": None,
                "state": "running",
            },
        ]
        text = agent.siblings_text(red, siblings)
        self.assertIn("Phare_Phare_BuildGithubPr: SUCCESS (gcc, Debug", text)
        self.assertIn("Phare_Phare_BuildGithubPrClang: running", text)
        self.assertNotIn("GccSs", text)  # the red build itself
        self.assertNotIn("PrChecksV1:", text)  # the pipeline's umbrella build

    def test_sandbox_denies_every_writing_or_network_tool(self):
        permission = agent.SANDBOX["permission"]
        for tool in (
            "edit",
            "bash",
            "webfetch",
            "websearch",
            "task",
            "external_directory",
        ):
            self.assertEqual(permission[tool], "deny", tool)


class Comment(unittest.TestCase):
    def test_clean_removes_mentions_links_html_and_owner_line(self):
        text = comment.clean(
            "Owner: pr\nWhat: ping @someone, see [here](https://evil.example) or https://x.y <img src=x>"
        )
        self.assertNotIn("@s", text)
        self.assertNotIn("http", text)
        self.assertNotIn("<img", text)
        self.assertNotIn("Owner:", text)
        self.assertIn("here", text)

    def test_render_has_marker_counts_and_fix(self):
        report = agent.parse(REPORT)
        others = [
            (
                build("Execution timeout"),
                triage.Verdict("timeout", "ci", "hit the build time limit"),
            )
        ]
        text = comment.render(
            1,
            "a" * 40,
            [([build("x")], triage.Verdict("compile", "pr", "compile error"), report)],
            others,
            "Written by m.",
        )
        self.assertTrue(text.startswith(comment.marker(1)))
        self.assertIn("2 red**, 1 from this PR, 1 not", text)
        self.assertIn("**Fix**, to paste into your coding agent:", text)
        self.assertIn("hit the build time limit", text)

    def test_failure_not_run_says_so(self):
        text = comment.render(
            1,
            "a" * 40,
            [
                (
                    [build("x")],
                    triage.Verdict("test", "pr", "failed test"),
                    comment.SKIPPED,
                )
            ],
            [],
            "",
        )
        self.assertIn("not explained: too many", text)


def status(context, state, build_id=1):
    url = f"https://hephaistos.lpp.polytechnique.fr/teamcity/buildConfiguration/Phare_Phare_BuildGithubPr/{build_id}"
    return {"context": context, "state": state, "target_url": url}


class Ci(unittest.TestCase):
    def test_waits_while_a_teamcity_build_is_pending(self):
        statuses = [
            status("gh_pr_gcc (Phare)", "failure", 100788),
            status("gh_pr_clang (Phare)", "pending", 100787),
        ]
        self.assertEqual(ci.decide(statuses), ("none", []))

    def test_red_once_all_done_other_tools_ignored(self):
        statuses = [
            status("gh_pr_gcc (Phare)", "failure", 100788),
            status("gh_pr_clang (Phare)", "success", 100787),
            status("gh_pr_gcc_samraisub (Phare)", "error", 100785),
            status("docs/readthedocs.org:phare", "pending"),
        ]
        self.assertEqual(ci.decide(statuses), ("red", [100785, 100788]))

    def test_green_and_nothing_from_teamcity(self):
        self.assertEqual(
            ci.decide([status("gh_pr_gcc (Phare)", "success")]), ("green", [])
        )
        self.assertEqual(
            ci.decide([status("docs/readthedocs.org:phare", "failure")]),
            ("none", []),
        )

    def test_build_id_from_both_link_shapes(self):
        self.assertEqual(
            ci.build_id(
                "https://hephaistos.lpp.polytechnique.fr/teamcity/build/100788"
            ),
            100788,
        )
        self.assertEqual(ci.build_id(status("x", "failure", 42)["target_url"]), 42)
        self.assertIsNone(ci.build_id(None))

    def test_a_commit_already_reported_is_not_redone(self):
        sha = "b" * 40
        red = comment.render(1, sha, [], [], "")
        green = comment.render_green(1, sha)
        self.assertTrue(ci.done(red, sha, "red"))
        self.assertFalse(ci.done(red, sha, "green"))
        self.assertTrue(ci.done(green, sha, "green"))
        self.assertFalse(ci.done(green, sha, "red"))
        self.assertFalse(ci.done(red, "c" * 40, "red"))


if __name__ == "__main__":
    unittest.main()
