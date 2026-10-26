"""Sorting, excerpts, report parsing and comment cleaning, on snippets of real TeamCity logs.

run: python3 -m unittest discover -s res/ci/explainer/tests
"""

import contextlib
import io
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import agent
import ci
import comment
import excerpt
import explain
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

    def test_exit_code_with_no_failed_test_found_goes_to_the_agent(self):
        verdict = triage.classify(
            build(
                "Tests passed: 1621; exit code 8 (Step: run mpi test (Command Line))"
            ),
            log("[10:00:00] :\t [Step 6/14] something the rules don't parse"),
        )
        self.assertEqual((verdict.kind, verdict.owner), ("unknown", "unknown"))

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
    # tricks that got past escaping the answer as markdown, each checked on GitHub's renderer
    TRICKS = (
        "What: ping @someone, see [here](https://evil.example) <img src=x> &#64;octocat",
        '\\`<a href="//evil.example">x</a>`',
        '```<a href="//evil.example">x</a>`',
        '`<a href="//evil.example">x</a> http://y`',
        'What: x\n```\ncode\n````\n<a href="//evil.example/a">x</a>\n```\n',
        '- item\n  ```\n<a href="//evil.example/b">x</a>\n  ```\n',
        "What: see [a\\]b](//evil.example/c) and [a [b] c](//evil.example/d)",
        "> [r]: //evil.example/e\n\nWhat: [QUOTEREF][r] HTTPS://evil.example/f",
        "`" * 7 + "\n<b>x</b>",
        "x" * (comment.MAX_REPORT - 10) + '`<a href="//evil.example">x</a>`',
    )

    def test_the_answer_is_one_code_block_it_cannot_close(self):
        for trick in self.TRICKS:
            text = comment.as_code(trick)
            fence = text.split("\n", 1)[0].removesuffix("text")
            self.assertTrue(set(fence) == {"`"} and len(fence) >= 3, trick)
            self.assertTrue(text.endswith("\n" + fence), trick)
            body = text[len(fence) + 5 : -len(fence) - 1]
            self.assertFalse(re.search("`" * len(fence), body), trick)

    def test_code_block_keeps_cpp_drops_owner_and_wraps(self):
        text = comment.as_code(
            "Owner: pr\nWhat: `std::vector<int>` and #include <cstdint>\nFix prompt:\n> "
            + "word " * 40
        )
        self.assertIn("`std::vector<int>` and #include <cstdint>", text)
        self.assertNotIn("Owner:", text)
        self.assertNotIn("> word", text)
        self.assertTrue(all(len(line) <= comment.WIDTH for line in text.splitlines()))

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
        self.assertIn("2 red**: 1 explained below, 1 put on the CI by the rules.", text)
        self.assertIn("```text\nWhat: src/simulator/simulator.hpp:395", text)
        self.assertIn("\nRemove the unused variable", text)
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
        self.assertIn("1 red**: 1 not explained.", text)

    def test_failures_without_a_test_name_are_not_merged(self):
        unknown = triage.Verdict("unknown", "unknown", "red")
        self.assertNotEqual(
            explain.failure_key(unknown, "segfault in a"),
            explain.failure_key(unknown, "out of memory in b"),
        )
        compile_error = triage.Verdict("compile", "pr", "compile error")
        self.assertEqual(  # gcc and clang builds print the same file differently
            explain.failure_key(compile_error, "src/core/a.hpp:196:5: error: x"),
            explain.failure_key(compile_error, "../src/core/a.hpp:196:12: error: y"),
        )
        named = triage.Verdict("test", "pr", "test failure", ("test-gridlayout",))
        self.assertEqual(
            explain.failure_key(named, "gcc words it so"),
            explain.failure_key(named, "clang words it otherwise"),
        )

    def test_the_agent_gets_no_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            repo.mkdir()
            (repo / "a.txt").write_text("a")
            (repo / "leak").symlink_to("/etc/hostname")
            (repo / "up").symlink_to("..")

            def git(*args):
                subprocess.run(
                    ["git", "-C", str(repo), *args], check=True, capture_output=True
                )

            git("init", "-q")
            git("add", ".")
            git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "x")
            dest = Path(tmp) / "src"
            explain.snapshot(repo, "HEAD")(dest)
            self.assertEqual(sorted(p.name for p in dest.iterdir()), ["a.txt"])


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

    def test_a_new_red_build_on_the_same_commit_is_redone(self):
        sha = "b" * 40
        body = comment.render(
            1,
            sha,
            [([build("x")], triage.Verdict("compile", "pr", "compile error"), None)],
            [],
            "",
        )
        self.assertTrue(ci.done(body, sha, "red", [1]))
        self.assertFalse(ci.done(body, sha, "red", [1, 2]))

    def test_a_red_build_passing_on_rerun_is_redone(self):
        sha = "b" * 40
        two = [dict(build("x"), id=101), dict(build("x"), id=102)]
        body = comment.render(
            1,
            sha,
            [(two, triage.Verdict("compile", "pr", "compile error"), None)],
            [],
            "",
        )
        self.assertTrue(ci.done(body, sha, "red", [101, 102]))
        self.assertFalse(ci.done(body, sha, "red", [102]))
        self.assertFalse(ci.done("**CI on `bbbbbbb`: 1 red** (102)", sha, "red", [102]))

    def test_nothing_posted_once_the_pr_moved_on(self):
        calls = []

        def api(path, method="GET", body=None):
            calls.append(method)
            return {"state": "open", "head": {"sha": "c" * 40}, "labels": []}

        real, ci.api = ci.api, api
        try:
            ci.post(1, "b" * 40, "text")
        finally:
            ci.api = real
        self.assertEqual(calls, ["GET"])

    def test_posted_while_the_commit_is_the_head(self):
        calls = []

        def api(path, method="GET", body=None):
            calls.append((method, path))
            return (
                []
                if "comments" in path
                else {"state": "open", "head": {"sha": "b" * 40}, "labels": []}
            )

        real, ci.api = ci.api, api
        try:
            ci.post(1, "b" * 40, "text")
        finally:
            ci.api = real
        self.assertEqual(calls[-1], ("POST", "issues/1/comments"))

    def test_no_comment_when_the_agent_answered_nothing_usable(self):
        red = (build("x"), triage.Verdict("unknown", "unknown", "red"), "excerpt")
        real = ci.tc.build, ci.tc.builds_of, ci.explain.sort, ci.explain.explain
        ci.tc.build = lambda i: build("x")
        ci.tc.builds_of = lambda sha: []
        ci.explain.sort = lambda builds, log: ([], {"k": [red]}, {})
        ci.explain.explain = lambda *args, **kwargs: None  # e.g. a mistyped model name
        try:
            with (
                tempfile.TemporaryDirectory() as out,
                contextlib.redirect_stdout(io.StringIO()),
            ):
                ci.explain_commit(1, "0" * 40, [1], Path(out))
                self.assertFalse((Path(out) / "comment.md").exists())
        finally:
            ci.tc.build, ci.tc.builds_of, ci.explain.sort, ci.explain.explain = real

    def test_no_comment_when_every_red_was_canceled(self):
        statuses = [status("gh_pr_gcc (Phare)", "failure", 7)]
        real = ci.api, ci.tc.build
        ci.api = lambda path, *_: {"statuses": statuses}
        ci.tc.build = lambda i: {
            "id": i,
            "pr": 1,
            "status": "UNKNOWN",
            "status_text": "Canceled",
        }
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                ci.check("b" * 40)
        finally:
            ci.api, ci.tc.build = real
        self.assertIn("action=none", out.getvalue())


if __name__ == "__main__":
    unittest.main()
