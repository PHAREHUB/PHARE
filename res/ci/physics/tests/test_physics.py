"""Nightly physics report scripts: what each must produce, and what it must survive.

run: python3 -m unittest discover -s res/ci/physics/tests
"""

import os
import stat
import sys
import tempfile
import textwrap
import unittest
import unittest.mock
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import decks


def script(path, body):
    """An executable Python script at path, returned as a str."""
    p = Path(path)
    p.write_text("#!/usr/bin/env python3\n" + body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return str(p)


class DecksTest(unittest.TestCase):
    CFG = """\
        model: m
        opencode_version: 1
        {page}
        jobs:
          hybrid:
            - {{name: x, source: s, workdir: w, figures: [], extractor: {extractor},
               theory: {{quantity: q, value: 1, deck_tolerance: null}}{extra}}}
        """

    def load(
        self, extra="", page="page: {branch: b, dots: 5}", extractor="figures_only"
    ):
        with tempfile.NamedTemporaryFile("w", suffix=".yml", delete=False) as f:
            f.write(
                textwrap.dedent(self.CFG).format(
                    page=page, extra=extra, extractor=extractor
                )
            )
        self.addCleanup(os.unlink, f.name)
        return decks.load(f.name)

    def test_repo_config_loads(self):
        cfg = decks.load()
        self.assertEqual(
            {j: [d["name"] for d in ds] for j, ds in cfg["jobs"].items()},
            {
                "hybrid": ["ion_ion_beam1d", "td1d"],
                "mhd": ["convergence", "orszag_tang"],
            },
        )
        self.assertIs(decks.job(cfg, "hybrid")[1]["judge"], False)  # td1d: vth = T(x)

    def test_bad_config_is_refused(self):
        deck = decks.job(self.load(", max_figures: 2, key_number: gamma"), "hybrid")[0]
        self.assertEqual((deck["max_figures"], deck["key_number"]), (2, "gamma"))
        extras = (
            ", max_figures: two",
            ", max_figures: 0",
            ", max_figures: true",  # a YAML true is an int to isinstance
            ", own_check: maybe",
            ", judge: no_thanks, judge_reason: r",
            ", judge: false",
            ", judge: false, judge_reason: ''",
            ", judge: false, judge_reason: 3",
            ", key_number: 3",
            ", key_number: ''",
        )
        # the page's branch reaches git: nothing that could pass for an option or a path
        pages = (
            "",
            "page: {branch: b, dots: 5, repo: x}",
            "page: {branch: b, dots: x}",
            "page: {branch: b, dots: 0}",
            "page: {branch: b, dots: true}",
            "page: {branch: '', dots: 5}",
            "page: {branch: -f, dots: 5}",
            "page: {branch: 'a b', dots: 5}",
            "page: {branch: 'a;rm', dots: 5}",
        )
        bad = (
            [{"extractor": "nope"}]
            + [{"extra": x} for x in extras]
            + [{"page": p} for p in pages]
        )
        for kw in bad:
            with self.subTest(**kw), self.assertRaises(ValueError):
                self.load(**kw)


if __name__ == "__main__":
    unittest.main()
