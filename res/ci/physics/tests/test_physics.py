"""Nightly physics report scripts: what each must produce, and what it must survive.

run: python3 -m unittest discover -s res/ci/physics/tests
"""

import json
import os
import shutil
import stat
import sys
import tempfile
import textwrap
import unittest
import unittest.mock
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import decks
import extract
import metrics


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


class MetricsTest(unittest.TestCase):
    def test_growth_rate_summary(self):
        s = metrics.growth_rate_summary(0.1, theory=0.09, tol=0.02)
        # |fft[1]| of a circularly polarised mode carries no frequency: no omega
        self.assertEqual(set(s), {"gamma", "abs_error", "within_deck_tolerance"})
        self.assertAlmostEqual(s["abs_error"], 0.01)
        self.assertTrue(s["within_deck_tolerance"])
        self.assertFalse(
            metrics.growth_rate_summary(0.2, 0.09, 0.02)["within_deck_tolerance"]
        )

    def test_exact_shift_gives_zero_error_and_wrong_shift_does_not(self):
        L = 200.0
        x = np.linspace(0, L, 400, endpoint=False)
        f0 = np.tanh((x - 50) / 2) - np.tanh((x - 150) / 2)
        f1 = np.interp(np.mod(x + 20.0, L), x, f0, period=L)  # moved by -20
        self.assertLess(metrics.shifted_profile_error(x, f0, f1, -20.0, L), 1e-2)
        self.assertGreater(metrics.shifted_profile_error(x, f0, f1, +20.0, L), 0.5)


class ExtractTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root, self.out = Path(tmp.name, "root"), Path(tmp.name, "out")
        self.plots = self.root / "build/ot/plots/4"
        self.plots.mkdir(parents=True)

    def deck(self, **kw):
        d = {
            "name": "orszag_tang",
            "source": "src.py",
            "workdir": "build/ot",
            "figures": ["plots/*/divb_t*.png"],
            "extractor": "figures_only",
            "theory": {"quantity": "div B", "value": 0.0, "deck_tolerance": None},
        }
        d.update(kw)
        return d

    def run_deck(self, deck=None, out="1"):
        return extract.run_deck(deck or self.deck(), self.root, self.out / out)

    def test_figures_are_copied_with_the_deck_prefix_unless_unsafe(self):
        secret = self.root.parent / "secret"
        secret.write_text("not public")
        names = (
            "divb_t0.png",
            "divb_t1](evil).png",
            "divb_t2 x.png",
            "divb_t3\n### x.png",
            "divb_t4.png",
            "divb_t5.png",
            "rho_t0.png",  # not a pattern of the deck
        )
        for n in names:
            (self.plots / n).write_bytes(b"png")
        (self.plots / "divb_t1.png").symlink_to(secret)
        r = self.run_deck(self.deck(max_figures=2))
        self.assertEqual(r["status"], "ok")
        self.assertEqual(
            r["figures"], ["orszag_tang__divb_t0.png", "orszag_tang__divb_t4.png"]
        )
        self.assertEqual(
            sorted(p.name for p in (self.out / "1/figures").iterdir()), r["figures"]
        )
        self.assertEqual(
            r["skipped"],
            [
                "divb_t1.png: symlink",
                "divb_t1](evil).png: unsafe file name",
                "divb_t2 x.png: unsafe file name",
                "divb_t3\n### x.png: unsafe file name",
                "divb_t5.png: over max_figures 2",
            ],
        )

    def test_orszag_tang_keeps_a_mid_time_and_the_last_one(self):
        """The deck plots every time of linspace(0, 1, 6); with an adaptive dt a dump for T is
        written at the first step past T - dt, so its name can be just below T."""
        deck = decks.job(decks.load(), "mhd")[1]
        p = self.root / deck["workdir"] / "phare_outputs/orszag_tang_plots/1"
        adaptive = ("0.0", "0.2", "0.3987654321", "0.5987654321", "0.8", "0.9987654321")
        exact = ("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")
        for times, kept in (
            (adaptive, ("0.5987654321", "0.9987654321")),
            (exact, ("0.6", "1.0")),
        ):
            shutil.rmtree(p, ignore_errors=True)
            p.mkdir(parents=True)
            for t in times:
                for q in ("divb", "rho", "Ranks", "vx", "bx", "p"):
                    (p / f"orszag_tang_{q}_t{t}.png").write_bytes(b"png")
            r = self.run_deck(deck, out=times[-1])
            want = [
                f"orszag_tang__orszag_tang_{q}_t{t}.png"
                for q in ("divb", "rho")
                for t in kept
            ]
            self.assertEqual(sorted(r["figures"]), sorted(want))

    def test_outputs_older_than_t0_are_ignored(self):
        old = self.plots / "divb_t0.png"
        old.write_bytes(b"old")
        os.utime(old, (1000, 1000))
        r = self.run_deck()  # no .t0: freshness unknown, the old figure still counts
        self.assertEqual(
            (r["freshness"], r["status"], r["figures"]),
            ("unknown", "ok", ["orszag_tang__divb_t0.png"]),
        )
        (self.root / ".t0").touch()
        new = self.plots / "divb_t1.png"
        new.write_bytes(b"new")
        r = self.run_deck(out="2")
        self.assertEqual(
            (r["freshness"], r["status"], r["figures"], r["skipped"]),
            ("ok", "ok", ["orszag_tang__divb_t1.png"], ["divb_t0.png: older than .t0"]),
        )
        os.utime(new, (1000, 1000))
        r = self.run_deck(out="3")
        self.assertEqual((r["status"], r["figures"]), ("stale output ignored", []))

    def test_stale_h5_is_not_read(self):
        run = self.root / "build/ib/ion_ion_beam1d"
        run.mkdir(parents=True)
        (run / "EM_B.h5").write_bytes(b"h5")
        os.utime(run / "EM_B.h5", (1000, 1000))
        (self.root / ".t0").touch()
        deck = self.deck(
            name="ion_ion_beam1d",
            workdir="build/ib",
            figures=[],
            extractor="ion_ion_beam",
            source="does/not/exist.py",  # reading the h5 would be an extract error
        )
        r = self.run_deck(deck)
        self.assertEqual(r["status"], "stale output ignored")
        self.assertEqual(r["numbers"], {})
        self.assertEqual(r["skipped"], ["ion_ion_beam1d/EM_B.h5: older than .t0"])

    def test_a_deck_failure_is_a_status_not_an_exception(self):
        r = self.run_deck(self.deck(workdir="build/none"))
        self.assertEqual(r["status"], "output missing")
        r = self.run_deck(
            self.deck(extractor="ion_ion_beam", source="does/not/exist.py")
        )
        self.assertTrue(r["status"].startswith("extract error:"), r["status"])
        # h5py's message holds the whole path over a line: only the file's name is kept
        h5 = FileNotFoundError(
            2, "Unable to open file (name = '/tmp/x/run/EM_B.h5', errno = 2"
        )
        bad = {"boom": unittest.mock.Mock(side_effect=h5)}
        with unittest.mock.patch.dict(extract.NUMBERS, bad):
            r = self.run_deck(self.deck(extractor="boom"))
        self.assertEqual(r["status"], "extract error: EM_B.h5 is missing")
        other = {
            "boom": unittest.mock.Mock(side_effect=ValueError("bad /tmp/x/run/a.h5"))
        }
        with unittest.mock.patch.dict(extract.NUMBERS, other):
            r = self.run_deck(self.deck(extractor="boom"))
        self.assertEqual(r["status"], "extract error: ValueError: bad a.h5")

    def test_main_writes_the_summary(self):
        self.out.mkdir()
        (self.out / "decks_rc").write_text("3\n")

        def summary():
            extract.main(
                ["--job", "mhd", "--root", str(self.root), "--out", str(self.out)]
            )
            return json.loads((self.out / "summary.json").read_text())

        s = summary()
        self.assertEqual(s["decks_exit_code"], 3)
        none = {"convergence": None, "orszag_tang": None}
        self.assertEqual({d["name"]: d["deck_passed"] for d in s["decks"]}, none)
        # deck_passed: a boolean from decks_results.json, anything else is unknown
        for results, want in (
            (
                {"convergence": False, "orszag_tang": "yes"},
                none | {"convergence": False},
            ),
            ({"orszag_tang": True}, none | {"orszag_tang": True}),
            ("{not json", none),
        ):
            text = results if isinstance(results, str) else json.dumps(results)
            (self.out / "decks_results.json").write_text(text)
            passed = {d["name"]: d["deck_passed"] for d in summary()["decks"]}
            self.assertEqual(passed, want)

    def test_find_logs(self):
        phlop = self.root / ".phlop/build"
        for f in (
            "ot/orszag_tang/OrszagTangTest/test_run.stdout",
            "other/orszag_tang.stdout",
            # a plain script's logs, named after its directory (nightly 101183, ionIonBeam)
            "tests/functional/ionIonBeam.stdout",
            "tests/functional/ionIonBeam.stderr",
            "tests/functional/td.stdout",
        ):
            (phlop / f).parent.mkdir(parents=True, exist_ok=True)
            (phlop / f).write_text("x")

        def names(workdir, stem):
            return [p.name for p in extract.find_logs(self.root, workdir, stem)]

        self.assertEqual(names("build/ot", "orszag_tang"), ["test_run.stdout"])
        self.assertEqual(
            names("build/tests/functional/ionIonBeam", "ion_ion_beam1d"),
            ["ionIonBeam.stdout", "ionIonBeam.stderr"],
        )


if __name__ == "__main__":
    unittest.main()
