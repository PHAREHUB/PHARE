"""Nightly physics report scripts: what each must produce, and what it must survive.

run: python3 -m unittest discover -s res/ci/physics/tests
"""

import base64
import contextlib
import io
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
import unittest.mock
import urllib.error
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import decks
import extract
import history
import judge
import metrics
import page
import publish
import tc


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


class HistoryTest(unittest.TestCase):
    SHA = "a" * 40

    def part(
        self, status="ok", agent="text", rc=0, passed=None, dstatus="ok", numbers=None
    ):
        deck = {
            "name": "d",
            "figures": [],
            "deck_passed": passed,
            "own_check": None,
            "numbers": numbers or {},
            "status": dstatus,
            "freshness": "fresh",
        }
        return {
            "job": "hybrid",
            "meta": {"agent_status": status, "agent_reason": "quota", "sha": self.SHA},
            "summary": {"decks_exit_code": rc, "decks": [deck]},
            "agent": agent,
        }

    def rec(self, part=None, cfg=None, model=("looks right", "fine")):
        part = part or self.part()
        deck = (part.get("summary") or {}).get("decks", [None])[0]
        return history.test_record(deck, cfg or {}, part, model)

    def hist(self, *statuses, passed=None, manual=False):
        out = []
        for i, s in enumerate(statuses, 1):
            rec = {"status": s, "job": "hybrid", "deck_passed": passed, "key": None}
            run = {"sha": self.SHA, "manual": manual, "time": "01:00", "builds": {}}
            out.append(history.night(f"2026-10-0{i}", run, {"t": rec}))
        return out

    def test_status_order(self):
        """The first matching rule wins (history.test_record)."""
        ok = ("looks right", "fine")
        for part, cfg, model, want in (
            (
                self.part(passed=False, dstatus="output missing"),
                {},
                ok,
                ("own check failed", None),
            ),
            (
                {"job": "hybrid", "missing": True},
                {},
                ok,
                ("not run", "the hybrid job left no results tonight"),
            ),
            (
                self.part(dstatus="output missing"),
                {},
                ok,
                ("not run", "output missing"),
            ),
            (
                self.part(dstatus="extract error <b>"),
                {},
                ok,
                ("no data", "extract error <b>"),
            ),
            (
                self.part(dstatus="stale output ignored"),
                {},
                ok,
                ("no data", "stale output ignored"),
            ),
            (
                self.part(),
                {"judge": False, "judge_reason": "why"},
                ok,
                ("not judged", "why"),
            ),
            (self.part(status="unavailable"), {}, ok, ("no analysis", "quota")),
            (self.part(agent=""), {}, ok, ("no analysis", "quota")),
            (self.part(), {}, ("suspicious", "x"), ("suspicious", None)),
            (self.part(), {}, None, ("cannot tell", None)),
            (self.part(), {}, ("bogus", "x"), ("cannot tell", None)),
            (self.part(), {}, ("own check failed", "x"), ("cannot tell", None)),
            (self.part(rc=1, passed=True), {}, ok, ("looks right", None)),
        ):
            r = self.rec(part, cfg, model)
            self.assertEqual((r["status"], r["reason"]), want)
        self.assertEqual(
            history.test_record(None, {}, self.part(), ok)["status"], "not run"
        )

    def test_a_failed_tests_step_is_never_green(self):
        for model in (("looks right", ""), ("cannot tell", ""), None, ("bogus", "")):
            r = self.rec(self.part(rc=137), model=model)
            self.assertEqual(
                (r["status"], r["reason"]),
                ("own check unknown", "tests step failed (exit 137)"),
                model,
            )
        r = self.rec(self.part(rc=137), model=("suspicious", "x"))
        self.assertEqual(
            (r["status"], r.get("own_check_unknown")), ("suspicious", True)
        )
        r = self.rec(self.part(status="unavailable", rc=137))
        self.assertEqual(
            (r["status"], r.get("own_check_unknown")), ("no analysis", True)
        )
        self.assertNotIn("own_check_unknown", self.rec(self.part(status="unavailable")))

    def test_since_new_tonight_and_prev_sha(self):
        for statuses, since, new in (
            (("looks right", "suspicious", "suspicious"), "2026-10-02", False),
            (("looks right", "suspicious"), "2026-10-02", True),
            (("suspicious",), "2026-10-01", False),  # a first night is not "new"
            # neutral nights neither break nor extend a run
            (
                ("looks right", "suspicious", "not run", "suspicious"),
                "2026-10-02",
                False,
            ),
            (("looks right", "no analysis", "suspicious"), "2026-10-03", True),
        ):
            s = history.state(self.hist(*statuses), "t", True)
            self.assertEqual(
                (s["abnormal"], s["since"], s["new_tonight"]),
                (True, since, new),
                statuses,
            )
        h = self.hist("looks right", "looks right", "suspicious")
        self.assertEqual(history.state(h, "t", True)["prev_sha"], self.SHA)
        h[0]["run"]["sha"], h[1]["run"]["sha"] = "b" * 40, "c" * 40
        self.assertEqual(history.state(h, "t", True)["prev_sha"], "c" * 40)
        h[1]["run"]["sha"] = None
        # the newest earlier night that has a sha
        self.assertEqual(history.state(h, "t", True)["prev_sha"], "b" * 40)
        h[0]["run"]["sha"] = None
        self.assertIsNone(history.state(h, "t", True)["prev_sha"])

    def test_carry_over_when_the_model_is_down(self):
        s = history.state(
            self.hist("suspicious", "no analysis", "no analysis"), "t", True
        )
        self.assertEqual(
            (s["abnormal"], s["since"], s["last_judged"], s["nights_since_judged"]),
            (True, "2026-10-01", "2026-10-01", 2),
        )

    def test_judge_false_rules(self):
        # not judged tonight: a run the model started is dropped...
        self.assertFalse(
            history.state(self.hist("suspicious", "not judged"), "t", False)["abnormal"]
        )
        self.assertTrue(
            history.state(self.hist("suspicious", "not judged"), "t", True)["abnormal"]
        )
        # ...and a run of failed own checks ends on a passing own check
        h = self.hist("own check failed", "not judged")
        h[1]["tests"]["t"]["deck_passed"] = True
        self.assertFalse(history.state(h, "t", False)["abnormal"])
        h = self.hist("suspicious", "no analysis")
        h[1]["tests"]["t"]["deck_passed"] = True
        self.assertTrue(history.state(h, "t", True)["abnormal"])

    def test_one_night_per_date(self):
        def n(date, manual, tag=None):
            return history.night(date, {"sha": tag, "manual": manual}, {})

        for old, new, kept, recorded in (
            (
                False,
                True,
                "old",
                False,
            ),  # a manual night never replaces a scheduled one
            (True, False, "new", True),
            (True, True, "new", True),
            (False, False, "new", True),
        ):
            h, ok = history.add_night(
                [n("2026-10-01", old, "old")], n("2026-10-01", new, "new")
            )
            self.assertEqual((len(h), h[0]["run"]["sha"], ok), (1, kept, recorded))
        h = history.merge_pending(
            [n("2026-10-01", False)],
            [
                n("2026-10-01", True),
                n("2026-10-03", True),
                n("2026-10-02", False),
                n("2026-10-03", False),
            ],
        )
        self.assertEqual(
            [(x["date"], x["run"]["manual"]) for x in h],
            [("2026-10-01", False), ("2026-10-02", False), ("2026-10-03", False)],
        )

    def test_key_number(self):
        cfg = {"key_number": "gamma"}
        for numbers, want in (
            ({"gamma": 0.0875}, 0.0875),
            ({"gamma": 1}, 1.0),
            ({}, None),
            ({"gamma": True}, None),
            ({"gamma": "1"}, None),
            ({"gamma": 10**400}, None),
        ):
            r = self.rec(self.part(numbers=numbers), cfg)
            self.assertEqual(r["key"], want)
            if want is not None:
                self.assertIs(type(r["key"]), float)

    def expiry_history(self, n=3):
        out = []
        for i in range(n):
            builds = {
                "hybrid": {"id": 100 + i, "type": "T_hybrid", "start": "2026-10-01"}
            }
            run = {"sha": self.SHA, "manual": False, "time": "01:00", "builds": builds}
            rec = {"status": "looks right", "job": "hybrid", "key": None}
            out.append(history.night(f"2026-10-{i + 1:02d}", run, {"t": rec}))
        return out

    def test_only_a_404_expires_a_night(self):
        h = self.expiry_history(2)
        h2, _ = history.expire(h, lambda b, j: None, True, 90)  # network error: unknown
        self.assertFalse(any("expired" in n["tests"]["t"] for n in h2))
        h3, _ = history.expire(h, lambda b, j: False, True, 90)
        self.assertTrue(all(n["tests"]["t"]["expired"] is True for n in h3))
        self.assertNotIn("expired", h[0]["tests"]["t"])  # the input is not changed
        calls = []
        history.expire(h3, lambda b, j: calls.append(b), True, 90)
        self.assertEqual(calls, [])  # an expired night is not asked again

    def test_canary_blocks_all_expiry(self):
        builds = [(100, "hybrid"), (101, "mhd")]
        self.assertFalse(history.canary(builds, lambda b, j: False))
        self.assertFalse(history.canary(builds, lambda b, j: None))
        self.assertTrue(history.canary(builds, lambda b, j: j == "mhd"))
        self.assertFalse(history.canary([], lambda b, j: True))
        h = self.expiry_history(2)
        calls = []
        h2, notes = history.expire(h, lambda b, j: calls.append(b) or False, False, 90)
        self.assertEqual((calls, h2, len(notes)), ([], h, 1))
        self.assertIsNot(h2, h)

    def test_what_is_asked(self):
        """Only the last `drawn` nights, each (build, job) once, never a job that is not a
        plain name (it would go into a URL)."""
        calls = []

        def exists(b, j):
            calls.append((b, j))
            return True

        history.expire(self.expiry_history(100), exists, True, 90)
        self.assertEqual(sorted(b for b, _ in calls), list(range(110, 200)))
        h = self.expiry_history(3)
        for n in h:
            n["run"]["builds"]["hybrid"]["id"] = 100  # a reused build
        h[0]["tests"]["u"] = {"status": "looks right", "job": "../x"}
        h[0]["run"]["builds"]["../x"] = {"id": 101, "type": "T", "start": "2026-10-01"}
        calls.clear()
        history.expire(h, exists, True, 90)
        self.assertEqual(calls, [(100, "hybrid")])

    def test_expire_stops_at_the_deadline(self):
        h = self.expiry_history(5)
        ticks = iter([0, 0, 1, 500, 500, 500])
        calls = []

        def exists(b, j):
            calls.append(b)
            return False

        _, notes = history.expire(h, exists, True, 90, 120, clock=lambda: next(ticks))
        self.assertLess(len(calls), 5)
        self.assertTrue(any("deadline" in n for n in notes))

    def test_load_validation(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "h.json"
            self.assertEqual(history.load(p), ([], [], True))  # no file yet
            good = self.hist("looks right", "suspicious", "suspicious")
            history.dump(good, p)
            self.assertEqual(history.load(p), (good, [], True))
            bad = dict(good[0], date=5)
            for text, kept in (
                ("{not json", 0),
                ("[]", 0),
                (json.dumps({"version": 2, "nights": []}), 0),
                (json.dumps({"version": 1, "nights": good + [bad]}), 3),
                (json.dumps({"version": 1, "nights": [bad, bad, good[0]]}), 0),
            ):
                with self.subTest(text=text[:60]):
                    p.write_text(text)
                    n, notes, clean = history.load(p)
                    self.assertEqual((len(n), clean, len(notes)), (kept, False, 1))
            for field, value in (
                ("job", ["x"]),
                ("job", ".."),
                ("job", "a/b"),
                ("job", 3),
                ("deck_passed", "yes"),
                ("deck_passed", 1),
                ("status", "s" * 41),
                ("key", True),
                ("key", "text"),
                ("key", 10**400),  # a huge int is a bad record, not an error
                ("reason", "r" * 201),
                ("expired", 1),
                ("own_check_unknown", "no"),
            ):
                with self.subTest(field=field, value=str(value)[:20]):
                    h = self.hist("looks right", "suspicious")
                    h[0]["tests"]["t"][field] = value
                    history.dump(h, p)
                    n, notes, clean = history.load(p)
                    self.assertFalse(clean)
                    self.assertEqual([len(x["tests"]) for x in n], [0, 1])
                    self.assertIn("1 bad record", notes[0])

    def test_unknown_status_is_kept(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "h.json"
            history.dump(self.hist("suspicious", "no output"), p)
            n, _, clean = history.load(p)
            self.assertTrue(clean)
            s = history.state(n, "t", True)
            self.assertEqual(s["status"], "unknown status no output")
            self.assertTrue(s["abnormal"])  # neutral: carries over the suspicious run


HYB_TYPE = "Phare_Next_NightlyV1gcc_virtual_hybrid_functional_V__1"
MHD_TYPE = "Phare_Next_NightlyV1gcc_virtual_mhd_functional_V__1"


class TcTest(unittest.TestCase):
    BASE = "http://tc"
    SHA = "b" * 40

    def fake(self, payload):
        seen = []

        def get(url, *a, **k):
            seen.append(url)
            return json.dumps(payload).encode()

        get.seen = seen
        return get

    def test_upstream_builds_parses_and_filters(self):
        def b(i, t, sha):
            return {
                "id": i,
                "buildTypeId": t,
                "startDate": "20261028T010012+0000",
                "defaultBranch": True,
                "triggered": {"type": "schedule"},
                "revisions": {"revision": [{"version": sha}]},
            }

        get = self.fake(
            {
                "build": [
                    b(1, HYB_TYPE, self.SHA),
                    b(2, "x/../y", self.SHA),
                    b(3, "T_mhd", "abc"),
                ]
            }
        )
        out = tc.upstream_builds(self.BASE, 7, get)
        self.assertEqual([x["id"] for x in out], [1, 3])
        self.assertEqual(out[0]["sha"], self.SHA)
        self.assertIsNone(out[1]["sha"])
        self.assertEqual(out[0]["start"].isoformat(), "2026-10-28T01:00:12+00:00")
        self.assertIs(out[0]["default_branch"], True)
        self.assertEqual(out[0]["triggered"], "schedule")
        self.assertTrue(
            get.seen[0].startswith(f"{self.BASE}/guestAuth/app/rest/builds?")
        )
        self.assertIn("snapshotDependency:(to:(id:7)),defaultFilter:false", get.seen[0])

    def test_head_of_picks_the_newest_head(self):
        get = self.fake(
            {
                "build": [
                    {"id": 101179, "buildTypeId": "Phare_Next_NightlyV1gcc"},
                    {"id": 101249, "buildTypeId": "Phare_Next_NightlyV1gcc"},
                    {"id": 101300, "buildTypeId": "Phare_Next_NightlyV1gcc_virtual_x"},
                ]
            }
        )
        self.assertEqual(tc.head_of(self.BASE, 101183, get), (101249, [101179]))
        self.assertIn("snapshotDependency:(from:(id:101183))", get.seen[0])

    def test_bad_answers_give_nothing(self):
        for payload in ([], [1], "x", {}, {"build": 3}, {"build": [1, "a"]}):
            get = self.fake(payload)
            self.assertEqual(tc.upstream_builds(self.BASE, 1, get), [])
            self.assertEqual(tc.head_of(self.BASE, 1, get), (None, []))
            self.assertIsNone(tc.build_info(self.BASE, 1, get))
        info = {"id": 5, "buildTypeId": "T_x", "startDate": "20261028T010012+0000"}
        b = tc.build_info(self.BASE, 5, self.fake(info))
        # absent is not the default branch, and not a schedule
        self.assertEqual(
            (b["default_branch"], b["sha"], b["triggered"]), (False, None, None)
        )
        for triggered in (5, "schedule", [1], {"type": 3}, {}):
            b = tc.build_info(self.BASE, 5, self.fake(info | {"triggered": triggered}))
            self.assertIsNone(b["triggered"], triggered)

    def test_physics_out_exists_codes(self):
        url = []

        def ok(u, timeout=30, headers=None):
            url.append((u, timeout))
            return b"[]"

        self.assertIs(tc.physics_out_exists(self.BASE, 9, "hybrid", ok), True)
        self.assertEqual(url[0][1], 10)
        self.assertTrue(
            url[0][0].endswith("/builds/id:9/artifacts/children/physics-out/hybrid")
        )

        def raiser(exc):
            def get(u, timeout=30, headers=None):
                raise exc

            return get

        def http(code):
            return urllib.error.HTTPError("u", code, "m", {}, None)

        self.assertIs(
            tc.physics_out_exists(self.BASE, 9, "h", raiser(http(404))), False
        )
        for exc in (http(500), urllib.error.URLError("x"), TimeoutError()):
            self.assertIsNone(tc.physics_out_exists(self.BASE, 9, "h", raiser(exc)))


AGENT_OK = """some chatter the agent added
<<<REPORT
### convergence
Verdict: looks right
Slopes 1.0, 2.0, 3.1, 5.0, 4.9 against expected 1, 2, 3, 5, 5.
REPORT>>>
more chatter
"""
AGENT_BLOCK = (
    "### convergence\nVerdict: looks right\n"
    "Slopes 1.0, 2.0, 3.1, 5.0, 4.9 against expected 1, 2, 3, 5, 5."
)


def fake_opencode(dirname, stdout, rc=0, version="1.18.32", stderr=""):
    return script(
        Path(dirname, "opencode"),
        f"import sys\n"
        f"if sys.argv[1:2] == ['--version']:\n    print({version!r}); sys.exit(0)\n"
        f"sys.stdout.write({stdout!r})\n"
        f"sys.stderr.write({stderr!r})\nsys.exit({rc})\n",
    )


class JudgeTest(unittest.TestCase):
    def work(self):
        """(a temporary dir, a work dir inside it holding PROMPT.md)."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        work = Path(tmp.name, "work")
        work.mkdir()
        (work / "PROMPT.md").write_text("x")
        return Path(tmp.name), work

    def test_version_is_one_word_from_the_version_output(self):
        for printed, want in (
            ("opencode 1.18.32\n", "1.18.32"),
            ("v1.19.0-beta.2 (abc)", "1.19.0-beta.2"),
            ("no version here", "unknown"),
            ("", "unknown"),
        ):
            self.assertEqual(judge.parse_version(printed), want, printed)

    def test_history_reads_finished_default_branch_nights_and_survives_errors(self):
        urls = []

        def get(url):
            urls.append(url)
            if "/app/rest/builds?" in url:
                return json.dumps({"build": [{"id": 9}, {"id": 8}, {"id": 7}]}).encode()
            if "/8:id/" in url:
                return json.dumps({"decks": []}).encode()
            raise OSError("404")

        h = judge.fetch_history("http://x", "bt", "mhd", 9, 7, get=get)
        # not 9 (tonight), not 7 (no summary.json)
        self.assertEqual([e["build_id"] for e in h], [8])
        for part in (
            "branch:(default:true)",
            "personal:false",
            "canceled:false",
            "failedToStart:false",
            "state:finished",
        ):
            self.assertIn(part, urls[0])

        def down(url):
            raise OSError("unreachable")

        self.assertEqual(
            judge.fetch_history("http://x", "bt", "mhd", 5, 7, get=down), []
        )

    def test_the_agent_runs_read_only_in_its_work_dir(self):
        """What opencode is given: its config file, its arguments, its environment."""
        tmp, work = self.work()
        exe = script(
            tmp / "opencode",
            "import json, os, sys\n"
            "config = json.load(open(os.environ['OPENCODE_CONFIG']))\n"
            "print(json.dumps({'env': dict(os.environ), 'argv': sys.argv[1:], "
            "'config': config}))\n",
        )

        def run():
            out, why = judge.run_agent(exe, work, "m/x", 30)
            self.assertIsNone(why)
            return json.loads(out)

        self.assertNotIn("-f", run()["argv"])  # no figures, no -f
        (work / "figures").mkdir()
        for n in ("d__b.png", "d__a.png", "d__c.mp4"):
            (work / "figures" / n).write_bytes(b"x")
        seen = run()
        argv, config = seen["argv"], seen["config"]
        self.assertEqual(argv[:2], ["run", "--pure"])
        self.assertEqual(argv[argv.index("-m") + 1], "m/x")
        self.assertEqual(argv[argv.index("--dir") + 1], str(work))
        self.assertLess(argv.index("x"), argv.index("-f"))  # the message comes first
        pngs = [str(work / "figures" / n) for n in ("d__a.png", "d__b.png")]
        self.assertEqual(argv[argv.index("-f") + 1 :], pngs)  # -f takes the rest
        self.assertEqual(config["model"], "m/x")
        self.assertIs(config["lsp"], False)
        self.assertIs(config["formatter"], False)
        perm = config["permission"]
        self.assertEqual(
            {k for k, v in perm.items() if v != "deny"},
            {"read", "glob", "grep", "list"},
        )
        for k in (
            "edit",
            "bash",
            "task",
            "external_directory",
            "webfetch",
            "websearch",
            "skill",
            "doom_loop",
        ):
            self.assertEqual(perm.get(k), "deny", k)

    def test_the_agent_gets_an_allow_listed_env_and_a_home_outside_its_work_dir(self):
        tmp, work = self.work()
        exe = script(
            tmp / "opencode", "import json, os\nprint(json.dumps(dict(os.environ)))\n"
        )
        dropped = {
            "SOME_PAT": "p",
            "GH_AUTH": "a",
            "DOCKER_CONFIG": "d",
            "OPENCODE_SERVER_PASSWORD": "o",  # OPENCODE_* only as judge sets them
            "SSL_KEY_PASSWORD": "k",  # SSL_* is not a wildcard
        }
        kept = {
            "https_proxy": "http://proxy:3128",
            "HTTPS_PROXY": "http://proxy:3128",
            "No_Proxy": "localhost",
            "all_proxy": "http://proxy:3128",
            "LC_ALL": "C.UTF-8",
            "TZ": "UTC",
            "SSL_CERT_FILE": "/etc/ssl/lab.pem",
            "SSL_CERT_DIR": "/etc/ssl/certs",
            "NODE_EXTRA_CA_CERTS": "/etc/ssl/lab.pem",
            "REQUESTS_CA_BUNDLE": "/etc/ssl/lab.pem",
            "CURL_CA_BUNDLE": "/etc/ssl/lab.pem",
        }
        with unittest.mock.patch.dict(os.environ, dropped | kept):
            out, why = judge.run_agent(exe, work, "m/x", 30)
        self.assertIsNone(why)
        env = json.loads(out)
        for k in dropped:
            self.assertNotIn(k, env)
        for k, v in kept.items():
            self.assertEqual(env.get(k), v, k)
        allowed = re.compile(
            r"PATH|LANG|TZ|LC_\w+|(?i:(https?|no|all)_proxy)|HOME|XDG_\w+|OPENCODE_\w+"
            r"|SSL_CERT_(FILE|DIR)|NODE_EXTRA_CA_CERTS|REQUESTS_CA_BUNDLE|CURL_CA_BUNDLE"
        )
        self.assertEqual([k for k in env if not allowed.fullmatch(k)], [])
        self.assertEqual(env["OPENCODE_CONFIG"], str(work / "opencode.json"))
        for k in ("HOME", "XDG_DATA_HOME", "XDG_CACHE_HOME", "XDG_STATE_HOME"):
            self.assertFalse(env[k].startswith(str(work)), k)
        self.assertFalse(Path(env["HOME"]).exists())  # removed after the run

    def test_symlinks_are_not_given_to_the_agent(self):
        tmp, _ = self.work()
        root, out = tmp / "root", tmp / "out"
        secret = tmp / "secret"
        secret.write_text("not public")
        (root / "decks").mkdir(parents=True)
        (root / "decks" / "real.py").write_text("x = 1")
        (root / "decks" / "link.py").symlink_to(secret)
        (out / "figures").mkdir(parents=True)
        (out / "figures" / "d__a.png").write_bytes(b"png")
        (out / "figures" / "d__b.png").symlink_to(secret)
        deck_list = [
            {"name": "d", "source": "decks/real.py"},
            {"name": "e", "source": "decks/link.py"},
        ]
        work, skipped = judge.prepare_workdir(
            root, out, {}, [], deck_list, tmp / "work/PROMPT.md"
        )
        self.addCleanup(shutil.rmtree, work, True)
        self.assertEqual(
            sorted(p.name for p in (work / "sources").iterdir()), ["real.py"]
        )
        self.assertEqual(
            sorted(p.name for p in (work / "figures").iterdir()), ["d__a.png"]
        )
        self.assertEqual(skipped, ["link.py", "d__b.png"])
        for p in work.rglob("*"):
            self.assertFalse(p.is_symlink(), p)

    def test_timeout_kills_children_that_hold_stdout(self):
        tmp, work = self.work()
        pidfile = tmp / "child.pid"
        # like a node/bun agent starting a worker: the child inherits stdout and outlives us
        exe = script(
            tmp / "opencode",
            "import subprocess, time\n"
            "c = subprocess.Popen(['sleep', '20'])\n"
            f"open({str(pidfile)!r}, 'w').write(str(c.pid))\n"
            "time.sleep(20)\n",
        )
        start = time.monotonic()
        out, why = judge.run_agent(exe, work, "m/x", 1)
        self.assertLess(time.monotonic() - start, 1 + 8)
        self.assertEqual((out, why), (None, "timeout after 1 s"))
        child = int(pidfile.read_text())
        time.sleep(0.2)
        try:
            os.kill(child, 0)
            alive = Path(f"/proc/{child}/stat").read_text().split()[2] != "Z"
        except (ProcessLookupError, FileNotFoundError):
            alive = False
        self.assertFalse(alive, "the agent's child survived the timeout")

    def test_spent_free_quota_stops_the_agent_at_once(self):
        # opencode 1.18.32 logs this, then retries without end and prints nothing else
        tmp, work = self.work()
        exe = script(
            tmp / "opencode",
            "import sys, time\n"
            'sys.stderr.write(\'level=ERROR message="stream error" error.error="AI_APICallError:'
            " Free usage limit reached. Free models work in OpenCode\"\\n')\n"
            "sys.stderr.flush()\ntime.sleep(30)\n",
        )
        start = time.monotonic()
        out, why = judge.run_agent(exe, work, "m/x", 20)
        self.assertLess(time.monotonic() - start, 8)
        self.assertEqual((out, why), (None, "the free usage limit of m/x is reached"))

    def test_main_never_fails(self):
        tmp, _ = self.work()
        out = tmp / "out"
        out.mkdir()
        (tmp / "p").write_text("teamcity.build.id=notanumber\n")
        argv = ["--job", "mhd", "--root", str(tmp), "--out", str(out), "--no-history"]
        argv += ["--opencode", "/nonexistent/opencode"]
        env = {"TEAMCITY_BUILD_PROPERTIES_FILE": str(tmp / "p")}
        with unittest.mock.patch.dict(os.environ, env):
            self.assertEqual(judge.main(argv), 0)  # no summary.json
            meta = json.loads((out / "meta.json").read_text())
            self.assertEqual(meta["agent_status"], "unavailable")
            (out / "summary.json").write_text(json.dumps({"job": "mhd", "decks": []}))
            self.assertEqual(judge.main(argv), 0)  # no agent executable
        meta = json.loads((out / "meta.json").read_text())
        self.assertEqual(
            (meta["agent_status"], meta["build_id"]), ("unavailable", None)
        )
        self.assertIn("could not", meta["agent_reason"])
        self.assertFalse((out / "agent.md").exists())

    def run_main(
        self,
        stdout,
        rc=0,
        prepare=None,
        decks_extra=None,
        decks_rc=0,
        props_extra="",
        version="1.18.32",
        stderr="",
        history_get=None,
    ):
        """decks_extra: {deck name: fields added to its summary.json entry}."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        out = Path(tmp.name, "out")
        out.mkdir()
        self.last_out = out
        if prepare:
            prepare(Path(tmp.name), out)
        entries = [
            {"name": n, "numbers": {}, "figures": [], "freshness": "ok"}
            | (decks_extra or {}).get(n, {})
            for n in ("convergence", "orszag_tang")
        ]
        (out / "summary.json").write_text(
            json.dumps({"job": "mhd", "decks_exit_code": decks_rc, "decks": entries})
        )
        exe = fake_opencode(tmp.name, stdout, rc, version, stderr)
        props = Path(tmp.name, "build.properties")
        props.write_text(
            "teamcity.build.id=123\nteamcity.buildType.id=BT\n" + props_extra
        )
        argv = [
            "--job",
            "mhd",
            "--root",
            tmp.name,
            "--out",
            str(out),
            "--opencode",
            exe,
        ]
        env = {"TEAMCITY_BUILD_PROPERTIES_FILE": str(props)}
        with (
            unittest.mock.patch.dict(os.environ, env),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            judge.main(
                argv + ([] if history_get else ["--no-history"]), get=history_get
            )
        agent = out / "agent.md"
        return (agent.read_text() if agent.exists() else None), json.loads(
            (out / "meta.json").read_text()
        )

    def test_main_writes_meta_json_and_the_agents_last_block(self):
        def plant(tmp, out):
            (out / "figures").mkdir()
            (out / "figures" / "convergence__x.png").symlink_to(tmp / "elsewhere")

        echo = "<<<REPORT\n### convergence\nVerdict: suspicious\ntemplate\nREPORT>>>\n"
        agent, meta = self.run_main(
            echo
            + AGENT_OK,  # the agent may echo the template first: the last block wins
            prepare=plant,
            version="opencode 1.18.32\n",
            props_extra="build.vcs.number=abc123\n",
        )
        # the agent's text alone, nothing of ours: publish fences the model's text by itself
        self.assertEqual(agent, AGENT_BLOCK)
        self.assertEqual(
            sorted(p.name for p in self.last_out.iterdir()),
            ["agent.md", "figures", "meta.json", "summary.json"],
        )
        self.assertEqual(
            meta,
            {
                "job": "mhd",
                "build_id": 123,
                "build_type_id": "BT",
                "sha": "abc123",
                "model": "opencode/deepseek-v4.1-flash-free",
                "opencode_version": "1.18.32",
                "agent_status": "ok",
                "agent_reason": None,
                "unavailable_nights": 0,
                "not_given": ["convergence__x.png"],
            },
        )

    def judge_custom(self, cfg_text, deck_names):
        """judge.main on a custom decks.yml (job 'j'); the fake agent records the deck names in
        the summary.json it was given and claims a verdict for every deck of the config."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        cfg = Path(tmp.name, "decks.yml")
        cfg.write_text(textwrap.dedent(cfg_text))
        out = Path(tmp.name, "out")
        out.mkdir()
        (out / "summary.json").write_text(
            json.dumps({"job": "j", "decks": [{"name": n} for n in deck_names]})
        )
        seen = Path(tmp.name, "seen.json")
        answer = (
            "<<<REPORT\n"
            + "".join(
                f"### {n}\nVerdict: looks right\nforged for {n}\n" for n in deck_names
            )
            + "REPORT>>>\n"
        )
        exe = script(
            Path(tmp.name, "opencode"),
            "import json, sys\n"
            "if sys.argv[1:2] == ['--version']:\n    print('1.18.32'); sys.exit(0)\n"
            "names = [d['name'] for d in json.load(open('summary.json'))['decks']]\n"
            f"open({str(seen)!r}, 'w').write(json.dumps(names))\n"
            f"sys.stdout.write({answer!r})\n",
        )
        argv = ["--job", "j", "--root", tmp.name, "--out", str(out), "--opencode", exe]
        env = {"TEAMCITY_BUILD_PROPERTIES_FILE": str(Path(tmp.name, "none"))}
        with unittest.mock.patch.dict(os.environ, env):
            judge.main(argv + ["--config", str(cfg), "--no-history"])
        seen_names = json.loads(seen.read_text()) if seen.exists() else None
        return json.loads((out / "meta.json").read_text()), seen_names

    CFG_AB = """\
        model: m
        opencode_version: 1.18.32
        page: {branch: b, dots: 5}
        jobs:
          j:
            - {name: a, source: a.py, workdir: w, figures: [], extractor: figures_only,
               theory: {quantity: q, value: 1, deck_tolerance: null}}
            - {name: b, source: b.py, workdir: w, figures: [], extractor: figures_only,
               theory: {quantity: q, value: 1, deck_tolerance: null},
               judge: false, judge_reason: "setup under review"}
        """

    def test_a_deck_with_judge_false_is_not_shown_to_the_agent(self):
        meta, seen = self.judge_custom(self.CFG_AB, ["a", "b"])
        self.assertEqual((seen, meta["agent_status"]), (["a"], "ok"))
        none_judged = self.CFG_AB.replace(
            "- {name: a, source: a.py,",
            "- {name: a, source: a.py, judge: false, judge_reason: x,",
        )
        meta, seen = self.judge_custom(none_judged, ["a", "b"])
        self.assertEqual((seen, meta["agent_status"]), (None, "not run"))

    GONE = "model opencode/deepseek-v4.1-flash-free is no longer available on OpenCode"

    def test_agent_failures_are_unavailable_with_a_plain_reason(self):
        # TODO(Task 7): replace by the exact message opencode 1.18.32 prints for a withdrawn model
        gone = "Error: ProviderModelNotFoundError: model not found"
        for kw, reason in (
            (
                {"stdout": "", "rc": 1, "stderr": "network unreachable"},
                "opencode exit code 1",
            ),
            ({"stdout": "", "rc": 1, "stderr": gone}, self.GONE),
            ({"stdout": gone}, self.GONE),
            ({"stdout": "no block"}, "no <<<REPORT block in the agent's answer"),
        ):
            with self.subTest(**kw):
                agent, meta = self.run_main(**kw)
                self.assertIsNone(agent)
                self.assertEqual(
                    (
                        meta["agent_status"],
                        meta["agent_reason"],
                        meta["unavailable_nights"],
                    ),
                    ("unavailable", reason, 1),
                )

    def test_consecutive_unavailable_nights_are_counted(self):
        statuses = {
            122: "unavailable",
            121: "unavailable",
            120: "ok",
            119: "unavailable",
        }

        def get(url):
            if "/app/rest/builds?" in url:
                return json.dumps({"build": [{"id": i} for i in statuses]}).encode()
            build = int(url.split(":id/")[0].rsplit("/", 1)[1])
            if url.endswith("/meta.json"):
                return json.dumps({"agent_status": statuses[build]}).encode()
            return json.dumps({"decks": []}).encode()

        _, meta = self.run_main("", rc=1, history_get=get)
        self.assertEqual(meta["unavailable_nights"], 3)
        _, meta = self.run_main(AGENT_OK, history_get=get)
        self.assertEqual(meta["unavailable_nights"], 0)
        statuses[122] = "ok"
        _, meta = self.run_main("", rc=1, history_get=get)
        self.assertEqual(meta["unavailable_nights"], 1)

    def test_a_judge_crash_says_when_and_keeps_a_failed_own_check(self):
        """Review Q39: deck_passed false, a failing job, and judge failing before the agent ran
        (an unreadable deck source makes prepare_workdir raise)."""
        source = decks.job(decks.load(), "mhd")[0]["source"]

        def unreadable_source(tmp, out):
            src = tmp / source
            src.parent.mkdir(parents=True)
            src.write_text("x = 1")
            src.chmod(0)

        _, meta = self.run_main(
            AGENT_OK,
            prepare=unreadable_source,
            decks_rc=1,
            decks_extra={"convergence": {"deck_passed": False}},
        )
        self.assertEqual(meta["agent_status"], "unavailable")
        self.assertEqual(meta["agent_reason"], "PermissionError before the agent ran")
        caps = {d["name"]: 6 for d in decks.job(decks.load(), "mhd")}
        part = publish.read_part(self.last_out, "mhd", caps)
        self.assertEqual(part["agent"], "")  # the agent never ran
        deck = next(d for d in part["summary"]["decks"] if d["name"] == "convergence")
        rec = history.test_record(deck, decks.job(decks.load(), "mhd")[0], part, None)
        self.assertEqual(rec["status"], "own check failed")
        with unittest.mock.patch.object(
            judge, "last_block", side_effect=RuntimeError("boom")
        ):
            _, meta = self.run_main(AGENT_OK)
        self.assertEqual(meta["agent_reason"], "RuntimeError after the agent ran")


class PageTest(unittest.TestCase):
    SHA_A, SHA_B = "a" * 40, "b" * 40
    BASE = "https://tc.example/teamcity"

    def cfg(self, dots=90, model="m1"):
        return {
            "model": model,
            "teamcity": {"url": self.BASE},
            "page": {"dots": dots},
        }

    def deck_cfg(self, names, **extra):
        return {
            "jobs": {
                "hybrid": [
                    dict(
                        {
                            "name": n,
                            "theory": {
                                "quantity": "q",
                                "value": 0.09,
                                "deck_tolerance": 0.02,
                            },
                            "key_number": "gamma",
                        },
                        **extra.get(n, {}),
                    )
                    for n in names
                ]
            }
        }

    def rec(self, status, **kw):
        return dict(
            {
                "status": status,
                "job": "hybrid",
                "deck_passed": None,
                "key": None,
                "reason": None,
            },
            **kw,
        )

    def night(self, date, tests, sha=None, manual=False, build=7, start=None):
        return {
            "date": date,
            "run": {
                "sha": sha or self.SHA_A,
                "manual": manual,
                "time": "03:15",
                "builds": {
                    "hybrid": {"id": build, "type": "Hy", "start": start or date}
                },
            },
            "tests": tests,
        }

    def part(self, names, sections=None, exit_code=0, figures=None):
        return {
            "job": "hybrid",
            "missing": False,
            "meta": {},
            "summary": {
                "decks_exit_code": exit_code,
                "decks": [
                    {"name": n, "figures": (figures or {}).get(n, [f"{n}__a.png"])}
                    for n in names
                ],
            },
            "sections": sections or {},
        }

    def render(
        self,
        hist,
        names,
        parts=None,
        recorded=True,
        notes=(),
        tonight=None,
        cfg=None,
        deck_cfg=None,
    ):
        return page.render(
            hist,
            tonight or hist[-1],
            recorded,
            (cfg or self.cfg()) | (deck_cfg or self.deck_cfg(names)),
            parts if parts is not None else [self.part(names)],
            list(notes),
            lambda job, build, name: (
                build and f"https://fig.example/{build['id']}/{job}/{name}"
            ),
        )

    def test_abnormal_tests_come_first_and_large_the_others_as_cards(self):
        names = ["a", "b", "c", "d"]
        hist = [
            self.night(
                "2026-10-28",
                {
                    "a": self.rec("looks right", key=0.0875, reason="fine"),
                    "b": self.rec("suspicious", key=0.5),
                    "c": self.rec("not judged", reason="why not"),
                    "d": self.rec("not run"),
                },
            )
        ]
        parts = [
            self.part(
                names,
                sections={n: ("looks right", f"theory of {n}") for n in names},
                figures={n: [f"{n}__1.png", f"{n}__2.png"] for n in names},
            )
        ]
        dc = self.deck_cfg(names, a={"own_check": "none"}, b={"own_check": "none"})
        out = self.render(hist, names, parts, deck_cfg=dc)
        self.assertEqual(
            re.findall(r'<section id="(\w+)" class="(\w+) (\w+)"', out),
            [
                ("b", "test", "abnormal"),
                ("a", "card", "nominal"),
                ("c", "card", "neutral"),
                ("d", "card", "notrun"),
            ],
        )
        self.assertIn('<a class="tile abnormal" href="#b">', out)
        quiet = out[out.index('<div class="quiet">') : out.index("<footer")]
        loud = out[: out.index('<div class="quiet">')]
        self.assertIn('class="big"', loud)
        self.assertNotIn('class="big"', quiet)
        self.assertIn("theory of b", loud)
        self.assertNotIn("theory of a", out)  # no model text on a test that looks right
        self.assertIn('<div class="analysis">theory of c</div>', quiet)
        self.assertEqual(quiet.count("<img"), 3)  # a card shows its first figure only
        for want in (
            "gamma = 0.5",
            "theory q: 0.09 ± 0.02",
            "the test script checks nothing itself",
        ):
            self.assertIn(want, loud)
        for want in (
            "gamma = 0.0875",
            "fine",
            "why not",
            'class="dots"',
            "the test script checks nothing itself",
        ):
            self.assertIn(want, quiet)
        self.assertNotIn("theory q", quiet)

    def test_a_deck_not_judged_tonight_is_not_carried_over(self):
        hist = [
            self.night("2026-10-27", {"a": self.rec("suspicious")}),
            self.night("2026-10-28", {"a": self.rec("not judged")}),
        ]
        out = self.render(
            hist, ["a"], deck_cfg=self.deck_cfg(["a"], a={"judge": False})
        )
        self.assertIn('<a class="tile neutral" href="#a">', out)
        self.assertNotIn("last judged", out)

    def test_the_page_escapes_what_it_is_given(self):
        """page.render escapes on its own, whatever clean_part or history.load let through."""
        evil = '"><img src="//evil.example/x.png" onload="x"><a href="//evil.example">'
        hist = [
            self.night("2026-10-27", {"a": self.rec(evil, reason=evil)}),
            self.night(
                "2026-10-28",
                {"a": self.rec("suspicious", reason=evil), "b": self.rec(evil)},
                manual=True,
            ),
        ]
        hist[1]["run"]["time"] = evil
        figures = [f"a__{evil}.png", f"a__{evil}.mp4"]
        parts = [
            self.part(
                ["a", "b"],
                sections={"a": ("suspicious", evil)},
                figures={"a": figures, "b": figures},
            )
        ]
        out = self.render(
            hist, ["a", "b"], parts, notes=[evil], cfg=self.cfg(model=evil)
        )
        links = re.findall(r'(?:src|href)="([^"]*)"', out)
        self.assertTrue(links)
        for url in links:
            self.assertTrue(
                url.startswith(
                    (
                        self.BASE + "/",
                        "https://fig.example/",
                        "https://github.com/",
                        "#",
                    )
                ),
                url,
            )
        for raw in ('<img src="//evil', 'onload="', '<a href="//evil'):
            self.assertNotIn(raw, out)
        self.assertEqual(out.count("<script>"), 1)

    def test_model_text_cap(self):
        hist = [self.night("2026-10-28", {"a": self.rec("suspicious")})]
        parts = [self.part(["a"], sections={"a": ("suspicious", "y" * 10000)})]
        out = self.render(hist, ["a"], parts)
        run = max(len(m) for m in re.findall(r"y+(?: \[…\])?", out))
        self.assertLessEqual(run, 4000 + len(" […]"))
        self.assertIn("[…]", out)

    def test_dots(self):
        hist = [
            self.night("2026-10-24", {"a": self.rec("looks right")}),  # beyond `dots`
            self.night("2026-10-25", {}),  # no record that night: no dot
            self.night("2026-10-26", {"a": self.rec("looks right")}),
            self.night("2026-10-27", {"a": self.rec("looks right", expired=True)}),
            self.night(
                "2026-10-28", {"a": self.rec("suspicious")}, build=9, start="2026-10-27"
            ),
        ]
        out = self.render(hist, ["a"], cfg=self.cfg(dots=4))
        self.assertEqual(out.count('class="dot '), 3)
        self.assertNotIn('title="2026-10-24', out)
        self.assertIn(
            f'<a class="dot s-ok" title="2026-10-26 · looks right" '
            f'href="{self.BASE}/viewLog.html?buildId=7&amp;tab=artifacts">',
            out,
        )
        self.assertIn('<span class="dot ring s-ok" title="2026-10-27', out)  # expired
        self.assertRegex(
            out,
            r'class="dot s-suspicious" title="2026-10-28[^"]*reused from 2026-10-27',
        )
        self.assertLess(out.index("2026-10-26 ·"), out.index("2026-10-27 ·"))

    def test_new_tonight_links_the_compare_view_between_two_known_commits(self):
        A, B = self.SHA_A, self.SHA_B
        for prev, now, link, pill in (
            (
                A,
                B,
                f"https://github.com/PHAREHUB/PHARE/compare/{A}...{B}",
                "new tonight",
            ),
            (A, A, None, "new tonight, no new commit"),
            (None, B, None, "new tonight"),
        ):
            hist = [
                self.night("2026-10-27", {"a": self.rec("looks right")}),
                self.night("2026-10-28", {"a": self.rec("suspicious")}, sha=now),
            ]
            hist[0]["run"]["sha"] = prev
            out = self.render(hist, ["a"])
            self.assertIn(pill, out)
            if link:
                self.assertIn(link, out)
            else:
                self.assertNotIn("compare/", out)

    def test_carry_over_says_last_judged(self):
        hist = [
            self.night("2026-10-26", {"a": self.rec("suspicious")}),
            self.night("2026-10-27", {"a": self.rec("suspicious")}),
            self.night("2026-10-28", {"a": self.rec("cannot tell")}),
        ]
        out = self.render(hist, ["a"])
        self.assertIn("since 2026-10-26", out)
        self.assertIn("last judged 2026-10-27, 1 night ago", out)
        hist[2]["tests"]["a"] = self.rec("suspicious")
        self.assertNotIn("last judged", self.render(hist, ["a"]))

    def test_movie_is_a_video(self):
        hist = [self.night("2026-10-28", {"shock": self.rec("suspicious")})]
        parts = [
            self.part(
                ["shock"],
                figures={"shock": ["shock__x.mp4", "shock__y.png", "shock__z.png"]},
            )
        ]
        out = self.render(hist, ["shock"], parts)
        self.assertIn('<video controls preload="none"', out)
        self.assertIn("https://fig.example/7/hybrid/shock__x.mp4", out)
        self.assertIn('class="zoom" data-group="shock"', out)
        self.assertEqual(out.count('loading="lazy"'), 2)
        hist[0]["tests"]["shock"] = self.rec("looks right")  # a card shows one figure
        out = self.render(hist, ["shock"], parts)
        self.assertEqual(out.count("<video"), 1)
        self.assertNotIn('class="zoom" data-group', out)

    def test_unknown_build_draws_no_broken_figure(self):
        hist = [self.night("2026-10-28", {"a": self.rec("looks right")})]
        hist[0]["run"]["builds"] = {}
        out = self.render(hist, ["a"])
        self.assertNotIn('class="zoom"', out)
        self.assertNotIn("https://fig.example", out)
        self.assertIn("figures unavailable", out)

    def test_failed_tests_step_banner(self):
        hist = [
            self.night(
                "2026-10-28",
                {"a": self.rec("own check unknown", own_check_unknown=True)},
            )
        ]
        said = {"a": ("looks right", "the model's reading")}
        out = self.render(hist, ["a"], [self.part(["a"], said, exit_code=3)])
        self.assertIn("hybrid tests step failed (exit 3): own checks unknown", out)
        self.assertIn('<div class="analysis">the model&#x27;s reading</div>', out)
        self.assertEqual(out.count(">own check unknown</span>"), 1)
        ok = self.render(hist, ["a"], [self.part(["a"], exit_code=None)])
        self.assertNotIn("tests step failed", ok)

    def test_counts_follow_the_table(self):
        statuses = [
            "own check failed",
            "suspicious",
            "looks right",
            "cannot tell",
            "not judged",
            "no analysis",
            "own check unknown",
            "no data",
            "not run",
        ]
        names = [f"t{i}" for i in range(9)]
        hist = [
            self.night("2026-10-28", {n: self.rec(s) for n, s in zip(names, statuses)})
        ]
        out = self.render(hist, names)
        self.assertIn("2 to look at · 1 fine · 5 can't tell · 1 not run", out)

    def test_unrecorded_manual_run_draws_tonight(self):
        sched = self.night("2026-10-28", {"a": self.rec("looks right")})
        manual = self.night("2026-10-28", {"a": self.rec("suspicious")}, manual=True)
        out = self.render([sched], ["a"], recorded=False, tonight=manual)
        self.assertIn('<a class="tile abnormal" href="#a">', out)
        self.assertIn("not recorded in the history", out)
        self.assertIn("manual run 03:15", out)
        self.assertEqual(out.count('class="dot '), 1)


SHA_ABC, SHA_DEF = "abc".ljust(40, "0"), "def".ljust(40, "0")
META = {
    "build_id": 5,
    "build_type_id": "BT",
    "model": "m",
    "opencode_version": "1",
    "agent_status": "ok",
}


# a meta.json written by hostile code in the functional job's pod (review Q44)
EVIL_META = {
    "build_id": 5,
    "sha": "../../issues/1/comments",
    "opencode_version": '1.18.32<img src="//evil.example/v.png">',
    "agent_status": "ok",
    "agent_reason": None,
    "not_given": ['<img src="//evil.example/n.png">'],
}


TC_BASE = "https://hephaistos.lpp.polytechnique.fr/teamcity"
TONIGHT = "20261028T010012+0000"  # 02:00 in Paris, same day
YESTERDAY = "20261027T010012+0000"
SECRET = "ghp_SECRET123"


def tc_build(i, type_, start=TONIGHT, default=True, trigger="schedule", sha=SHA_ABC):
    return {
        "id": i,
        "buildTypeId": type_,
        "startDate": start,
        "defaultBranch": default,
        "triggered": {"type": trigger},
        "revisions": {"revision": [{"version": sha}]},
    }


def tc_nightly(*upstream):
    return list(upstream) or [
        tc_build(11, HYB_TYPE),
        tc_build(12, MHD_TYPE),
    ]


def fake_teamcity(
    upstream,
    own_start=TONIGHT,
    head_trigger="schedule",
    heads=(9,),
    head_default=True,
    own_default=True,
    head_broken=False,
):
    """A get for the TeamCity reads of publish: the page writer is build 5, the head 9."""
    urls = []

    def get(url, *a, **k):
        urls.append(url)
        if "/artifacts/children/" in url:
            return b"[]"
        if "snapshotDependency:(from:" in url:
            return json.dumps(
                {"build": [{"id": h, "buildTypeId": "T_head"} for h in heads]}
            ).encode()
        if "snapshotDependency:(to:" in url:
            return json.dumps({"build": upstream}).encode()
        i = int(re.search(r"builds/id:(\d+)\?", url).group(1))
        if i == 9 and head_broken:
            return b"{}"
        start, trigger = (own_start, "schedule") if i == 5 else (TONIGHT, head_trigger)
        default = own_default if i == 5 else head_default
        return json.dumps(
            tc_build(i, "T_x", start, trigger=trigger, default=default)
        ).encode()

    get.urls = urls
    return get


def down(url, *a, **k):
    raise OSError("down")


class PublishTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = Path(tmp.name)
        self.site = self.d / "site"
        self.site.mkdir()
        env = unittest.mock.patch.dict(
            os.environ, {"TEAMCITY_BUILD_PROPERTIES_FILE": str(self.d / "none")}
        )
        env.start()
        self.addCleanup(env.stop)
        os.environ.pop("NIGHTLY_PAGE_TOKEN", None)

    def write_job(self, job, decks_in, agent="", meta=None, status=None):
        j = self.d / "in" / job
        j.mkdir(parents=True, exist_ok=True)
        (j / "meta.json").write_text(
            json.dumps(
                dict(
                    META,
                    build_type_id=HYB_TYPE,
                    build_id=11,
                    opencode_version="1.18.32",
                )
                | (meta or {})
            )
        )
        decks_out = [
            {
                "name": n,
                "figures": [f"{n}__a.png"],
                "numbers": {"gamma": 0.09},
                "status": status,
            }
            for n in decks_in
        ]
        (j / "summary.json").write_text(
            json.dumps({"decks_exit_code": 0, "decks": decks_out})
        )
        (j / "agent.md").write_text(agent)

    def write_both(self, hybrid_meta=None, mhd_meta=None):
        self.write_job(
            "hybrid",
            ("ion_ion_beam1d", "td1d"),
            "### ion_ion_beam1d\nVerdict: looks right\nall good\n",
            hybrid_meta,
        )
        self.write_job(
            "mhd",
            ("convergence", "orszag_tang"),
            "### convergence\nVerdict: suspicious\norder 1 not 3\n"
            "### orszag_tang\nVerdict: looks right\nfine\n",
            mhd_meta,
        )

    def edit_summary(self, job, change):
        p = self.d / "in" / job / "summary.json"
        s = json.loads(p.read_text())
        change(s)
        p.write_text(json.dumps(s))

    def main(self, get, extra=(), runner=None, jobs="hybrid,mhd", own=True):
        argv = ["--in", str(self.d / "in"), "--jobs", jobs, "--site", str(self.site)]
        argv += ["--own-build-id", "5"] if own else []
        out = io.StringIO()
        kwargs = {"run": runner} if runner else {}
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            rc = publish.main(argv + list(extra), get=get, **kwargs)
        self.printed = out.getvalue()
        return rc

    def index(self):
        return (self.site / "index.html").read_text()

    def nights(self, name="history.json"):
        return json.loads((self.site / name).read_text())["nights"]

    def test_writes_index_and_history(self):
        self.write_both()
        self.write_job(
            "mhd",
            ("convergence", "orszag_tang"),
            "### convergence\nVerdict: suspicious\norder 1 not 3\n"
            "### orszag_tang\nVerdict: own check failed\n",  # only the scripts can say so
        )
        get = fake_teamcity(tc_nightly())
        self.assertEqual(self.main(get), 0)
        (night,) = self.nights()
        self.assertEqual(night["date"], "2026-10-28")
        self.assertEqual(night["run"]["sha"], SHA_ABC)
        self.assertIs(night["run"]["manual"], False)
        self.assertEqual(
            night["run"]["builds"],
            {
                "hybrid": {"id": 11, "type": HYB_TYPE, "start": "2026-10-28"},
                "mhd": {"id": 12, "type": MHD_TYPE, "start": "2026-10-28"},
            },
        )
        status = {n: r["status"] for n, r in night["tests"].items()}
        self.assertEqual(
            status,
            {
                "ion_ion_beam1d": "looks right",
                "td1d": "not judged",
                "convergence": "suspicious",
                "orszag_tang": "cannot tell",
            },
        )
        self.assertEqual(night["tests"]["ion_ion_beam1d"]["key"], 0.09)
        html = self.index()
        self.assertIn("order 1 not 3", html)  # the abnormal test shows the model's text
        self.assertIn(
            f"{TC_BASE}/guestAuth/repository/download/{HYB_TYPE}/11:id/"
            "physics-out/hybrid/figures/ion_ion_beam1d__a.png",
            html,
        )
        self.assertFalse((self.site / "history.pending.json").exists())
        self.assertEqual(self.printed, "")

    def test_reused_upstream_build_keeps_tonights_date(self):
        self.write_both()
        old = tc_build(11, HYB_TYPE, YESTERDAY, trigger="user")
        get = fake_teamcity(tc_nightly(old, tc_build(12, MHD_TYPE)))
        self.main(get)
        (night,) = self.nights()
        self.assertEqual(night["date"], "2026-10-28")
        self.assertIs(night["run"]["manual"], False)
        self.assertEqual(night["run"]["builds"]["hybrid"]["start"], "2026-10-27")

    def test_one_night_per_date_and_a_manual_one_never_replaces_a_scheduled_one(self):
        self.write_both()
        self.main(fake_teamcity(tc_nightly(), head_trigger="user"))
        (night,) = self.nights()
        self.assertEqual(
            (night["run"]["manual"], night["run"]["time"]), (True, "02:00")
        )
        self.assertIn("manual run 02:00", self.index())
        self.main(fake_teamcity(tc_nightly()))  # the scheduled night replaces it
        self.main(fake_teamcity(tc_nightly(), head_trigger="user"))
        # an unreadable head: unknown trigger, so manual
        self.main(fake_teamcity(tc_nightly(), head_broken=True))
        (night,) = self.nights()
        self.assertIs(night["run"]["manual"], False)
        self.assertIn("not recorded in the history", self.index())

    def test_teamcity_is_asked_once_per_build(self):
        self.write_both()
        get = fake_teamcity(tc_nightly(), heads=())  # no head: the writer is the start
        self.main(get)
        own = [u for u in get.urls if "builds/id:5?" in u]
        self.assertEqual(len(own), 1)
        asked = [u for u in get.urls if "/artifacts/children/" in u]
        self.assertTrue(asked)
        self.assertEqual(len(asked), len(set(asked)))

    def test_the_commit_comes_from_teamcity(self):
        self.write_both(hybrid_meta={"sha": SHA_DEF})
        self.main(fake_teamcity(tc_nightly()))
        self.assertEqual(self.nights()[0]["run"]["sha"], SHA_ABC)
        self.assertIn("differs from TeamCity", self.index())
        disagree = tc_nightly(
            tc_build(11, HYB_TYPE), tc_build(12, MHD_TYPE, sha=SHA_DEF)
        )
        self.main(fake_teamcity(disagree))
        self.assertIsNone(self.nights()[0]["run"]["sha"])

    def test_the_builds_are_found_from_the_writer_and_its_pipeline_head(self):
        self.write_both()
        props = self.d / "props"
        props.write_text("teamcity.build.id=5\n")
        # undone by setUp's patch.dict
        os.environ["TEAMCITY_BUILD_PROPERTIES_FILE"] = str(props)
        get = fake_teamcity(tc_nightly(), heads=(8, 9))
        self.main(get, own=False)
        urls = "\n".join(get.urls)
        self.assertIn("snapshotDependency:(from:(id:5))", urls)
        self.assertIn("snapshotDependency:(to:(id:9))", urls)  # the newest head
        self.assertEqual(len(self.nights()[0]["run"]["builds"]), 2)
        self.assertIn("other pipeline heads found: 8", self.index())
        get = fake_teamcity(tc_nightly(), heads=(8, 9))
        self.main(get, ["--head-build-id", "77"])
        urls = "\n".join(get.urls)
        self.assertNotIn("snapshotDependency:(from:", urls)
        self.assertIn("snapshotDependency:(to:(id:77))", urls)

    def test_without_teamcity_the_night_is_today_and_manual_with_no_builds(self):
        self.write_both()
        self.assertEqual(self.main(down), 0)
        (night,) = self.nights()
        self.assertRegex(night["date"], r"\d{4}-\d{2}-\d{2}")
        self.assertEqual(night["run"]["builds"], {})
        self.assertIsNone(night["run"]["sha"])
        self.assertIs(night["run"]["manual"], True)  # an unknown trigger is manual

    def test_no_tzdata_means_utc_and_a_note(self):
        self.write_both()
        with unittest.mock.patch.object(publish, "ZoneInfo", side_effect=KeyError("x")):
            self.main(fake_teamcity(tc_nightly()))
        self.assertIn("UTC", self.index())
        self.assertEqual(self.nights()[0]["date"], "2026-10-28")

    def test_never_a_gate(self):
        """An internal error: exit 0, one printed line, yesterday's page left as it was."""
        self.write_both()
        (self.site / "index.html").write_text("yesterday")
        with unittest.mock.patch.object(
            page, "render", side_effect=RuntimeError("boom")
        ):
            self.assertEqual(self.main(fake_teamcity(tc_nightly())), 0)
        self.assertEqual(self.index(), "yesterday")
        self.assertEqual(len(self.printed.strip().splitlines()), 1)
        self.assertIn("boom", self.printed)
        self.assertFalse((self.site / "history.json").exists())
        self.assertEqual(self.main(down, ["--nope"]), 0)  # bad arguments

    def test_unreadable_history_is_not_overwritten(self):
        self.write_both()
        (self.site / "history.json").write_text("{not json")
        get = fake_teamcity(tc_nightly())
        self.assertEqual(self.main(get), 0)
        self.assertEqual((self.site / "history.json").read_text(), "{not json")
        (night,) = self.nights("history.pending.json")
        self.assertEqual(night["date"], "2026-10-28")
        self.assertIn("history.json has problems: 1 nights kept aside", self.index())
        # fixed by hand: the next run merges what was kept aside, and deletes it
        history.dump([], self.site / "history.json")
        self.main(get)
        self.assertEqual([n["date"] for n in self.nights()], ["2026-10-28"])
        self.assertFalse((self.site / "history.pending.json").exists())

    def test_notes_from_the_parts(self):
        self.write_both(
            hybrid_meta={
                "agent_status": "unavailable",
                "agent_reason": "timeout",
                "unavailable_nights": 3,
                "opencode_version": "1.19.0",
                "not_given": ["a.png"],
            }
        )

        def change(s):
            s["decks"][0]["freshness"] = "unknown"
            s["decks"][0]["numbers"] = {"bad key": 1}

        self.edit_summary("hybrid", change)
        self.main(fake_teamcity(tc_nightly()))
        html = self.index()
        for note in (
            "agent unavailable 3 nights in a row (hybrid)",
            "freshness unknown: figures may be from an earlier run (hybrid)",
            "OpenCode 1.19.0 differs from the pinned 1.18.32",
            "not given to the agent (symlink or not a regular file): a.png",
            "number ion_ion_beam1d.bad key",
        ):
            self.assertIn(note, html)

    def test_hostile_input_renders_no_foreign_link(self):
        """The pod's meta.json, summary.json and agent.md, and a hand-edited history.json,
        are untrusted: every src and href on the page goes to TeamCity, GitHub or an anchor."""
        evil = "https://evil.example/x.png"
        prev = history.night(
            "2026-10-27",
            {
                "sha": SHA_ABC,
                "manual": False,
                "time": "01:00",
                "builds": {"mhd": {"id": 12, "type": MHD_TYPE, "start": "2026-10-27"}},
            },
            {
                "convergence": {
                    "status": '<a href="//evil.example">x</a>',
                    "job": "mhd",
                    "reason": f'"><img src="{evil}" onload="x">',
                },
                "orszag_tang": {"status": 'x y" onload="alert(1)', "job": "mhd"},
            },
        )
        # nights with a bad build id, type or start are dropped whole, never linked
        bad_builds = [
            {"id": "13ZZbad", "type": MHD_TYPE, "start": "2026-10-20"},
            {"id": 13, "type": "ZZbad/../x", "start": "2026-10-20"},
            {"id": 13, "type": MHD_TYPE, "start": "0000ZZbad"},
        ]
        nights = []
        for i, b in enumerate(bad_builds + [None] * 3):
            n = dict(
                prev,
                date=f"2026-10-2{i}",
                tests={"convergence": {"status": "looks right", "job": "mhd"}},
            )
            n["run"] = dict(
                prev["run"], builds={"mhd": b or prev["run"]["builds"]["mhd"]}
            )
            nights.append(n)
        history.dump(nights + [prev], self.site / "history.json")
        self.write_both()
        self.write_job(
            "mhd",
            ("convergence", "orszag_tang"),
            f'### convergence\nVerdict: suspicious\n<img src={evil}>\n<a href="{evil}">x</a>',
            EVIL_META
            | {"sha": SHA_DEF, "build_id": f'5"><img src="{evil}'}
            | {"agent_reason": f'<img src="{evil}">'},
        )

        def change(s):
            s["decks"][0]["figures"] = [
                f'convergence__x.png"><img src="{evil}',
                "convergence__ok.svg",
                "convergence__ok.png",
            ]
            s["decks"][1]["status"] = f'extract error "><img src="{evil}">'
            s["decks"].append({"name": f'<a href="{evil}">click</a>', "figures": []})

        self.edit_summary("mhd", change)
        self.main(fake_teamcity(tc_nightly()))
        html = self.index()
        links = re.findall(r'(?:src|href)="([^"]*)"', html)
        self.assertTrue(links)
        for url in links:
            self.assertTrue(
                url.startswith(
                    (TC_BASE + "/", "https://github.com/PHAREHUB/PHARE/", "#")
                ),
                url,
            )
        for raw in (
            "<img src=https",
            '<img src="https://evil',
            'onload="',
            '<a href="//evil',
        ):
            self.assertNotIn(raw, html)
        self.assertEqual(html.count("<script>"), 1)
        self.assertNotIn("ZZbad", html)
        self.assertIn("history.json: dropped 3 bad night(s)", html)
        self.assertIn("convergence__ok.png", "".join(links))
        self.assertNotIn("convergence__ok.svg", "".join(links))

    def test_clean_part_meta(self):
        meta = EVIL_META | {
            "build_id": '5"><img src="//evil.example/b.png">',
            "agent_status": 'ok"><img src="//evil.example/s.png">',
            "unavailable_nights": '9"><img src="//evil.example/u.png">',
            "sha": "A" * 40,  # not lower-case hex
            "agent_reason": "r" * 500,
            "not_given": ["n" * 500] + [f"f{i}" for i in range(30)],
        }
        out = publish.clean_part({"job": "mhd", "meta": meta, "summary": {}}, {})
        clean = out["meta"]
        self.assertEqual(
            [
                clean[k]
                for k in ("build_id", "sha", "agent_status", "unavailable_nights")
            ],
            [None, None, None, 0],
        )
        self.assertEqual(clean["opencode_version"], "unknown")
        self.assertEqual(len(clean["agent_reason"]), 200)
        self.assertEqual(
            (len(clean["not_given"]), len(clean["not_given"][0])), (20, 200)
        )
        self.assertTrue([d for d in out["dropped"] if "not_given" in d])
        good = META | {"sha": SHA_ABC, "opencode_version": "1.19.0-beta.2"}
        clean = publish.clean_part({"job": "mhd", "meta": good, "summary": {}}, {})[
            "meta"
        ]
        self.assertEqual(
            [clean[k] for k in ("build_id", "sha", "opencode_version", "agent_status")],
            [5, SHA_ABC, "1.19.0-beta.2", "ok"],
        )

    def test_clean_part_decks(self):
        figures = [
            "ion_ion_beam1d__a.png",
            "ion_ion_beam1d__b.mp4",
            "ion_ion_beam1d__a.svg",
            "td1d__x.png",  # another deck's prefix
            'ion_ion_beam1d__c.png"><img',
            5,
        ] + [f"ion_ion_beam1d__{i}.png" for i in range(3)]
        numbers = {
            "gamma": 0.1,
            "ok": True,
            "n_int": 3,
            "bad key!": 1,
            "n": float("nan"),
            "inf": float("inf"),
            "huge": 10**400,
            "s": "x",
            "k" * 41: 1,
        }
        summary = {
            "decks_exit_code": "1",
            "decks": [
                {
                    "name": "ion_ion_beam1d",
                    "figures": figures,
                    "numbers": numbers,
                    "status": "s" * 500,
                    "freshness": "f" * 500,
                    "deck_passed": "yes",
                },
                {"name": "td1d", "numbers": [1], "status": 5, "deck_passed": False},
                {"name": "td1d"},  # a second time
                {"name": "other"},  # not a deck of this job
                "x",
            ],
        }
        out = publish.clean_part(
            {"job": "hybrid", "meta": META, "summary": summary},
            {"ion_ion_beam1d": 3, "td1d": 6},
        )
        ib, td = out["summary"]["decks"]
        self.assertEqual(
            ib["figures"],
            ["ion_ion_beam1d__a.png", "ion_ion_beam1d__b.mp4", "ion_ion_beam1d__0.png"],
        )
        self.assertEqual(ib["numbers"], {"gamma": 0.1, "ok": True, "n_int": 3})
        self.assertEqual(
            (len(ib["status"]), len(ib["freshness"]), ib["deck_passed"]),
            (200, 200, None),
        )
        self.assertEqual(
            (td["figures"], td["numbers"], td["status"], td["deck_passed"]),
            ([], {}, None, False),
        )
        self.assertIsNone(out["summary"]["decks_exit_code"])
        dropped = " | ".join(out["dropped"])
        for bad in (
            "figure ion_ion_beam1d__a.svg",
            "figure td1d__x.png",
            "figure ion_ion_beam1d__1.png over max_figures",
            "number ion_ion_beam1d.bad key!",
            "numbers of td1d",
            "status of td1d",
            "deck other",
        ):
            self.assertIn(bad, dropped)

    def test_a_part_without_readable_meta_and_summary_is_missing(self):
        with tempfile.TemporaryDirectory() as t:
            j = Path(t)
            for broken in ("meta.json", "summary.json"):
                for text in (None, "{not json", "[1]"):
                    (j / "meta.json").write_text(json.dumps({"build_id": 5}))
                    (j / "summary.json").write_text(json.dumps({"decks": []}))
                    self.assertFalse(publish.read_part(j, "mhd", {}).get("missing"))
                    if text is None:
                        (j / broken).unlink()
                    else:
                        (j / broken).write_text(text)
                    self.assertEqual(
                        publish.read_part(j, "mhd", {}), {"job": "mhd", "missing": True}
                    )

    def test_sections_keep_only_known_decks_and_strict_verdicts(self):
        names = {"td1d", "orszag_tang"}
        pre, secs = publish.sections(
            "intro\n### td1d  \nVerdict: suspicious  \nx\n### other\nVerdict: looks right\n"
            "### orszag_tang\nsomething else\n### td1d\nVerdict: looks right",
            names,
        )
        self.assertEqual(pre, "intro")
        self.assertEqual(
            secs["td1d"], ("suspicious", "x\n### other\nVerdict: looks right")
        )
        self.assertEqual(secs["orszag_tang"][0], "cannot tell")
        self.assertIn("something else", secs["orszag_tang"][1])
        # only the scripts can say "own check failed"
        _, secs = publish.sections("### td1d\nVerdict: own check failed\nx", names)
        self.assertEqual(secs["td1d"][0], "cannot tell")
        for sep in ("\u2028", "\x85", "\x0b", "\x1c"):  # str.splitlines splits on these
            _, secs = publish.sections(
                f"### td1d\nVerdict: suspicious\nx {sep}### orszag_tang{sep}Verdict: looks right",
                names,
            )
            self.assertEqual(list(secs), ["td1d"], repr(sep))

    # push

    def push_setup(self):
        """A bare repo as origin and the site as a clone of it."""
        bare = self.d / "origin.git"
        subprocess.run(["git", "init", "-q", "--bare", str(bare)], check=True)
        subprocess.run(["git", "init", "-q", "-b", "main", str(self.site)], check=True)
        subprocess.run(
            ["git", "-C", str(self.site), "remote", "add", "origin", str(bare)],
            check=True,
        )
        return bare

    @staticmethod
    def b64(token):
        return base64.b64encode(f"x-access-token:{token}".encode()).decode()

    def recorder(self):
        calls, envs = [], []

        def run(argv, **kw):
            calls.append(list(argv))
            envs.append(kw.get("env") or {})
            return subprocess.run(argv, check=False, **kw)

        run.calls, run.envs = calls, envs
        return run

    def test_push_never_exposes_the_token(self):
        bare = self.push_setup()
        self.write_both()
        run = self.recorder()
        with unittest.mock.patch.dict(os.environ, {"NIGHTLY_PAGE_TOKEN": SECRET}):
            self.assertEqual(
                self.main(fake_teamcity(tc_nightly()), ["--push"], runner=run), 0
            )
        shown = subprocess.run(
            ["git", "-C", str(bare), "show", "gh-pages:history.json"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        self.assertEqual(json.loads(shown)["nights"][0]["date"], "2026-10-28")
        log = subprocess.run(
            ["git", "-C", str(bare), "log", "--format=%s|%an|%ae", "gh-pages"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        self.assertEqual(
            log.strip(),
            "nightly 2026-10-28|phare-nightly|phare-nightly@nightly.invalid",
        )
        # the credentials travel in git's environment, scoped to GitHub
        env = run.envs[0]
        self.assertEqual(
            env["GIT_CONFIG_KEY_0"], "http.https://github.com/.extraheader"
        )
        self.assertNotIn("NIGHTLY_PAGE_TOKEN", env)
        seen = " ".join(" ".join(c) for c in run.calls)
        for secret in (SECRET, self.b64(SECRET)):
            self.assertNotIn(secret, seen)
            self.assertNotIn(secret, self.printed)
            self.assertNotIn(secret, (self.site / ".git/config").read_text())
            self.assertNotIn(secret, (bare / "config").read_text())
        self.assertEqual(
            [c[3] for c in run.calls if c[0] == "git"],
            ["add", "diff", "commit", "push"],
        )
        # a second run with nothing new commits nothing but still pushes
        run.calls.clear()
        with unittest.mock.patch.dict(os.environ, {"NIGHTLY_PAGE_TOKEN": SECRET}):
            self.main(fake_teamcity(tc_nightly()), ["--push"], runner=run)
        self.assertNotIn("commit", [c[3] for c in run.calls])

    def fake_run(self, returncode=0, stdout="", stderr=""):
        calls = []

        def run(argv, **kw):
            calls.append((list(argv), kw.get("env") or {}))
            rc = 1 if argv[3:5] == ["diff", "--cached"] else returncode
            return subprocess.CompletedProcess(argv, rc, stdout, stderr)

        run.calls = calls
        return run

    def test_a_failed_push_is_one_line_and_exit_zero(self):
        self.write_both()
        leak = f"denied {SECRET} and {self.b64(SECRET)}"
        run = self.fake_run(returncode=1, stdout=leak, stderr=leak)
        with unittest.mock.patch.dict(os.environ, {"NIGHTLY_PAGE_TOKEN": SECRET}):
            self.assertEqual(
                self.main(fake_teamcity(tc_nightly()), ["--push"], runner=run), 0
            )
        lines = self.printed.strip().splitlines()
        self.assertEqual(len(lines), 1)
        self.assertIn("denied ***", lines[0])
        for secret in (SECRET, self.b64(SECRET)):
            self.assertNotIn(secret, self.printed)
        self.assertIn("PHARE nightly physics", self.index())

    def test_push_refusals(self):
        """Only a night whose every known build is on the default branch reaches the page."""
        self.write_both()
        off = {"default": False}
        cases = {
            "a job build off the default branch": (
                fake_teamcity([tc_build(11, HYB_TYPE, **off), tc_build(12, MHD_TYPE)]),
                [],
                SECRET,
            ),
            "no build of a decks.yml job": (
                fake_teamcity([tc_build(13, "Other_virtual_build_V__1")]),
                [],
                SECRET,
            ),
            # the job builds are master builds, reused by a branch-run pipeline
            "the head off the default branch": (
                fake_teamcity(tc_nightly(), head_default=False),
                [],
                SECRET,
            ),
            "the writer off the default branch": (
                fake_teamcity(tc_nightly(), own_default=False),
                [],
                SECRET,
            ),
            "TeamCity down": (down, [], SECRET),
            "no token": (fake_teamcity(tc_nightly()), [], None),
            "local figures": (
                fake_teamcity(tc_nightly()),
                ["--figure-base", "../physics-out"],
                SECRET,
            ),
        }
        for why, (get, extra, token) in cases.items():
            with self.subTest(why):
                run = self.fake_run()
                env = {"NIGHTLY_PAGE_TOKEN": token} if token else {}
                with unittest.mock.patch.dict(os.environ, env):
                    self.assertEqual(self.main(get, ["--push", *extra], runner=run), 0)
                self.assertEqual(run.calls, [])
                self.assertEqual(len(self.printed.strip().splitlines()), 1)
                self.assertIn("push refused", self.printed)
        # the local preview links its figures without any TeamCity build
        self.main(fake_teamcity(tc_nightly()), ["--figure-base", "../in"], own=False)
        self.assertIn('src="../in/hybrid/figures/ion_ion_beam1d__a.png"', self.index())

    def test_an_unreadable_pending_file_is_set_aside_and_pushed(self):
        self.write_both()
        (self.site / "history.pending.json").write_text("{not json")
        run = self.fake_run()
        with unittest.mock.patch.dict(os.environ, {"NIGHTLY_PAGE_TOKEN": SECRET}):
            self.main(fake_teamcity(tc_nightly()), ["--push"], runner=run)
        self.assertFalse((self.site / "history.pending.json").exists())
        broken = "history.pending.broken-2026-10-28.json"
        self.assertEqual((self.site / broken).read_text(), "{not json")
        self.assertEqual(len(self.nights()), 1)  # history.json is still written
        self.assertIn("set aside", self.index())
        files = run.calls[0][0][6:]  # git -C site add -A -- <files>
        self.assertEqual(
            sorted(files),
            sorted(["index.html", "history.json", "history.pending.json", broken]),
        )


if __name__ == "__main__":
    unittest.main()
