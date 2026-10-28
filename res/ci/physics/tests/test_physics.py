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
import urllib.error
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import decks
import extract
import history
import metrics
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


if __name__ == "__main__":
    unittest.main()
