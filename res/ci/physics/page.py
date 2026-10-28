"""The nightly physics page: one self-contained static HTML string.

Pure module: no network, no import of publish. Every dynamic string goes through e();
the only <script> is the page's own lightbox.
"""

import html

import history as hist

COMMIT = "https://github.com/PHAREHUB/PHARE/commit/"
COMPARE = "https://github.com/PHAREHUB/PHARE/compare/"
MODEL_TEXT_CAP = 4000
SLUGS = {
    "own check failed": "failed",
    "suspicious": "suspicious",
    "looks right": "ok",
    "cannot tell": "cannot-tell",
    "not judged": "not-judged",
    "no analysis": "no-analysis",
    "own check unknown": "own-check-unknown",
    "no data": "no-data",
    "not run": "not-run",
}
ICONS = {"abnormal": "⚠", "nominal": "✓", "neutral": "?", "notrun": "-"}

CSS = """
:root{--bg:#fbfbfa;--fg:#1d232a;--muted:#6b7280;--card:#fff;--line:#e5e7eb;
--ok-bg:#eef8f1;--ok-line:#a7d7b5;--ok-fg:#276b3b;--warn-bg:#fff3dc;--warn-line:#d97706;
--warn-fg:#92400e;--grey-bg:#f3f4f6;--grey-fg:#9aa1a9;--note:#fffaf0;
--dg:#3f9d5a;--do:#d97706;--dx:#9aa1a9;--dr:#c2410c}
@media (prefers-color-scheme: dark){:root{--bg:#14181d;--fg:#e5e7eb;--muted:#9aa1a9;
--card:#1c2229;--line:#303842;--ok-bg:#16291d;--ok-line:#2f6b43;--ok-fg:#8fd3a3;
--warn-bg:#33260f;--warn-line:#d97706;--warn-fg:#f5c27a;--grey-bg:#222931;
--grey-fg:#7d8590;--note:#2a2415;--dg:#5ec17b;--do:#f0a23b;--dx:#6b7280;--dr:#f97316}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,sans-serif}
main{max-width:900px;margin:0 auto;padding:16px}
h1{margin:0 0 2px;font-size:20px}h2{margin:0 0 6px;font-size:15px}
.muted{color:var(--muted);font-size:12px}
a{color:inherit}
.banner{background:var(--warn-bg);border:2px solid var(--warn-line);color:var(--warn-fg);
border-radius:6px;padding:8px 12px;margin:8px 0}
ul.notes{margin:6px 0;padding-left:18px;font-size:12px;color:var(--muted)}
.tiles{display:flex;gap:6px;flex-wrap:wrap;margin:12px 0 16px}
.tile{padding:6px 10px;border-radius:6px;font-size:12px;border:1px solid var(--line);
background:var(--card);text-decoration:none}
.tile.nominal{border-color:var(--ok-line);background:var(--ok-bg)}
.tile.abnormal{border-color:var(--warn-line);background:var(--warn-bg);font-weight:600}
.tile.notrun{color:var(--grey-fg);background:var(--grey-bg)}
section.test,section.card{background:var(--card);border:1px solid var(--line);border-radius:6px;
padding:12px;margin-bottom:10px;overflow-wrap:anywhere}
section.abnormal{border:2px solid var(--warn-line);padding:14px}
.quiet{display:grid;grid-template-columns:repeat(auto-fill,minmax(180px,1fr));gap:10px;
align-items:start}
.quiet section.card{margin:0;padding:10px}.quiet h2{font-size:13px}
.quiet .figs .thumb{width:100%;max-width:100%}
.pill{display:inline-block;font-size:11px;padding:1px 8px;border-radius:10px;
margin-left:4px;background:var(--grey-bg);color:var(--muted);vertical-align:1px}
.pill.status.nominal{background:var(--ok-bg);color:var(--ok-fg)}
.pill.status.abnormal,.pill.warn{background:var(--warn-bg);color:var(--warn-fg)}
.figs{display:flex;gap:8px;flex-wrap:wrap;align-items:flex-start;margin:8px 0}
.figs .big{max-width:100%;width:520px}.figs .thumb{width:150px;max-width:45%}
.figs img,.figs video{width:100%;height:auto;border:1px solid var(--line);border-radius:3px;background:#fff}
.zoom{cursor:zoom-in}
.analysis{background:var(--note);border-radius:4px;padding:10px 12px;margin:8px 0;
white-space:pre-wrap;font-size:13px}
.facts{font-size:13px;font-variant-numeric:tabular-nums}
.dots{margin-top:8px;line-height:1}
.dot{display:inline-block;width:8px;height:8px;border-radius:50%;margin-right:2px;
background:var(--dx);border:1px solid var(--dx)}
.dot.ring{background:transparent}
.s-ok{background:var(--dg);border-color:var(--dg)}
.s-suspicious,.s-failed{background:var(--do);border-color:var(--do)}
.s-failed{background:var(--dr);border-color:var(--dr)}
.dot.ring.s-ok,.dot.ring.s-suspicious,.dot.ring.s-failed{background:transparent}
#lb{display:none;position:fixed;inset:0;background:rgba(15,18,22,.88);z-index:9999;
align-items:center;justify-content:center;cursor:zoom-out}
#lb img{max-width:94vw;max-height:92vh;background:#fff;border-radius:4px}
footer{margin-top:18px}
@media (min-width:1000px){body{font-size:16px}main{max-width:1500px;padding:24px 32px}
h1{font-size:26px}h2{font-size:18px}.quiet h2{font-size:16px}.muted,ul.notes{font-size:14px}
.tile{font-size:14px;padding:8px 12px}.pill{font-size:13px}
.quiet{grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:14px}
.figs .big{width:900px}.figs .thumb{width:240px}
.analysis,.facts{font-size:15px}.dot{width:11px;height:11px}
section.test{display:grid;grid-template-columns:minmax(0,900px) minmax(300px,1fr);
column-gap:24px;align-items:start}section.test h2{grid-column:1/-1}
section.test .figs{margin:8px 0}section.test .side{min-width:0}}
"""

JS = """
(function(){var lb=document.getElementById('lb'),im=document.getElementById('lbimg'),
g=[],i=0;
function show(){im.src=g[i].src;lb.style.display='flex'}
function hide(){lb.style.display='none';im.removeAttribute('src')}
document.addEventListener('click',function(ev){var t=ev.target;
if(t.classList&&t.classList.contains('zoom')){
g=[].slice.call(document.querySelectorAll('img.zoom')).filter(function(x){
return x.dataset.group===t.dataset.group});i=g.indexOf(t);show()}
else if(lb.style.display==='flex'){hide()}});
document.addEventListener('keydown',function(ev){if(lb.style.display!=='flex')return;
if(ev.key==='Escape')hide();
else if(ev.key==='ArrowRight'){i=(i+1)%g.length;show()}
else if(ev.key==='ArrowLeft'){i=(i-1+g.length)%g.length;show()}});
})();
"""


def e(x):
    return html.escape(str(x), quote=True)


def _build_id(b):
    return b["id"] if isinstance(b, dict) and type(b.get("id")) is int else None


def _decks(cfg):
    """[(job, deck)] in decks.yml order; cfg is the loaded config."""
    return [(j, d) for j, ds in cfg["jobs"].items() for d in ds]


def _kind(st):
    if st["abnormal"]:
        return "abnormal"
    if st["status"] in hist.NOMINAL:
        return "nominal"
    return "notrun" if st["status"] in hist.NOT_RUN else "neutral"


def _with_tonight(history, tonight, recorded):
    if recorded:
        return list(history)
    out = [n for n in history if n["date"] != tonight["date"]] + [tonight]
    return sorted(out, key=lambda n: n["date"])


def _pills(st, rec, tonight, kind):
    p = [f'<span class="pill status {kind}">{e(st["status"])}</span>']
    if rec and rec.get("own_check_unknown") and st["status"] != "own check unknown":
        p.append('<span class="pill warn">own check unknown</span>')
    if st["abnormal"]:
        p.append(f'<span class="pill warn">since {e(st["since"])}</span>')
    if st["new_tonight"]:
        a, b = st["prev_sha"], tonight["run"].get("sha")
        if a and b and a == b:
            p.append('<span class="pill warn">new tonight, no new commit</span>')
        elif a and b:
            p.append(f'<a class="pill warn" href="{COMPARE}{a}...{b}">new tonight</a>')
        else:
            p.append('<span class="pill warn">new tonight</span>')
    carried = st["abnormal"] and (rec or {}).get("status") not in hist.ABNORMAL
    if carried and st["last_judged"]:
        night_word = "night" if st["nights_since_judged"] == 1 else "nights"
        p.append(
            f'<span class="pill">last judged {e(st["last_judged"])}, '
            f"{e(st['nights_since_judged'])} {night_word} ago</span>"
        )
    return "".join(p)


def _figures(name, job, part, build, figure_url, compact=False):
    """The deck's figures: the first large and the rest as thumbnails, or, compact, the
    first one alone as a thumbnail."""
    if part is None:
        return ""
    names = next(
        (d.get("figures", []) for d in part["summary"]["decks"] if d["name"] == name),
        [],
    )
    names = names[:1] if compact else names
    if not names:
        return ""
    urls = [figure_url(job, build, f) for f in names]
    if None in urls:
        return '<p class="muted">figures unavailable: the build is unknown</p>'
    out = []
    for k, (f, url) in enumerate(zip(names, urls, strict=True)):
        url = e(url)
        size = "big" if k == 0 and not compact else "thumb"
        if f.endswith(".mp4"):
            out.append(
                f'<div class="{size}"><video controls preload="none" src="{url}">'
                f'</video> <a class="muted" href="{url}" download>download</a></div>'
            )
        elif f.endswith(".png"):
            out.append(
                f'<div class="{size}"><img loading="lazy" class="zoom" '
                f'data-group="{e(name)}" src="{url}" alt="{e(f)}"></div>'
            )
    return f'<div class="figs">{"".join(out)}</div>'


def _text(part, name):
    text = ((part or {}).get("sections") or {}).get(name)
    text = text[1] if text else ""
    return text[:MODEL_TEXT_CAP] + " […]" if len(text) > MODEL_TEXT_CAP else text


def _facts(deck, rec, compact=False):
    th = {} if compact else deck.get("theory") or {}
    bits = []
    key = rec.get("key") if rec else None
    if key is not None:
        bits.append(f"{e(deck.get('key_number', 'key'))} = {key:.4g}")
    if "value" in th:
        tol = th.get("deck_tolerance")
        theory = e(th["value"]) + (f" ± {e(tol)}" if tol is not None else "")
        bits.append(f"theory {e(th.get('quantity', ''))}: {theory}")
    lines = [f'<div class="facts">{" · ".join(bits)}</div>'] if bits else []
    if rec and rec.get("reason"):
        lines.append(f'<div class="facts muted">{e(rec["reason"])}</div>')
    if deck.get("own_check") == "none":
        lines.append(
            '<div class="facts muted">the test script checks nothing itself</div>'
        )
    return "".join(lines)


def _dots(name, nights, base, limit):
    out = []
    for n in nights[-limit:] if limit > 0 else []:
        rec = n["tests"].get(name)
        if rec is None:
            continue
        status = rec["status"]
        slug = SLUGS.get(status, "unknown")
        build = n["run"]["builds"].get(rec.get("job"))
        title = [n["date"], status] + (
            [f"{rec['key']:.4g}"] if rec.get("key") is not None else []
        )
        if build and build["start"] < n["date"]:
            title.append(f"reused from {build['start']}")
        title = e(" · ".join(title))
        if build is None or rec.get("expired"):
            out.append(f'<span class="dot ring s-{slug}" title="{title}"></span>')
        else:
            href = e(f"{base}/viewLog.html?buildId={build['id']}&tab=artifacts")
            out.append(f'<a class="dot s-{slug}" title="{title}" href="{href}"></a>')
    return f'<div class="dots">{"".join(out)}</div>' if out else ""


def _section(name, deck, job, st, rec, ctx):
    """A test: the abnormal ones large, with figures and the model's text, every other one
    as a compact card (first figure, key number, reason, dots)."""
    kind = _kind(st)
    compact = kind != "abnormal"
    part = ctx["parts"].get(job)
    build = ctx["tonight"]["run"]["builds"].get(job)
    if _build_id(build) is None:
        build = None
    side = []
    # the model's text wherever it can help: not on a test that looks right or didn't run
    text = _text(part, name)
    if text and kind in ("abnormal", "neutral"):
        side.append(f'<div class="analysis">{e(text)}</div>')
    side.append(_facts(deck, rec, compact))
    side.append(_dots(name, ctx["nights"], ctx["base"], ctx["dots"]))
    body = [
        f"<h2>{ICONS[kind]} {e(name)} {_pills(st, rec, ctx['tonight'], kind)}</h2>",
        _figures(name, job, part, build, ctx["figure_url"], compact),
        f'<div class="side">{"".join(side)}</div>',
    ]
    cls = "card" if compact else "test"
    return f'<section id="{e(name)}" class="{cls} {kind}">{"".join(body)}</section>'


def _header(ctx, cfg, counts, notes):
    t, run = ctx["tonight"], ctx["tonight"]["run"]
    when = e(t["date"])
    if run.get("manual"):
        when += f" · manual run {e(run.get('time', ''))}".rstrip()
    if not ctx["recorded"]:
        when += " · not recorded in the history"
    sha = run.get("sha")
    commit = (
        f'<a href="{COMMIT}{sha}">commit {sha[:7]}</a>' if sha else "commit unknown"
    )
    links = " · ".join(
        f'<a href="{e(ctx["base"])}/viewLog.html?buildId={b["id"]}">{e(j)} build {b["id"]}</a>'
        for j, b in run["builds"].items()
        if _build_id(b) is not None
    )
    head = [
        "<h1>PHARE nightly physics</h1>",
        f'<div class="muted">{when} · {commit}{" · " + links if links else ""}'
        f" · read by {e(cfg.get('model', 'unknown'))} · advice only, never a gate</div>",
        f'<div class="muted">{counts}</div>',
    ]
    if notes:
        head.append(
            '<ul class="notes">' + "".join(f"<li>{e(n)}</li>" for n in notes) + "</ul>"
        )
    for part in ctx["parts"].values():
        rc = part["summary"].get("decks_exit_code")
        if rc not in (0, None):
            head.append(
                f'<div class="banner">{e(part["job"])} tests step failed '
                f"(exit {e(rc)}): own checks unknown</div>"
            )
    return "".join(head)


def render(history, tonight, recorded, cfg, parts, notes, figure_url):
    nights = _with_tonight(history, tonight, recorded)
    ctx = {
        "tonight": tonight,
        "recorded": recorded,
        "nights": nights,
        "parts": {p["job"]: p for p in parts},
        "figure_url": figure_url,
        "base": cfg["teamcity"]["url"],
        "dots": cfg["page"]["dots"],
    }
    rows = []
    for job, deck in _decks(cfg):
        name = deck["name"]
        rec = tonight["tests"].get(name)
        st = hist.state(nights, name, deck.get("judge") is not False)
        rows.append((name, deck, job, st, rec))
    kinds = [_kind(r[3]) for r in rows]
    counts = (
        f"{kinds.count('abnormal')} to look at · {kinds.count('nominal')} fine · "
        f"{kinds.count('neutral')} can't tell · {kinds.count('notrun')} not run"
    )
    tiles = "".join(
        f'<a class="tile {_kind(st)}" href="#{e(name)}">{ICONS[_kind(st)]} {e(name)}</a>'
        for name, _, _, st, _ in rows
    )
    shown = [
        (_section(name, deck, job, st, rec, ctx), st)
        for name, deck, job, st, rec in rows
    ]
    sections = "".join(h for h, st in shown if st["abnormal"])
    quiet = "".join(h for h, st in shown if not st["abnormal"])
    if quiet:
        sections += f'<div class="quiet">{quiet}</div>'
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>PHARE nightly physics</title><style>{CSS}</style></head><body><main>"
        f"{_header(ctx, cfg, counts, notes)}"
        f'<div class="tiles">{tiles}</div>{sections}'
        '<footer class="muted">Rewritten every night by the TeamCity shadow nightly. '
        "Figures load from TeamCity: if it is down, they are missing.</footer></main>"
        '<div id="lb"><img id="lbimg" alt=""></div>'
        f"<script>{JS}</script></body></html>"
    )
