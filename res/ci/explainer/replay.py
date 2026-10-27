#!/usr/bin/env python3
"""Replay the explainer on past reds, offline: nothing is posted, nothing in CI changes.

For every red PR build since SINCE: sort it, cut the excerpt, and for the reds a PR caused run the
agent once per distinct failure (the same error on 3 configurations is explained once). Writes, in OUT:
  would_post.md   the comments it would have posted, newest commit of each PR
  grading.md      one row per explained failure, to mark right / useless / wrong
  summary.json    every red with its verdict
  cache/, work/   logs, and the agent's folders (context.md, src/, raw output)

usage: replay.py SINCE OUT [--model M] [--limit N] [--only ID,ID] [--no-agent] [--repo PATH] [--remote NAME]
       SINCE like 20261019 (00:00 UTC)
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import agent
import comment
import explain
import tc


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("since")
    parser.add_argument("out", type=Path)
    parser.add_argument("--model", default=agent.MODEL)
    parser.add_argument(
        "--limit", type=int, default=0, help="explain at most N failures (0: all)"
    )
    parser.add_argument(
        "--only",
        type=lambda s: {int(x) for x in s.split(",")},
        default=set(),
        help="explain only the failures containing these build ids",
    )
    parser.add_argument(
        "--no-agent", action="store_true", help="sort and cut excerpts only"
    )
    parser.add_argument(
        "--repo", type=Path, default=Path(__file__).resolve().parents[3]
    )
    parser.add_argument("--remote", default="upstream")
    args = parser.parse_args()

    cache = args.out / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    builds = tc.red_builds(f"{args.since}T000000+0000")

    def log_of(build_id):
        path = cache / f"{build_id}.log"
        if not path.exists():
            path.write_text(tc.log(build_id))
        return path.read_text()

    rows, groups, others = explain.sort(builds, log_of)
    (args.out / "summary.json").write_text(json.dumps(rows, indent=1))

    signature = explain.signature(args.model)
    explained = defaultdict(list)  # (pr, sha) -> [(builds, verdict, report)]
    grading = [
        "| PR | commit | builds | kind | owner (agent) | what | grade (right / useless / wrong) |",
        "|---|---|---|---|---|---|---|",
    ]
    for n, ((pr, sha, _), members) in enumerate(
        sorted(groups.items(), key=lambda item: item[1][0][0]["id"])
    ):
        first, verdict, _ = members[0]
        report = None
        wanted = not args.only or any(b["id"] in args.only for b, _, _ in members)
        if not args.no_agent and wanted and (not args.limit or n < args.limit):
            report = explain.explain(
                members,
                args.repo,
                args.remote,
                args.out / "work" / f"{pr}-{sha[:7]}-{first['id']}",
                args.model,
                tc.builds_of(sha),
            )
            print(
                f"{'ok ' if report else 'BAD'} #{pr} {sha[:7]} {verdict.kind} {first['id']}",
                file=sys.stderr,
            )
        explained[(pr, sha)].append(([b for b, _, _ in members], verdict, report))
        what = next(
            (
                line[5:].strip()
                for line in (report or "").splitlines()
                if line.startswith("What:")
            ),
            "",
        )
        owner = next(
            (
                line[6:].strip()
                for line in (report or "").splitlines()
                if line.startswith("Owner:")
            ),
            "",
        )
        grading.append(
            f"| #{pr} | {sha[:7]} | {', '.join(str(b['id']) for b, _, _ in members)} | {verdict.kind} | {owner} | {what.replace('|', '/')} | |"
        )

    newest = {}
    for pr, sha in list(explained) + list(others):
        ids = [b["id"] for b in builds if b["pr"] == pr and b["sha"] == sha]
        if pr not in newest or max(ids) > newest[pr][1]:
            newest[pr] = (sha, max(ids))
    posts = [
        comment.render(
            pr, sha, explained.get((pr, sha), []), others.get((pr, sha), []), signature
        )
        for pr, (sha, _) in sorted(newest.items())
    ]
    (args.out / "would_post.md").write_text(
        "\n\n---\n\n".join(
            f"## PR #{pr}\n\n{post}"
            for (pr, _), post in zip(sorted(newest.items()), posts)
        )
    )
    (args.out / "grading.md").write_text("\n".join(grading) + "\n")
    print(
        f"{len(builds)} reds, {len(groups)} distinct failures from PRs, {len(newest)} PRs -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
