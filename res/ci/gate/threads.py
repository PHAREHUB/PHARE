#!/usr/bin/env python3
"""Check the PR's review threads: every thread a human started must be resolved.

Threads started by bots (CodeRabbit, Copilot, ...) are ignored: they are advice, not review.

GitHub does not run workflows when a thread is resolved: this status refreshes on the next
push, review, label or edit of the PR, or when the job is re-run.

usage: threads.py OWNER/REPO PR_NUMBER     (token in GITHUB_TOKEN; read access is enough)
exit:  0 no unresolved human thread, 1 otherwise
"""

import json
import os
import sys
import urllib.request

from gitdiff import Finding, report

CHECK = "review-threads"

QUERY = """
query($owner: String!, $repo: String!, $pr: Int!, $after: String) {
  repository(owner: $owner, name: $repo) {
    pullRequest(number: $pr) {
      reviewThreads(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes {
          isResolved
          path
          line
          originalLine
          comments(first: 1) { nodes { url author { login __typename } } }
        }
      }
    }
  }
}
"""


def graphql(query, variables, token):
    request = urllib.request.Request(
        "https://api.github.com/graphql",
        data=json.dumps({"query": query, "variables": variables}).encode(),
        headers={"Authorization": f"bearer {token}", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        answer = json.load(response)
    if answer.get("errors"):
        raise RuntimeError(answer["errors"])
    return answer["data"]


def fetch(owner, repo, pr, token):
    """List of thread nodes."""
    threads, after = [], None
    while True:
        data = graphql(QUERY, {"owner": owner, "repo": repo, "pr": pr, "after": after}, token)
        pull = data["repository"]["pullRequest"]
        page = pull["reviewThreads"]
        threads += page["nodes"]
        if not page["pageInfo"]["hasNextPage"]:
            return threads
        after = page["pageInfo"]["endCursor"]


def is_bot(author):
    return author is None or author.get("__typename") == "Bot" or author.get("login", "").endswith("[bot]")


def evaluate(threads):
    findings = []
    for t in threads:
        first = (t["comments"]["nodes"] or [{}])[0]
        starter = first.get("author")
        if is_bot(starter):
            continue
        where = (t.get("path") or "", t.get("line") or t.get("originalLine") or 0)
        if not t["isResolved"]:
            findings.append(Finding(*where, CHECK, f"unresolved review thread started by {starter['login']}: {first.get('url', '')}"))
    return findings


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    owner, repo = sys.argv[1].split("/")
    threads = fetch(owner, repo, int(sys.argv[2]), os.environ["GITHUB_TOKEN"])
    sys.exit(report(evaluate(threads), CHECK))
