"""Size formatting for reports of what a read touches (see ``Run.plan``)."""

from __future__ import annotations


def format_bytes(n):
    """Short human-readable size, e.g. ``'1.50 GiB'``."""
    n = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.2f} {unit}"
        n /= 1024
    return f"{n:.2f} TiB"
