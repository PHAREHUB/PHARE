#!/usr/bin/env python3
import re
import ast
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "perms"

records = []
with open(path) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        m = re.search(
            r"N particles:(\d+).*?MEM:(\d+)mb.*?ppp time:([\d.]+)ns\s+(\{.*\})", line
        )
        if not m:
            continue
        n_particles = int(m.group(1))
        mem_mb = int(m.group(2))
        ppp_ns = float(m.group(3))
        params = ast.literal_eval(m.group(4))
        records.append(
            {
                "tile_size": params["tile_size"],
                "tag_buffer": params["tag_buffer"],
                "tagging_threshold": params["tagging_threshold"],
                "ppp_ns": ppp_ns,
                "mem_mb": mem_mb,
                "n_particles": n_particles,
            }
        )

tile_sizes = sorted(set(r["tile_size"] for r in records))
tag_buffers = sorted(set(r["tag_buffer"] for r in records))
thresholds = sorted(set(r["tagging_threshold"] for r in records))

# ── Figure 1: ppp_time vs tile_size, lines per tag_buffer, subplots per threshold ──
fig, axes = plt.subplots(1, len(thresholds), figsize=(14, 4), sharey=True)
fig.suptitle("Push time (ns/particle) vs tile_size", fontsize=13)

colors = cm.tab10(np.linspace(0, 0.4, len(tag_buffers)))

for ax, thresh in zip(axes, thresholds):
    for color, tb in zip(colors, tag_buffers):
        pts = sorted(
            [
                r
                for r in records
                if r["tag_buffer"] == tb and r["tagging_threshold"] == thresh
            ],
            key=lambda r: r["tile_size"],
        )
        ax.plot(
            [r["tile_size"] for r in pts],
            [r["ppp_ns"] for r in pts],
            marker="o",
            color=color,
            label=f"tb={tb}",
        )
    ax.set_title(f"threshold={thresh}")
    ax.set_xlabel("tile_size")
    ax.set_xticks(tile_sizes)
    ax.grid(True, alpha=0.3)

axes[0].set_ylabel("ppp time (ns)")
axes[-1].legend(title="tag_buffer", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("perms_ppp_lines.png", dpi=120, bbox_inches="tight")
print("saved perms_ppp_lines.png")

# ── Figure 2: heatmaps ppp_time [tile_size × tag_buffer] per threshold ──
fig, axes = plt.subplots(1, len(thresholds), figsize=(14, 3.5))
fig.suptitle("Heatmap: ppp time (ns) — tile_size × tag_buffer", fontsize=13)

for ax, thresh in zip(axes, thresholds):
    mat = np.zeros((len(tile_sizes), len(tag_buffers)))
    for i, ts in enumerate(tile_sizes):
        for j, tb in enumerate(tag_buffers):
            row = next(
                r
                for r in records
                if r["tile_size"] == ts
                and r["tag_buffer"] == tb
                and r["tagging_threshold"] == thresh
            )
            mat[i, j] = row["ppp_ns"]

    im = ax.imshow(mat, aspect="auto", cmap="RdYlGn_r")
    ax.set_xticks(range(len(tag_buffers)))
    ax.set_xticklabels(tag_buffers)
    ax.set_yticks(range(len(tile_sizes)))
    ax.set_yticklabels(tile_sizes)
    ax.set_xlabel("tag_buffer")
    ax.set_ylabel("tile_size")
    ax.set_title(f"threshold={thresh}")
    for i in range(len(tile_sizes)):
        for j in range(len(tag_buffers)):
            ax.text(j, i, f"{mat[i,j]:.0f}", ha="center", va="center", fontsize=7)
    plt.colorbar(im, ax=ax, label="ns")

plt.tight_layout()
plt.savefig("perms_ppp_heatmap.png", dpi=120, bbox_inches="tight")
print("saved perms_ppp_heatmap.png")

# ── Figure 3: MEM and N particles vs tile_size ──
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
fig.suptitle("Memory & particle count vs tile_size", fontsize=13)

for color, tb in zip(colors, tag_buffers):
    for thresh, ls in zip(thresholds, ["-", "--", "-.", ":"]):
        pts = sorted(
            [
                r
                for r in records
                if r["tag_buffer"] == tb and r["tagging_threshold"] == thresh
            ],
            key=lambda r: r["tile_size"],
        )
        label = f"tb={tb} tt={thresh}" if thresh == thresholds[0] else None
        axes[0].plot(
            [r["tile_size"] for r in pts],
            [r["mem_mb"] for r in pts],
            marker=".",
            color=color,
            linestyle=ls,
            alpha=0.7,
            label=label,
        )
        axes[1].plot(
            [r["tile_size"] for r in pts],
            [r["n_particles"] / 1e6 for r in pts],
            marker=".",
            color=color,
            linestyle=ls,
            alpha=0.7,
        )

axes[0].set_xlabel("tile_size")
axes[0].set_ylabel("MEM (MB)")
axes[0].set_title("Memory")
axes[0].grid(True, alpha=0.3)
axes[0].set_xticks(tile_sizes)
axes[1].set_xlabel("tile_size")
axes[1].set_ylabel("N particles (M)")
axes[1].set_title("Particle count")
axes[1].grid(True, alpha=0.3)
axes[1].set_xticks(tile_sizes)
axes[0].legend(fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig("perms_mem_npart.png", dpi=120, bbox_inches="tight")
print("saved perms_mem_npart.png")

# ── Print best configs ──
print("\nTop 10 fastest configs (lowest ppp_ns):")
top = sorted(records, key=lambda r: r["ppp_ns"])[:10]
for r in top:
    print(
        f"  {r['ppp_ns']:.1f} ns  tile_size={r['tile_size']}  tag_buffer={r['tag_buffer']}  threshold={r['tagging_threshold']}  MEM={r['mem_mb']}MB  N={r['n_particles']:,}"
    )
