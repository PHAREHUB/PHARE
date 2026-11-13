#!/usr/bin/env python3

import unittest

import matplotlib

matplotlib.use("Agg")

import numpy as np

from pyphare.pharesee.particles import Particles
from pyphare.pharesee.plotting import dist_plot


def particles(n=1000, seed=13337):
    rng = np.random.default_rng(seed)
    return Particles(
        icells=rng.integers(0, 10, size=(n, 1)),
        deltas=rng.random((n, 1)),
        v=rng.normal(loc=(0.5, -0.25, 0.0), size=(n, 3)),
        weights=rng.random((n, 1)),  # (N, 1) as from C++
        dl=np.full((n, 1), 0.1),
    )


def bulk(parts, ivel):
    return np.average(parts.v[:, ivel], weights=parts.weights[:, 0])


class DistPlotBulkTest(unittest.TestCase):
    def test_bulk_lines_two_velocity_axes(self):
        parts = particles()
        _, ax = dist_plot(parts, axis=("Vx", "Vy"), bulk=True)
        vline, hline = ax.lines
        self.assertAlmostEqual(vline.get_xdata()[0], bulk(parts, 0))
        self.assertAlmostEqual(hline.get_ydata()[0], bulk(parts, 1))

    def test_bulk_line_position_and_velocity_axes(self):
        parts = particles()
        _, ax = dist_plot(parts, axis=("x", "Vz"), bulk=True)
        (hline,) = ax.lines
        self.assertAlmostEqual(hline.get_ydata()[0], bulk(parts, 2))

    def test_bulk_line_one_velocity_axis(self):
        parts = particles()
        _, ax = dist_plot(parts, axis=("Vx",), bulk=True)
        vline = ax.lines[-1]  # after the histogram
        self.assertAlmostEqual(vline.get_xdata()[0], bulk(parts, 0))


if __name__ == "__main__":
    unittest.main()
