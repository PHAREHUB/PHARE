#!/usr/bin/env python3
"""Numbers quoted next to the 'what each method does at a jump' figure of numerics/mhd_reconstruction.
run from doc/source/_ext:  python3 ../validation/scripts/mhd_reconstruction_numbers.py"""
import sys
sys.path.insert(0, ".")
import numpy as np
import mhd_recon as m

n = 30
x = np.arange(n) + 0.5
u = 0.2 + 0.6 * np.exp(-((x - 8) / 3.0) ** 2) + np.where(x > 20, 0.5, 0.0)
print("highest cell value %.4f" % u.max())
for name, (L, R, g) in m.METHODS.items():
    faces = [v for j in range(4, n - 4) for v in (L(u, j), R(u, j))]
    jump = [v for j in range(17, 25) for v in (L(u, j), R(u, j))]
    print(f"{name:18s} highest face value {max(faces):.4f} | max |L-R| on the bump {max(abs(L(u, j) - R(u, j)) for j in range(4, 14)):.4f}"
          f" | at the jump {max(abs(L(u, j) - R(u, j)) for j in range(18, 24)):.4f} | jump-region faces in [{min(jump):.4f}, {max(jump):.4f}]")
