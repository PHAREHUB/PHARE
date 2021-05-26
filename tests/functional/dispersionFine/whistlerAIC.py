#!/usr/bin/env python3

from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np



def omega(k, p):
    k2 = k*k
    return 0.5*k2*(np.sqrt(1+4/k2)+p)


rc('text', usetex = True)

fig, ax = plt.subplots(figsize=(4,3), nrows=1)

# theoretical branches
k_the = np.arange(0.04, 20, 0.001)
w_thR = omega(k_the, +1)
w_thL = omega(k_the, -1)

ax.plot(k_the, w_thR, '-k')
ax.plot(k_the, w_thL, '-k')

# numerical branches for 2d
k_numR2d = [0.6275, 1.2551, 2.5101,  5.0203, 10.0406]
w_numR2d = [0.8168, 2.2305, 7.1000, 25.6354, 93.7137]

k_numL2d = [0.6275, 1.2551, 2.5101, 5.0203, 10.0406]
w_numL2d = [0.4398, 0.6597, 0.8482, 0.9111,  0.9739]

ax.plot(k_numR2d, w_numR2d, 'b+', label='R mode', markersize=8)
ax.plot(k_numL2d, w_numL2d, 'rx', label='L mode', markersize=8)

# numerical branches for 1d
k_numR1d = [0.0628, 0.1257, 0.2513, 0.5027, 1.0053, 2.0106,  4.0212,  8.0425]
w_numR1d = [0.0628, 0.1257, 0.2827, 0.6597, 1.6336, 4.8381, 16.8075, 62.1721]

k_numL1d = [0.0628, 0.1257, 0.2513, 0.5027, 1.0053, 2.0106, 4.0212, 8.0425]
w_numL1d = [0.0628, 0.1257, 0.2199, 0.3770, 0.6283, 0.8168, 0.9111, 0.9739]


ax.plot(k_numR1d, w_numR1d, 'b+', markersize=5)
ax.plot(k_numL1d, w_numL1d, 'rx', markersize=5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$k_{\parallel} c / \omega_p$')
ax.set_ylabel('$\omega / \Omega_p$')

ax.legend(loc='upper left', frameon=False)

fig.tight_layout()
fig.savefig("disp.eps")

