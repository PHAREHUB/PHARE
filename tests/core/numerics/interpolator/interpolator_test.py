
from scipy.signal import bspline
import numpy as np

import sys
import os


def start_index(pos, order):
    return int(pos - (float(order)-1.)/2.)

def get_nodes(pos, order):
    return start_index(pos, order) + np.arange(order+1)


def main():

    orders = [1,2,3]
    particle_positions = 3. + np.arange(0,1.,0.1)

    bsplines = []

    path = sys.argv[1]

    for order in orders:
        filename = path+os.path.sep + "bsplines_{}.dat".format(order)

        node_splines = np.zeros((order+1, particle_positions.size), dtype=np.int32)
        splines_val  = np.zeros((order+1, particle_positions.size), dtype=np.float64)

        for ipos, pos in enumerate(particle_positions):
            node_splines[:,ipos] = get_nodes(pos, order)

        # loop on all nodes and calculate the spline at the node
        # from the current particle position

        with open(filename, 'wb') as f:
            for ipos, pos in enumerate(particle_positions):
                splines_val[:,ipos] = bspline(node_splines[:, ipos] - pos, order)
                print(ipos, node_splines[:, ipos], splines_val[:,ipos])
                node_splines[:, ipos].tofile(f)
                splines_val[:, ipos].tofile(f)



if __name__ == "__main__":
    main()
