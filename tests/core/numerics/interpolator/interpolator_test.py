import os
import sys
import numpy as np


def start_index(pos, order, centering):
    assert centering in ["primal", "dual"]

    cell = int(pos)
    delta = pos - cell
    if order == 1:
        if centering == "primal":
            return cell
        else:
            less = 1 if delta < 0.5 else 0
            return cell - less

    if order == 2:
        if centering == "primal":
            less = 1 if delta < 0.5 else 0
            return cell - less
        else:
            return cell - 1

    if order == 3:
        if centering == "primal":
            return cell - 1
        else:
            less = 2 if delta < 0.5 else 1
            return cell - less


def get_nodes(pos, order, centering):
    return start_index(pos, order, centering) + np.arange(order + 1)


def bspline(x, n):
    from scipy.interpolate import BSpline

    knots = np.arange(-(n + 1) / 2, (n + 3) / 2)
    out = BSpline.basis_element(knots)(x)
    out[(x < knots[0]) | (x > knots[-1])] = 0.0
    return out


def main():
    orders = [1, 2, 3]
    particle_positions = 3.0 + np.arange(0, 1.0, 0.1)

    bsplines = []
    path = sys.argv[1]

    for centering in ["primal", "dual"]:
        for order in orders:
            filename = os.path.join(path, f"bsplines_{order}_{centering}.dat")

            node_splines = np.zeros(
                (order + 1, particle_positions.size), dtype=np.int32
            )
            splines_val = np.zeros(
                (order + 1, particle_positions.size), dtype=np.float64
            )

            for ipos, pos in enumerate(particle_positions):
                node_splines[:, ipos] = get_nodes(pos, order, centering)

            # loop on all nodes and calculate the spline at the node
            # from the current particle position

            with open(filename, "wb") as f:
                for ipos, pos in enumerate(particle_positions):
                    if centering == "dual":
                        pos -= 0.5
                    splines_val[:, ipos] = bspline(node_splines[:, ipos] - pos, order)
                    print(ipos, node_splines[:, ipos], splines_val[:, ipos])
                    node_splines[:, ipos].tofile(f)
                    splines_val[:, ipos].tofile(f)


def time_interpolate(before_time, after_time, interp_time, before_data, after_data):
    from pyphare.core.phare_utilities import FloatingPoint_comparator as FP_cmp

    assert before_data.shape == after_data.shape
    assert before_time < after_time
    assert FP_cmp(before_time) <= FP_cmp(interp_time) <= FP_cmp(after_time)

    alpha = (interp_time - before_time) / (after_time - before_time)
    return (1.0 - alpha) * before_data + alpha * after_data


if __name__ == "__main__":
    main()
