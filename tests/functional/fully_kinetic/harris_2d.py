#!/usr/bin/env python3
"""
Fully kinetic 2D Harris sheet, one level, order 1, PHARE units (B0, n0, V_A, d_i, 1/Omega_ci).
Not a test (not registered): the setup and its cost, before a production run.

Periodic in x and y (the FullyKineticModel takes periodic domains only,
fully_kinetic_model.hpp:87-94), so two sheets, at y = Ly/4 and 3 Ly/4, with opposite currents,
as the hybrid tests/functional/harris/harris_2d.py.

Exact Harris equilibrium, four populations:
- sheet ions and electrons, density s1 + s2 (s = sech^2((y - y0) / L)), uniform temperatures,
  drifting along z: V_i / T_i = -V_e / T_e and V_i - V_e = 1 / L (J = curl B, n0 = 1), so the
  ions carry T_i / (T_i + T_e) of the current. The drift flips with the sheet.
- background ions and electrons, density n_b, at rest, same temperatures.
Pressure balance n0 (T_i + T_e) = B0^2 / 2 sets T_i + T_e = 1/2. E = 0 is the equilibrium
field in this frame and the loader starts with E = 0.

B_x = tanh((y - y1) / L) - tanh((y - y2) / L) - 1 plus a localized X-point perturbation on
each sheet at x = Lx / 2 (as the hybrid script): dB = z x grad(psi), psi = +-db sigma
exp(-r^2 / sigma^2) summed over its periodic images, divergence free and periodic. Its curl
has no particle current: E answers it with a light and plasma wave transient at t = 0.

dt = courant dx / (c sqrt(2)) rounded down to 1e-6, the light CFL c dt sqrt(sum 1/dx^2) of the setup checks
(pyphare/pharein/fully_kinetic.py).

Checks, at step 0 and at the last step n, the run doing n + 1 steps (particle dumps of every
population at both, E at steps n and n + 1):
- energy W = sum_p m w |v|^2 / 2 dV + sum_nodes (E^{n-1/2} . E^{n+1/2} / c^2 + |B|^2) / 2 dV
  (energy.py), relative drift
- Gauss: F = div E - c^2 F^N rho, rho of the particles at x - (dt/2) v (E is dumped half a step
  behind), filtered as the current; with cleaning 0 the deposit conserves charge, so F stays at
  its t = 0 value (-c^2 rho of the loading noise): reports max |F(end) - F(0)| / c^2 max |F^N rho(0)|
  as tests/simulator/test_fully_kinetic.py.

Timing: wall time of every advance, a progress line every 100 steps; steps that dump are listed
apart at the end. Peak memory per rank.
Runs under MPI (mpirun -n N); the checks run on rank 0 after the simulation is reset.
"""

import argparse
import resource
import time

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

LX, LY = 25.6, 12.8
HALF_WIDTH = 0.5  # L, in d_i
N_BACKGROUND = 0.2
TI_OVER_TE = 5.0
PERTURBATION = 0.1  # db, in B0
SIGMA = 1.0  # d_i
END_TIME = 30.0
NBR_DUMPS = 10
PROGRESS = 100  # steps between progress lines on rank 0


def parse_args():
    parser = argparse.ArgumentParser(description="FullyKineticModel 2D Harris sheet")
    parser.add_argument("--steps", type=int, help="time steps (default: up to t = 30)")
    parser.add_argument("--ppc", type=int, default=100, help="particles per cell and population")
    parser.add_argument("--mass-ratio", type=float, default=25.0, help="m_i / m_e")
    parser.add_argument("--dx", type=float, default=0.05, help="cell size, both directions, d_i")
    parser.add_argument("--c-over-va", type=float, default=16.0)
    parser.add_argument("--courant", type=float, default=0.75, help="c dt sqrt(2) / dx")
    parser.add_argument("--filter-passes", type=int, default=1)
    parser.add_argument("--cleaning", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--no-dumps", action="store_true", help="no field and moment dumps (the checks still dump)"
    )
    parser.add_argument("--no-checks", action="store_true", help="no energy and Gauss checks")
    parser.add_argument("--diag-dir", default="phare_outputs/fully_kinetic_harris_2d")
    return parser.parse_args()


class Setup:
    def __init__(self, args):
        self.args = args
        self.dl = args.dx
        self.cells = (round(LX / args.dx), round(LY / args.dx))
        self.c = args.c_over_va
        # rounded down to 1e-6: the dump times must be multiples of dt to 1e-9
        self.dt = np.floor(1e6 * args.courant * args.dx / (self.c * np.sqrt(2))) / 1e6
        self.steps = args.steps or round(END_TIME / self.dt)
        self.mass_e = 1.0 / args.mass_ratio
        self.Te = 0.5 / (1 + TI_OVER_TE)
        self.Ti = 0.5 - self.Te
        # V_i - V_e = 1 / L, split as the temperatures
        self.Vi = self.Ti / (self.Ti + self.Te) / HALF_WIDTH
        self.Ve = -self.Te / (self.Ti + self.Te) / HALF_WIDTH
        self.masses = {"protons": 1.0, "electrons": self.mass_e}
        self.masses = {
            f"{kind}_{name}": m
            for kind in ("sheet", "background")
            for name, m in self.masses.items()
        }
        self.charges = {pop: (1.0 if "protons" in pop else -1.0) for pop in self.masses}

    def time(self, step):
        # as the simulation computes its final time, never past it
        return step * self.dt

    def check_steps(self):
        if self.args.no_checks:
            return []
        return [0, 1, self.steps, self.steps + 1]

    def print_table(self):
        n_peak = 1 + N_BACKGROUND
        vthe = np.sqrt(self.Te / self.mass_e)
        wpe = lambda n: self.c * np.sqrt(n / self.mass_e)  # noqa: E731
        debye = vthe / wpe(n_peak)
        print(
            f"cells {self.cells}, dx {self.dl}, dt {self.dt:.6g}, steps {self.steps}"
            f" (t = {self.steps * self.dt:.3g}), light CFL {self.c * self.dt * np.sqrt(2) / self.dl:.3f}\n"
            f"d_e {np.sqrt(self.mass_e / n_peak):.3f} (peak n {n_peak:g}),"
            f" {np.sqrt(self.mass_e / N_BACKGROUND):.3f} (lobes); lambda_De {debye:.4f} (peak),"
            f" {vthe / wpe(N_BACKGROUND):.4f} (lobes); dx / lambda_De {self.dl / debye:.2f}\n"
            f"w_pe dt {wpe(n_peak) * self.dt:.3f}, Omega_ce dt {(1 + PERTURBATION) / self.mass_e * self.dt:.3f},"
            f" T_i {self.Ti:.4f} T_e {self.Te:.4f}, V_i {self.Vi:+.3f} V_e {self.Ve:+.3f} (sheet 1),"
            f" v_th,e / c {vthe / self.c:.3f}"
        )


def config(setup):
    args = setup.args
    Lx, Ly = np.asarray(setup.cells) * setup.dl
    y1, y2 = 0.25 * Ly, 0.75 * Ly
    L = HALF_WIDTH

    # each sheet with its periodic images at +-Ly: density, current and curl B stay
    # consistent and periodic (without them the drift flips sign at y = 0)
    images = (-Ly, 0.0, Ly)

    def sech2(y, y0):
        return sum(1.0 / np.cosh((y - y0 - shift) / L) ** 2 for shift in images)

    def tanh(y, y0):
        return sum(np.tanh((y - y0 - shift) / L) for shift in images)

    def perturbation(x, y):
        """(dB_x, dB_y) = z x grad(psi), psi = db sigma exp(-r^2 / sigma^2) centred on (Lx/2, y1)
        minus the same on (Lx/2, y2), each with its periodic images at +-Lx and +-Ly: without
        them B_y jumps at the y border (exp(-(Ly/4)^2 / sigma^2) is 4e-5 here) and the copies
        of a shared border node differ from t = 0"""
        dbx, dby = 0.0 * x, 0.0 * x
        for y0, sign in ((y1, 1.0), (y2, -1.0)):
            for xshift in (-Lx, 0.0, Lx):
                for yshift in images:
                    x0, yy0 = x - 0.5 * Lx - xshift, y - y0 - yshift
                    g = sign * 2 * PERTURBATION / SIGMA * np.exp(-(x0**2 + yy0**2) / SIGMA**2)
                    dbx = dbx - yy0 * g  # dpsi/dy
                    dby = dby + x0 * g  # -dpsi/dx
        return dbx, dby

    def bx(x, y):
        return tanh(y, y1) - tanh(y, y2) - 1 + perturbation(x, y)[0]

    def by(x, y):
        return perturbation(x, y)[1]

    def bz(x, y):
        return 0.0 * x

    def sheet_density(x, y):
        return sech2(y, y1) + sech2(y, y2) + 0.0 * x

    def background_density(x, y):
        return N_BACKGROUND + 0.0 * x

    def drift(V):
        """sheet 1 carries J_z = -(1/L) sech^2, sheet 2 the opposite"""

        def fn(x, y):
            s1, s2 = sech2(y, y1), sech2(y, y2)
            return -V * (s1 - s2) / np.maximum(s1 + s2, 1e-300) + 0.0 * x

        return fn

    def uniform(value):
        return lambda x, y: value + 0.0 * x

    def population(kind, name, seed):
        mass = 1.0 if name == "protons" else setup.mass_e
        T = setup.Ti if name == "protons" else setup.Te
        vth = uniform(np.sqrt(T / mass))
        sheet = kind == "sheet"
        V = setup.Vi if name == "protons" else setup.Ve
        return {
            "charge": 1 if name == "protons" else -1,
            "mass": mass,
            "density": sheet_density if sheet else background_density,
            "vbulkx": uniform(0.0),
            "vbulky": uniform(0.0),
            "vbulkz": drift(V) if sheet else uniform(0.0),
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": args.ppc,
            "init": {"seed": seed + 10 * args.seed},
        }

    ph.global_vars.sim = None
    check_steps = setup.check_steps()
    sim = ph.Simulation(
        time_step=setup.dt,
        # E^{n+1/2} of the last check step is dumped one step later
        time_step_nbr=setup.steps + (1 if check_steps else 0),
        cells=setup.cells,
        dl=[setup.dl] * 2,
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=setup.c,
        current_filter_passes=args.filter_passes,
        divergence_cleaning=args.cleaning,
        diag_options={
            "format": "phareh5",
            "options": {"dir": args.diag_dir, "mode": "overwrite"},
        },
        # dx / lambda_De is 3.0 at the sheet peak at c/V_A = 16, on the guard: warned, not refused
        strict=False,
    )
    pops = {
        f"{kind}_{name}": population(kind, name, seed)
        for seed, (kind, name) in enumerate(
            (k, n) for k in ("sheet", "background") for n in ("protons", "electrons")
        )
    }
    ph.MaxwellianFluidModel(bx=bx, by=by, bz=bz, **pops)

    dump_steps = sorted(set(check_steps))
    if not args.no_dumps:
        dump_steps = sorted(set(dump_steps) | set(np.linspace(0, setup.steps, NBR_DUMPS + 1).round().astype(int)))
        moments = sorted(set(np.linspace(0, setup.steps, NBR_DUMPS + 1).round().astype(int)))
        moment_times = [setup.time(s) for s in moments]
        for pop in pops:
            for quantity in ("density", "flux"):
                ph.FluidDiagnostics(
                    quantity=quantity, write_timestamps=moment_times, population_name=pop
                )
    if dump_steps:
        for quantity in ("E", "B"):
            ph.ElectromagDiagnostics(
                quantity=quantity, write_timestamps=[setup.time(s) for s in dump_steps]
            )
    if check_steps:
        for pop in pops:
            ph.ParticleDiagnostics(
                quantity="domain",
                write_timestamps=[setup.time(s) for s in (0, setup.steps)],
                population_name=pop,
            )
    return sim, dump_steps


def grid(setup, quantity, component, time):
    """the periodic level 0 values of a field component, ghosts dropped, on a (2 nx, 2 ny) grid
    of half-cell keys (primal 2i, dual 2i + 1); returns the slice of its centering"""
    n = np.asarray(setup.cells)
    out = np.full(2 * n, np.nan)
    hier = hierarchy_from(
        h5_filename=f"{setup.args.diag_dir}/EM_{quantity}.h5", times=f"{time:.10f}"
    )
    for patch in hier.level(0, time).patches:
        pdata = patch.patch_datas[f"{quantity}{component}"]
        mesh = pdata.meshgrid()
        values = np.asarray(pdata.dataset[:])
        keys = [np.round(2 * x / setup.dl).astype(int) for x in mesh]
        inside = np.all([(k >= 0) & (k < 2 * m) for k, m in zip(keys, n)], axis=0)
        out[keys[0][inside], keys[1][inside]] = values[inside]
    # Yee: E_d dual along d, B_d primal along d and dual along the others
    parity = [int((component == c) == (quantity == "E")) for c in "xy"]
    field = out[parity[0] :: 2, parity[1] :: 2]
    assert not np.isnan(field).any(), (quantity, component, time)
    return field


def particles(setup, pop, time):
    """positions in cells (iCell + delta), velocities, weights of the domain particles"""
    import h5py

    xs, vs, ws = [], [], []
    with h5py.File(f"{setup.args.diag_dir}/ions_pop_{pop}_domain.h5", "r") as f:
        key = next(k for k in f["t"] if abs(float(k) - time) < 1e-9)
        for patch in f["t"][key]["pl0"].values():
            if "iCell" not in patch:
                continue
            xs.append(patch["iCell"][:, :2] + patch["delta"][:, :2].astype(np.float64))
            vs.append(patch["v"][:].astype(np.float64))
            ws.append(patch["weight"][:].ravel())
    return np.concatenate(xs), np.concatenate(vs), np.concatenate(ws)


def deposit(setup, pos, weights):
    """order 1 (cloud in cell) on the periodic primal nodes"""
    n = np.asarray(setup.cells)
    left = np.floor(pos).astype(np.int64)
    frac = pos - left
    rho = np.zeros(n[0] * n[1])
    for ox in (0, 1):
        for oy in (0, 1):
            w = weights * (frac[:, 0] if ox else 1 - frac[:, 0]) * (frac[:, 1] if oy else 1 - frac[:, 1])
            index = ((left[:, 0] + ox) % n[0]) * n[1] + (left[:, 1] + oy) % n[1]
            rho += np.bincount(index, weights=w, minlength=rho.size)
    return rho.reshape(n)


def binomial(field, passes):
    for _ in range(passes):
        for axis in (0, 1):
            field = 0.25 * np.roll(field, 1, axis) + 0.5 * field + 0.25 * np.roll(field, -1, axis)
    return field


def kinetic_and_rho(setup, time):
    kinetic, rho, count = 0.0, 0.0, 0
    for pop, mass in setup.masses.items():
        x, v, w = particles(setup, pop, time)
        count += len(w)
        kinetic += 0.5 * mass * np.sum(w * np.sum(v**2, axis=1))
        rho = rho + deposit(setup, x - 0.5 * setup.dt * v[:, :2] / setup.dl, setup.charges[pop] * w)
        del x, v, w
    return kinetic * setup.dl**2, binomial(rho, setup.args.filter_passes), count


def field_energy(setup, step):
    t0, t1 = setup.time(step), setup.time(step + 1)
    electric = sum(np.sum(grid(setup, "E", c, t0) * grid(setup, "E", c, t1)) for c in "xyz")
    magnetic = sum(np.sum(grid(setup, "B", c, t0) ** 2) for c in "xyz")
    return electric / setup.c**2 / 2 * setup.dl**2, magnetic / 2 * setup.dl**2


def gauss(setup, step, rho):
    t = setup.time(step)
    ex, ey = grid(setup, "E", "x", t), grid(setup, "E", "y", t)
    div = (ex - np.roll(ex, 1, 0)) / setup.dl + (ey - np.roll(ey, 1, 1)) / setup.dl
    return div - setup.c**2 * rho


def checks(setup):
    first, last = 0, setup.steps
    results = {}
    for step in (first, last):
        kinetic, rho, count = kinetic_and_rho(setup, setup.time(step))
        electric, magnetic = field_energy(setup, step)
        results[step] = (kinetic, electric, magnetic, gauss(setup, step, rho), rho, count)
        print(
            f"step {step}: particles {count}, kinetic {kinetic:.6e} electric {electric:.3e}"
            f" magnetic {magnetic:.6e} total {kinetic + electric + magnetic:.6e}"
        )
    W = {s: sum(r[:3]) for s, r in results.items()}
    scale = setup.c**2 * np.abs(results[first][4]).max()
    drift = np.abs(results[last][3] - results[first][3]).max() / scale
    print(
        f"energy W / W(0) - 1 = {(W[last] - W[first]) / W[first]:+.3e} over {last} steps;"
        f" Gauss residual: initial {np.abs(results[first][3]).max() / scale:.3e},"
        f" drift {drift:.3e} (of c^2 max |F^N rho(0)|)"
    )


def main():
    args = parse_args()
    startMPI()
    setup = Setup(args)
    if cpp.mpi_rank() == 0:
        setup.print_table()
    sim, dump_steps = config(setup)
    simulator = Simulator(sim)
    tick = time.perf_counter()
    simulator.initialize()
    init = time.perf_counter() - tick
    walls = []
    for step in range(1, sim.time_stepper.time_step_nbr + 1):
        tick = time.perf_counter()
        simulator.advance()
        walls.append(time.perf_counter() - tick)
        if step % PROGRESS == 0 and cpp.mpi_rank() == 0:
            last = np.array(walls[-PROGRESS:])
            print(
                f"step {step} t {setup.time(step):.4f}: last {PROGRESS} steps median {np.median(last):.4f} s,"
                f" mean {last.mean():.4f} s; elapsed {sum(walls) / 3600:.3f} h",
                flush=True,
            )
    walls = np.array(walls)
    simulator.reset()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
    print(f"rank {cpp.mpi_rank()}: peak RSS {peak:.2f} GB", flush=True)
    cpp.mpi_barrier()
    if cpp.mpi_rank() == 0:
        quiet = np.ones(len(walls), bool)
        quiet[[s - 1 for s in dump_steps if s > 0]] = False
        line = f"ranks {cpp.mpi_size()}: init {init:.1f} s; {len(walls)} steps {walls.sum():.1f} s"
        if quiet.any():
            q = walls[quiet]
            line += (
                f"; steps without dump ({quiet.sum()}): mean {q.mean():.4f} s, median"
                f" {np.median(q):.4f} s, min {q.min():.4f} s, max {q.max():.4f} s"
            )
        print(f"{line}; dump steps {walls[~quiet].round(2).tolist()}")
        if not args.no_checks:
            checks(setup)
    cpp.mpi_barrier()


if __name__ == "__main__":
    main()
