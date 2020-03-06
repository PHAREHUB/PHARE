# formatted with black code formatter

import os, h5py, phare.pharein as ph
from .diagnostic.dtype import _EM, _Fluid, _Particle, Particles
from .diagnostic.patch import Patch, PatchLevel


class H5FileRegexer:
    """Used to Validates h5 input file to ignore files which are not to be read/extracted"""

    @staticmethod
    def ignore_population_name_in(h5files):
        import re

        pop = "${POP}"
        regexing = {}
        for file in h5files:
            if pop not in file:
                continue
            bits = file.split(pop)
            left, right = bits[0], bits[1]
            regexing[file] = re.compile("^" + left + "(.*)" + right + "$")
        return regexing


class Diagnostic:
    REFINEMENT_RATIO = 2
    """12 decimal places is arbitrary, but hopefully sufficient"""
    DL_PRECISION = 12

    """List of all possible input HDF5 Files"""
    h5Files = [
        "EM_B.h5",
        "EM_E.h5",
        "ions_bulkVelocity.h5",
        "ions_density.h5",
        "ions_pop_ions_${POP}_density.h5",
        "ions_pop_ions_${POP}_flux.h5",
        "ions_pop_ions_${POP}_domain.h5",
        "ions_pop_ions_${POP}_patchGhost.h5",
        "ions_pop_ions_${POP}_levelGhost.h5",
    ]
    """Compiled Regexs to ignore population name in HDF5 input files"""
    h5FilesRegexes = H5FileRegexer.ignore_population_name_in(h5Files)

    def __init__(self, sim: ph.Simulation, file):
        def truncate(number, digits) -> float:
            import math

            stepper = 10.0 ** digits
            return math.trunc(stepper * number) / stepper

        self.sim = sim
        self.file = file
        self.hifile = h5py.File(file, "r")
        self.levels = {}
        self.dl = [truncate(x, Diagnostic.DL_PRECISION) for x in sim.dl]
        self.dim = len(sim.dl)
        self.initialize_patch_hierarchy_from_hdf5_file()

    def initialize_patch_hierarchy_from_hdf5_file(self):
        for timestamp in self.hifile.keys():
            nbrCells = self.sim.cells
            for patch_level in self.hifile[timestamp].keys():
                lvl = int(patch_level[2:])  # pl123 - get 123
                self.levels[lvl] = PatchLevel(self, lvl, nbrCells)
                for patch in self.hifile[timestamp][patch_level].keys():
                    extracted = Diagnostics._extract(
                        self, self.hifile[timestamp][patch_level][patch]
                    )
                    if extracted is not None:
                        self.levels[lvl].patches.append(
                            Patch(
                                self.levels[lvl],
                                self.hifile[timestamp][patch_level][patch],
                                extracted,
                            )
                        )
                self.levels[lvl].patchDict = {p.id: p for p in self.levels[lvl].patches}
                nbrCells = [v * Diagnostic.REFINEMENT_RATIO for v in nbrCells]

    @staticmethod
    def is_input_file_valid(file):
        return file in Diagnostic.h5Files or any(
            [regex.match(file) for regex in Diagnostic.h5FilesRegexes.values()]
        )


class Diagnostics:
    @staticmethod
    def _extract_type_recursive(diag, h5item, Type):  # "type" is a python3 keyword
        assert isinstance(diag, Diagnostic)

        for key in list(h5item.keys()):
            for k in Type.end_keys:
                if key.endswith(k):
                    return Type(diag, k, h5item[key])
            if not isinstance(h5item[key], h5py.Dataset):
                return Diagnostics._extract_type_recursive(diag, h5item[key], Type)

    @staticmethod
    def _extract_type(diag, h5item, Type):
        assert isinstance(diag, Diagnostic)

        keys = list(h5item.keys())
        if len(keys) and keys[0].startswith(Type.start_key):
            return Diagnostics._extract_type_recursive(diag, h5item, Type)

    @staticmethod
    def _extract(diag, h5item):
        assert isinstance(diag, Diagnostic)

        types = [
            Diagnostics._extract_type(diag, h5item, Type)
            for Type in [_EM, _Fluid, _Particle]
        ]
        types = [i for i in types if i]
        return types[0] if len(types) else None

    @staticmethod
    def _get_h5_files(sim: ph.Simulation):
        assert isinstance(sim, ph.Simulation)

        h5files = []
        for root, dirs, files in os.walk(sim.diag_options["options"]["dir"]):
            for file in files:
                if Diagnostic.is_input_file_valid(file):
                    h5files.append(os.path.join(root, file))
                else:
                    print("Warning: Input file not recognized, ignored: ", file)
        return h5files

    @staticmethod
    def _validate_simulator(sim: ph.Simulation):
        """Fails if the input Simulation object has no diagnostic output"""

        assert isinstance(sim, ph.Simulation)

        if any(
            [
                sim.diag_options is None,
                sim.diag_options["options"] is None
                or sim.diag_options["options"]["dir"] is None,
            ]
        ):
            raise ValueError(
                "Cannot load diagnostics from specific simulation: no diagnostics made"
            )

    @staticmethod
    def _merge_particle_pops(diags):
        """Takes per HDF5 file Diagnostic objects and coallates the particle population
            types into one DType object."""

        ndiags, particle_pops = [], {}
        for d in diags:
            if not isinstance(d.type, _Particle):
                ndiags.append(d)
            else:
                if d.type.pop not in particle_pops:
                    particle_pops[d.type.pop] = {}
                particle_pops[d.type.pop][d.type.name] = d
        for p, t in particle_pops.items():
            ndiags.append(Particles(p, t["domain"], t["patchGhost"], t["levelGhost"]))
        return ndiags

    @staticmethod
    def extract(sim: ph.Simulation):
        assert isinstance(sim, ph.Simulation)

        Diagnostics._validate_simulator(sim)
        diags = [Diagnostic(sim, f) for f in Diagnostics._get_h5_files(sim)]
        if len(diags) == 0:
            raise ValueError(
                "Cannot load diagnostics from specific simulation: no diagnostics found"
            )
        dic = {clazz.__name__: [] for clazz in [_EM, _Fluid, _Particle]}
        [
            dic[type(d.type).__name__].append(d)
            for d in Diagnostics._merge_particle_pops(diags)
        ]
        return dic

    @staticmethod
    def get_type(diags, Type):
        nDiags = []

        if not Type.endswith(".h5"):
            Type = Type + ".h5"

        for diag in diags:
            file = os.path.basename(diag.file)
            if file == Type:
                nDiags.append(diag)
            else:
                for key, regex in Diagnostic.h5FilesRegexes.items():
                    if Type == key and regex.match(file):
                        nDiags.append(diag)

        assert len(nDiags)

        return nDiags
