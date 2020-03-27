# formatted with black code formatter

import os, h5py, phare.pharein as ph, itertools
from .diagnostic.patch_data import _EM, _Fluid, _Particle, Particles
from .diagnostic.patch import Patch, PatchLevel


class H5FileRegexer:
    """Used to Validates h5 input file to ignore files which are not to be read/extracted"""

    @staticmethod
    def ignore_population_name_in(h5FilesPerType):
        import re

        pop = "${POP}"
        regexing = {Type: [] for Type in h5FilesPerType.keys()}
        for Type, files in h5FilesPerType.items():
            for file in files:
                if pop not in file:
                    continue
                bits = file.split(pop)
                left, right = bits[0], bits[1]
                regexing[Type].append(
                    (file, re.compile("^" + left + "(.*)" + right + "$"))
                )
        return regexing


class Diagnostic:
    REFINEMENT_RATIO = 2
    """12 decimal places is arbitrary, but hopefully sufficient"""
    DL_PRECISION = 12

    """List of all possible input HDF5 Files"""
    h5FilesPerType = {
        _EM: ["EM_B.h5", "EM_E.h5"],
        _Fluid: [
            "ions_bulkVelocity.h5",
            "ions_density.h5",
            "ions_pop_ions_${POP}_density.h5",
            "ions_pop_ions_${POP}_flux.h5",
        ],
        _Particle: [
            "ions_pop_ions_${POP}_domain.h5",
            "ions_pop_ions_${POP}_patchGhost.h5",
            "ions_pop_ions_${POP}_levelGhost.h5",
        ],
    }
    h5Files = [f for files in itertools.chain(h5FilesPerType.values()) for f in files]

    """Compiled Regexs to ignore population name in HDF5 input files"""
    h5FilesRegexesPerType = H5FileRegexer.ignore_population_name_in(h5FilesPerType)
    h5FilesRegexes = {
        r[0]: r[1]
        for rgxs in itertools.chain(h5FilesRegexesPerType.values())
        for r in rgxs
    }

    def __init__(self, sim: ph.Simulation, file):
        def truncate(number, digits) -> float:
            import math

            stepper = 10.0 ** digits
            return math.trunc(stepper * number) / stepper

        self.sim = sim
        self.file = file
        self.patch_data_type = Diagnostics.get_type_from_file(file)
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
    def __init__(self, input):
        if isinstance(input, ph.Simulation):
            self.diags = Diagnostics.extract(input)
        self.diags = input

    def getPopDensities(self):
        return Diagnostics.get(self.diags, "ions_pop_ions_${POP}_density")

    def getDensity(self, pop=""):
        if len(pop):
            return Diagnostics.get(self.diags, "ions_pop_ions_" + pop + "_density")[0]
        return Diagnostics.get(self.diags, "ions_density")[0]

    @staticmethod
    def _extract_recursive(diag, h5item, pre_key=""):
        assert isinstance(diag, Diagnostic)

        kinder = list(h5item.keys())
        n_keys = len(kinder)

        if n_keys == 1 and isinstance(h5item[kinder[0]], h5py.Dataset):
            """Single dataset per patch hdf5 file, e.g. ion density"""
            return diag.patch_data_type(diag, kinder[0], h5item[kinder[0]])

        if n_keys > 1 and isinstance(h5item[kinder[0]], h5py.Dataset):
            """Multiple dataset per patch hdf5 file, e.g. ion flux"""
            return diag.patch_data_type(diag, pre_key, h5item)

        if n_keys:
            return Diagnostics._extract_recursive(diag, h5item[kinder[0]], kinder[0])

    @staticmethod
    def _extract(diag, h5item):
        assert isinstance(diag, Diagnostic)
        return Diagnostics._extract_recursive(diag, h5item)

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
    def get_type_from_file(file):
        file = os.path.basename(file)
        for Type, files in Diagnostic.h5FilesPerType.items():
            if file in files:
                return Type
        for Type, regexs in Diagnostic.h5FilesRegexesPerType.items():
            for h5File, regex in regexs:
                if regex.match(file):
                    return Type
        assert False

    @staticmethod
    def get(diags, Type):
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
