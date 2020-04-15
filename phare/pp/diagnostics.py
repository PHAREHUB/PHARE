# formatted with black code formatter

import os, h5py, itertools
from .diagnostic.patch_data import (
    _EMPatchData,
    _FluidPatchData,
    _ParticlePatchData,
    Particles,
)
from .diagnostic.patch import Patch, PatchLevel


def _regex_to_ignore_population_name_in(h5FilesPerPatchDataType):
    """Validates h5 input file to ignore files which are not to be read/extracted"""
    import re

    pop = "${POP}"
    regexing = {Type: [] for Type in h5FilesPerPatchDataType.keys()}
    for Type, files in h5FilesPerPatchDataType.items():
        for file in files:
            """Currently only files with a population are handled here"""
            if pop not in file:
                continue
            bits = file.split(pop)
            left, right = bits[0], bits[1]

            """example append = "^ions_pop_ions_(.*)_domain.h5$"
                which is a regex string to ignore the population name in such a file
            """
            regexing[Type].append((file, re.compile("^" + left + "(.*)" + right + "$")))
    return regexing


"""Dict per type with list of possible input HDF5 Files/file quantity strings"""
h5FilesPerPatchDataType = {
    _EMPatchData: ["EM_B.h5", "EM_E.h5"],
    _FluidPatchData: [
        "ions_bulkVelocity.h5",
        "ions_density.h5",
        "ions_pop_ions_${POP}_density.h5",
        "ions_pop_ions_${POP}_flux.h5",
    ],
    _ParticlePatchData: [
        "ions_pop_ions_${POP}_domain.h5",
        "ions_pop_ions_${POP}_patchGhost.h5",
        "ions_pop_ions_${POP}_levelGhost.h5",
    ],
}
"""List of all files with no associated type"""
h5Files = [
    f for files in itertools.chain(h5FilesPerPatchDataType.values()) for f in files
]

"""PatchDataType to tuple (quantity_file_string, regex_for_quantity_file_string"""
h5FilesRegexesPerPatchDataType = _regex_to_ignore_population_name_in(
    h5FilesPerPatchDataType
)
"""File quantity string to compiled regex for said quauntity string"""
h5FilesRegexesPerQuantity = {
    PatchDataType: regex_quantity_list
    for rgxs in itertools.chain(h5FilesRegexesPerPatchDataType.values())
    for PatchDataType, regex_quantity_list in rgxs
}


class Diagnostic:
    ROUNDING_PLACES = 5
    """
    Convenience class for accessing quantities within contained Diagnostic objects

    Parameters:
    -------------------

    file : str of HDF5 file to use for extracting contained quantity data

    """

    def __init__(self, file):
        self.file = file
        self.patch_data_type = patchDataTypeFrom(file)
        self.h5file = h5py.File(file, "r")
        self.levels = {}
        self.dl = [float(v) for v in self.h5file.attrs["meshSize"].split(",")]
        self.cells = [int(v) for v in self.h5file.attrs["cells"].split(",")]
        self.origin = [float(v) for v in self.h5file.attrs["origin"].split(",")]
        self.domain_upper = [
            round(dl * n + ori, Diagnostic.ROUNDING_PLACES)
            for dl, n, ori in zip(self.dl, self.cells, self.origin)
        ]
        self.dim = len(self.dl)
        self.initialize_patch_hierarchy_from_hdf5_file()

    def initialize_patch_hierarchy_from_hdf5_file(self):
        for timestamp in self.h5file.keys():
            for patch_level in self.h5file[timestamp].keys():
                lvl = int(patch_level[2:])  # pl123 - get 123
                self.levels[lvl] = PatchLevel(self, lvl)
                patches = []
                for patch in self.h5file[timestamp][patch_level].keys():
                    patch_data = _extract_patch_data(
                        self, self.h5file[timestamp][patch_level][patch]
                    )
                    # levelGhost for level0 non-border patches can be None
                    if patch_data is None:
                        continue
                    patches.append(
                        Patch(
                            self.levels[lvl],
                            self.h5file[timestamp][patch_level][patch],
                            _extract_patch_data(
                                self, self.h5file[timestamp][patch_level][patch]
                            ),
                        )
                    )
                self.levels[lvl].patches = {p.id: p for p in patches}


class Diagnostics:
    """
    Convenience class for accessing quantities within contained Diagnostic objects

    Parameters:
    -------------------

    input        : str of directory ot extract Diagnostics, or list[Diagnostic]

    """

    def __init__(self, input):

        self.diags = input
        if isinstance(input, str):
            self.diags = extract_diagnostics(input)

        assert len(self.diags) and all(
            [isinstance(diag, Diagnostic) for diag in self._diags_as_list()]
        )

    def ionDensity(self):
        return self._selectDiag("ions_density")

    def ionVelocity(self):
        return self._selectDiag("ions_bulkVelocity")

    def popDensities(self):
        return self._selectDiags("ions_pop_ions_${POP}_density")

    def popFluxes(self):
        return self._selectDiags("ions_pop_ions_${POP}_flux")

    def E(self):
        return self._selectDiag("EM_E")

    def B(self):
        return self._selectDiag("EM_B")

    def _selectDiags(self, file):
        return selectDiags(self._diags_as_list(), file)

    def _selectDiag(self, file):
        return selectDiag(self._diags_as_list(), file)

    def _diags_as_list(self):
        """full list of diagsnostics with no associated PatchDataType"""
        return [
            diag for diags in itertools.chain(self.diags.values()) for diag in diags
        ]


def extract_diagnostics(diag_dir):
    assert isinstance(diag_dir, str)

    diags = [Diagnostic(file) for file in _h5_files_in(diag_dir)]
    if len(diags) == 0:
        raise ValueError(
            "Cannot load diagnostics from specific simulation: no diagnostics found"
        )
    diagDic = {
        clazz.__name__: []
        for clazz in [_EMPatchData, _FluidPatchData, _ParticlePatchData]
    }
    for diag in diags:
        diagDic[diag.patch_data_type.__name__].append(diag)
    return diagDic


def patchDataTypeFrom(file):
    assert isinstance(file, str)

    file = os.path.basename(file)
    for PatchDataType, files in h5FilesPerPatchDataType.items():
        if file in files:
            return PatchDataType
    for PatchDataType, regexs in h5FilesRegexesPerPatchDataType.items():
        for h5File, regex in regexs:
            if regex.match(file):
                return PatchDataType
    raise ValueError("Unable to deduce PatchDataType from file " + file)


def selectDiags(diags, quantity):
    selectedDiags = []

    if not quantity.endswith(".h5"):
        quantity = quantity + ".h5"

    for diag in diags:
        file = os.path.basename(diag.file)
        if file == quantity:
            selectedDiags.append(diag)
        else:
            for key, regex in h5FilesRegexesPerQuantity.items():
                if quantity == key and regex.match(file):
                    selectedDiags.append(diag)

    assert len(selectedDiags)

    return selectedDiags


def selectDiag(diags, quantity):
    """When you know only one diagnostic should exist for specified quantity"""
    diags = selectDiags(diags, quantity)
    assert len(diags) == 1
    return diags[0]


def _extract_patch_data_recursive(diag, h5item, pre_key=""):
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
        return _extract_patch_data_recursive(diag, h5item[kinder[0]], kinder[0])


def _extract_patch_data(diag, h5item):
    assert isinstance(diag, Diagnostic)
    return _extract_patch_data_recursive(diag, h5item)


def _is_input_file_valid(file):
    return file in h5Files or any(
        [regex.match(file) for regex in h5FilesRegexesPerQuantity.values()]
    )


def _h5_files_in(directory: str):

    h5files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if _is_input_file_valid(file):
                h5files.append(os.path.join(root, file))
            else:
                print("Warning: Input file not recognized, ignored: ", file)
    return h5files
