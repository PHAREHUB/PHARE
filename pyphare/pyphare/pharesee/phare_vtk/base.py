#
#
#

import vtk


class PhaseOutput:
    def __init__(self, **kwargs):
        self.kwargs = {**kwargs}

    def __iter__(self):
        return self.kwargs.items().__iter__()


@staticmethod
def surface_filter(output, **kwargs):
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputConnection(output.GetOutputPort())
    return PhaseOutput(surface=surface)


@staticmethod
def composite_data_geometry_filter(output, **kwargs):
    geom = vtk.vtkCompositeDataGeometryFilter()
    geom.SetInputConnection(output.GetOutputPort())
    return PhaseOutput(geom=geom)


def poly_data_mapper(output, **kwargs):
    array_name = kwargs.get("array_name", "data")
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(output.GetOutputPort())
    mapper.SelectColorArray(array_name)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(
        output.GetOutput().GetPointData().GetArray(array_name).GetRange(0)
    )
    return PhaseOutput(mapper=mapper)


def all_times_in(reader):
    info = reader.GetOutputInformation(0)

    time_key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()

    if not info.Has(time_key):
        raise RuntimeError("VTK Error: cannot get times from file")

    times = info.Get(time_key)
    return times


def bounds_in(reader):
    amr = reader.GetOutput()
    bbox = vtk.vtkBoundingBox()

    for level in range(amr.GetNumberOfLevels()):
        for idx in range(amr.GetNumberOfDataSets(level)):
            ds = amr.GetDataSet(level, idx)
            if ds is not None:
                b = [0.0] * 6
                ds.GetBounds(b)
                bbox.AddBounds(b)

    global_bounds = [0.0] * 6
    bbox.GetBounds(global_bounds)
    return global_bounds


class VtkFile:
    @staticmethod
    def _phases():
        return [surface_filter, composite_data_geometry_filter, poly_data_mapper]

    def __init__(self, filename, time=None, array_name="data", phases=_phases()):
        if len(phases) == 0:
            raise RuntimeError("Error: Zero phases!")

        self.array_name = array_name
        self.reader = vtk.vtkHDFReader(file_name=filename)
        self.reader.UpdateInformation()

        self.times = all_times_in(self.reader)
        self.reader.UpdateTimeStep(time if time else self.times[-1])

        _in = self.reader
        out = [None] * len(phases)
        for i in range(len(phases)):
            ret = phases[i](_in, array_name=array_name)
            for key, val in ret:
                if hasattr(val, "Update"):
                    val.Update()
                setattr(self, key, val)
                _in = val

        self.mapper = _in
        self.bounds = bounds_in(self.reader)
        self.spacing = self.reader.GetOutput().GetDataSet(0, 0).GetSpacing()

    def _count_type(self, endswith):
        return sum([1 for out in self.outputs if type(out).__name__.endswith(endswith)])


class VtkFieldFile(VtkFile):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class VtkTensorFieldFile(VtkFile):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
