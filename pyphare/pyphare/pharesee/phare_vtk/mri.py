#
#
#

import vtk

from . import base  # impor VtkTensorFieldFile, PhaseOutput


def plane_cutter(output, **kwargs):
    plane = vtk.vtkPlane()
    plane.SetNormal(0, 0, 1)

    cutter = vtk.vtkCutter()
    cutter.SetInputConnection(output.GetOutputPort())
    cutter.SetCutFunction(plane)
    return base.PhaseOutput(plane=plane, cutter=cutter)


def mri_phases():
    return [
        # surface_filter,
        base.plane_cutter,
        base.composite_data_geometry_filter,
        base.poly_data_mapper,
    ]


def scan(vtk_file, out_file="vtk.mp4", slices="all"):
    if isinstance(slices, str):
        if slices not in ["all"]:
            raise RuntimeError('Invalid slices param, must be int or str("all")')

    if isinstance(vtk_file, str):
        vtk_file = base.VtkTensorFieldFile(vtk_file, phases=mri_phases())

    vtk_file.geom.GetOutput().GetPointData().SetActiveScalars("data")

    colors = vtk.vtkNamedColors()
    renderer = vtk.vtkRenderer()
    renderer.AddActor(vtk.vtkActor(mapper=vtk_file.mapper))
    renderer.SetBackground(colors.GetColor3d("Silver"))  # .SetBackground(0, 0, 0)

    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetSize(800, 600)
    ren_win.OffScreenRenderingOn()  # do not flash image on screen

    ds = vtk_file.reader.GetOutput()
    # geom = vtk_file.geom
    bounds = vtk_file.bounds
    print(ds.GetClassName())
    print("Points:", ds.GetNumberOfPoints())
    print("Cells:", ds.GetNumberOfCells())
    print("Bounds:", bounds)

    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.SetParallelProjection(True)
    # camera.SetPosition(camera.GetPosition())
    # camera.SetFocalPoint(camera.GetFocalPoint())
    # camera.SetViewUp(camera.GetViewUp())

    # bounds = vtk_file.reader.GetOutput().GetBounds(0)
    zmin, zmax = bounds[4], bounds[5]
    dz = vtk_file.reader.GetOutput().GetDataSet(0, 0).GetSpacing()[2]
    n_slices = int((zmax - zmin) / dz) + 1
    print("n_slices", n_slices, zmax, zmin, bounds, dz)

    for i in range(n_slices):
        znew = zmin + i * dz
        print("znew", znew)
        vtk_file.plane.SetOrigin(0, 0, znew)

        # vtk_file.cutter.Modified()
        # vtk_file.cutter.Update()
        # vtk_file.geom.Update()

        # vtk_file.mapper.SetInputConnection(vtk_file.geom.GetOutputPort())
        # vtk_file.mapper.SelectColorArray("data")
        # vtk_file.mapper.SetScalarRange(
        #     vtk_file.geom.GetOutput().GetPointData().GetArray("data").GetRange()
        # )

        # renderer.ResetCameraClippingRange()
        ren_win.Render()

        # Save PNG
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(ren_win)
        w2if.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(f"slice_{i:04d}.png")
        writer.SetInputConnection(w2if.GetOutputPort())
        writer.Write()
