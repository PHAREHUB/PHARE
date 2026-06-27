import sys
import vtk
import h5py
import logging
import matplotlib.pyplot as plt

from pyphare.logger import getLogger

logging.basicConfig(level=logging.INFO)

logger = getLogger(__name__)
data_path = "/home/p/git/phare/master/phare_outputs/harris_3d/EM_B.vtkhdf"
# data_path = "phare_outputs/vtk_diagnostic_test/test_vtk/test_dump_diags_1/1/2/1/ions_pop_alpha_flux.vtkhdf"
array_name = "data"
fps = 5

if __name__ == "__main__":
    # from pyphare.pharesee.phare_vtk import show as vtk_show
    # from pyphare.pharesee.phare_vtk import plot as vtk_plot
    # from pyphare.pharesee.phare_vtk import scan as vtk_scan

    # vtk_plot(data_path)
    # # vtk_show(data_path)
    # vtk_scan(data_path)

    component = None  # 0/1/2 for x/y/z, or None for magnitude
    output_pattern = "slice_{:04d}.png"
    image_size = (1024, 1024)

    # -----------------------------
    # 1. Reader
    # -----------------------------
    reader = vtk.vtkHDFReader(file_name=data_path)
    reader.UpdateInformation()
    reader.UpdateTimeStep(1.0)
    print("time", reader.GetTimeValue())

    # -----------------------------
    # 4. Plane and cutter
    # -----------------------------
    plane = vtk.vtkPlane()
    plane.SetNormal(0, 0, 1)

    cutter = vtk.vtkCutter()
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.SetCutFunction(plane)

    geom = vtk.vtkCompositeDataGeometryFilter()
    geom.SetInputConnection(cutter.GetOutputPort())

    geom_all = vtk.vtkCompositeDataGeometryFilter()
    geom_all.SetInputConnection(reader.GetOutputPort())
    geom_all.Update()
    poly_all = geom_all.GetOutput()

    # -----------------------------
    # 5. Vector magnitude / component
    # -----------------------------
    if component is None:
        # magnitude
        calc = vtk.vtkArrayCalculator()
        calc.SetInputConnection(geom.GetOutputPort())
        calc.AddVectorArrayName(array_name)
        calc.SetFunction(f"mag({array_name})")
        calc.SetResultArrayName("mag")
        calc.Update()
        scalar_name = "mag"
        mapper_input = calc.GetOutputPort()
        assert calc.GetOutput().GetPointData().GetArray("mag")
    else:
        # single component
        geom.Update()
        # pd = geom.GetOutput().GetPointData()
        scalar_name = array_name
        mapper_input = geom.GetOutputPort()

    # -----------------------------
    # 6. Mapper and actor
    # -----------------------------
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(mapper_input)
    mapper.SetScalarModeToUsePointData()
    mapper.SelectColorArray(scalar_name)

    # Set fixed scalar range
    if component is None:
        arr = calc.GetOutput().GetPointData().GetArray("mag")
        assert arr
    else:
        arr = geom.GetOutput().GetPointData().GetArray(array_name)
        assert arr
    mapper.SetScalarRange(arr.GetRange(component or 0))

    amr = reader.GetOutput()

    # -----------------------------
    # 2. Compute global bounds across all AMR levels
    # -----------------------------
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
    xmin, xmax, ymin, ymax, zmin, zmax = global_bounds

    # -----------------------------
    # 3. Determine slice spacing
    # -----------------------------
    # Use level 0 spacing for dz
    ds0 = amr.GetDataSet(0, 0)
    dz = ds0.GetSpacing()[2]
    n_slices = int((zmax - zmin) / dz) + 1
    z_values = [zmin + i * dz for i in range(n_slices)]

    print(
        f"Generating {n_slices} slices from z={zmin:.3f} to z={zmax:.3f} with dz={dz:.3f}"
    )

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetRepresentationToSurface()

    # -----------------------------
    # 7. Renderer and offscreen window
    # -----------------------------
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.2, 0.2, 0.2)

    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetSize(*image_size)
    ren_win.OffScreenRenderingOn()

    # Lock camera
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.SetParallelProjection(True)

    # -----------------------------
    # 8. Slice loop
    # -----------------------------
    points = poly_all.GetPoints()
    z_coords = [points.GetPoint(i)[2] for i in range(points.GetNumberOfPoints())]
    z_values = sorted(set(z_coords))

    for i in range(n_slices):
        z = z_values[i]
        plane.SetOrigin(0, 0, z)
        cutter.Modified()  # mark cutter dirty
        geom.Update()  # propagate
        if component is None:
            calc.Update()  # recompute magnitude

        ren_win.Render()

        # Capture PNG
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(ren_win)
        w2if.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(w2if.GetOutputPort())
        writer.SetFileName(output_pattern.format(i))
        writer.Write()

        print(
            f"Saved slice {i+1}/{n_slices} at z={z:.3f}, cells={cutter.GetOutputDataObject(0).GetNumberOfCells()}",
            z,
        )

    import shlex
    import subprocess

    cmd = shlex.split(
        f"ffmpeg -y -framerate {fps} -i slice_%04d.png -c:v libx264 -pix_fmt yuv420p mri.mp4"
    )
    subprocess.call(cmd)

else:
    colors = vtk.vtkNamedColors()
    reader = vtk.vtkHDFReader(file_name=data_path)
    # reader.SetInputArrayToProcess(
    #     0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "data"
    # )

    # f0 = vtk.vtkOutlineFilter()
    # f0.SetInputConnection(reader.GetOutputPort())
    # f0.Update()

    # plane = vtk.vtkPlane()
    # plane.SetNormal(0, 0, 1)
    # plane.SetOrigin(0, 0, 0)

    # cutter = vtk.vtkCutter()
    # cutter.SetCutFunction(plane)
    # cutter.SetInputConnection(reader.GetOutputPort())

    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputConnection(reader.GetOutputPort())
    surface.Update()

    print("surface", surface)
    print("surface.GetOutput", surface.GetOutput())
    # poly = surface.GetOutput()

    # cell_data = poly.GetCellData()
    # array_name = cell_data.GetArrayName(0)  # or explicit name
    # cell_data.SetActiveScalars(array_name)

    # calc = vtk.vtkArrayCalculator()
    # calc.SetInputConnection(surface.GetOutputPort())
    # calc.AddVectorArrayName(array_name)
    # calc.SetFunction(f"mag({array_name})")
    # calc.SetResultArrayName("magnitude")
    # calc.Update()

    geom = vtk.vtkCompositeDataGeometryFilter()
    geom.SetInputConnection(surface.GetOutputPort())
    geom.Update()

    poly = geom.GetOutput()
    print(poly.GetNumberOfCells(), poly.GetNumberOfPoints())

    pd = poly.GetPointData()
    for i in range(pd.GetNumberOfArrays()):
        print(pd.GetArrayName(i))

    print(poly.GetPointData().GetArray(array_name).GetRange())

    pd.SetActiveScalars("data")

    vec_array = pd.GetArray(array_name)
    print("Number of components:", vec_array.GetNumberOfComponents())

    # calc = vtk.vtkArrayCalculator()
    # calc.SetInputConnection(geom.GetOutputPort())
    # calc.AddVectorArrayName(array_name)
    # calc.SetFunction(f"mag({array_name})")
    # calc.SetResultArrayName("mag")
    # calc.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geom.GetOutputPort())
    mapper.SelectColorArray("data")
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(geom.GetOutput().GetPointData().GetArray("data").GetRange(0))
    # mapper.SetScalarRange(calc.GetOutput().GetPointData().GetArray("mag").GetRange())

    actor = vtk.vtkActor(mapper=mapper)
    # actor.GetProperty().SetDiffuse(0)
    # actor.GetProperty().SetAmbient(1)
    # actor.GetProperty().SetInterpolationToFlat()
    # actor.GetProperty().SetEdgeVisibility(True)
    actor.GetProperty().SetRepresentationToSurface()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d("Silver"))  # .SetBackground(0, 0, 0)

    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetSize(800, 600)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    ren_win.Render()
    # iren.Start()

    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(ren_win)
    w2if.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName("phare_surface.png")
    writer.SetInputConnection(w2if.GetOutputPort())
    writer.Write()
