#
#
#

from .base import *


def show(vtk_file):
    if isinstance(vtk_file, str):
        vtk_file = VtkTensorFieldFile(vtk_file)
    vtk_file.geom.GetOutput().GetPointData().SetActiveScalars("data")

    colors = vtk.vtkNamedColors()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(vtk.vtkActor(mapper=vtk_file.mapper))
    renderer.SetBackground(colors.GetColor3d("Silver"))  # .SetBackground(0, 0, 0)

    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetSize(800, 600)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    ren_win.Render()
    iren.Start()
