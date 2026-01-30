#
#
#


import vtk
from .base import VtkTensorFieldFile


def plot(vtk_file, out_file="vtk.png"):
    if isinstance(vtk_file, str):
        vtk_file = VtkTensorFieldFile(vtk_file)
    vtk_file.geom.GetOutput().GetPointData().SetActiveScalars("data")

    colors = vtk.vtkNamedColors()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(vtk.vtkActor(mapper=vtk_file.mapper))
    renderer.SetBackground(colors.GetColor3d("Silver"))

    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetSize(800, 600)
    ren_win.OffScreenRenderingOn()  # do not flash image on screen

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    ren_win.Render()

    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(ren_win)
    w2if.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(out_file)
    writer.SetInputConnection(w2if.GetOutputPort())
    writer.Write()
