from pathlib import Path

from paraview.simple import ColorBy, GetActiveViewOrCreate, LegacyVTKReader, Render, Show


vtk_folder = Path(__file__).parent / "build" / "output" / "vtk"
vtk_files = sorted(vtk_folder.glob("exercise1_*.vtk"))

if not vtk_files:
    raise RuntimeError("Run learnSPH_app first.")

reader = LegacyVTKReader(FileNames=[str(path) for path in vtk_files])
view = GetActiveViewOrCreate("RenderView")
display = Show(reader, view)
display.Representation = "Points"
display.PointSize = 8
ColorBy(display, ("POINTS", "scalar"))
display.RescaleTransferFunctionToDataRange(True, False)
view.ResetCamera()
Render(view)
