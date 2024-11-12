import os
from paraview.simple import *

# Config
dataset = '~/vofflow/data/jet-collision-ds2/jet-collision-ds2.0152.vtr'
plugin_path = '~/vofflow/install/lib/paraview-5.12/plugins/VofFlow/VofFlow.so'
output_path = '~/vofflow/output/jet_plic.png'

# Load plugin
LoadPlugin(os.path.expanduser(plugin_path), ns=globals())

LoadPalette(paletteName='WhiteBackground')

renderView1 = GetActiveViewOrCreate('RenderView')

# Read dataset
reader = XMLRectilinearGridReader(registrationName='dataset', FileName=os.path.expanduser(dataset))
reader.CellArrayStatus = ['f3-function[-]', 'n_c_3ph[1]', 'vof-function[-]']

# PLIC3 Filter
pLIC31 = PLIC3(registrationName='PLIC31', InputGrid=reader)
pLIC31.VoF3 = ['CELLS', 'f3-function[-]']
pLIC31.VoF = ['CELLS', 'vof-function[-]']
pLIC31.VoFNorm = ['CELLS', 'n_c_3ph[1]']
pLIC31.Epsilon = 1e-05
pLIC31.NumIterations = 15

# Filter for surface of inner phase
threshold1 = Threshold(registrationName='Threshold1', Input=pLIC31)
threshold1.Scalars = ['CELLS', 'PhaseType']
threshold1.LowerThreshold = 1.0
threshold1.UpperThreshold = 1.0

# Filter for surface of outer phase
threshold2 = Threshold(registrationName='Threshold2', Input=pLIC31)
threshold2.Scalars = ['CELLS', 'PhaseType']
threshold2.LowerThreshold = 2.0
threshold2.UpperThreshold = 2.0

# Mirror the inner surface
reflect1 = Reflect(registrationName='Reflect1', Input=threshold1)
reflect1.Plane = 'Z Min'
reflect1.CopyInput = 0
reflect1.FlipAllInputArrays = 0

# Mirror the outer surface
reflect2 = Reflect(registrationName='Reflect2', Input=threshold2)
reflect2.Plane = 'Z Min'
reflect2.CopyInput = 0
reflect2.FlipAllInputArrays = 0

# Surface colors
color_inner = [0.0, 1.0, 0.0]
color_outer = [0.8784313725490196, 0.8784313725490196, 0.8784313725490196]

# Display inner phase
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')
threshold1Display.Representation = 'Surface'
threshold1Display.AmbientColor = color_inner
threshold1Display.ColorArrayName = [None, '']
threshold1Display.DiffuseColor = color_inner
threshold1Display.Specular = 0.5
threshold1Display.SpecularPower = 25.0

# Display outer phase
threshold2Display = Show(threshold2, renderView1, 'UnstructuredGridRepresentation')
threshold2Display.Representation = 'Surface'
threshold2Display.AmbientColor = color_outer
threshold2Display.ColorArrayName = [None, '']
threshold2Display.DiffuseColor = color_outer
threshold2Display.Opacity = 0.5
threshold2Display.Specular = 0.5
threshold2Display.SpecularPower = 25.0

# Display mirrored inner phase
reflect1Display = Show(reflect1, renderView1, 'UnstructuredGridRepresentation')
reflect1Display.Representation = 'Surface'
reflect1Display.AmbientColor = color_inner
reflect1Display.ColorArrayName = [None, '']
reflect1Display.DiffuseColor = color_inner
reflect1Display.Specular = 0.5
reflect1Display.SpecularPower = 25.0

# Display mirrored outer phase
reflect2Display = Show(reflect2, renderView1, 'UnstructuredGridRepresentation')
reflect2Display.Representation = 'Surface'
reflect2Display.AmbientColor = color_outer
reflect2Display.ColorArrayName = [None, '']
reflect2Display.DiffuseColor = color_outer
reflect2Display.Opacity = 0.5
reflect2Display.Specular = 0.5
reflect2Display.SpecularPower = 25.0

# Camera setup
renderView1.ResetActiveCameraToPositiveZ()
renderView1.ResetCamera(False)
renderView1.CameraPosition = [0.34060654044151306, -0.03483813628554344, 0.4401733070104946]
renderView1.CameraFocalPoint = [0.34060654044151306, -0.03483813628554344, 0.005889910257715658]
renderView1.CameraViewUp = [0.0, 1.0, 0.0]
renderView1.OrientationAxesVisibility = 0
renderView1.Update()

# Save screenshot
output_path = os.path.expanduser(output_path)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
SaveScreenshot(output_path, renderView1, ImageResolution=[2400, 800])
