import os
from paraview.simple import *

# Config
dataset = '~/vofflow/output/result/bound_gap_smoothed.vtp'
plugin_path = '~/vofflow/install/lib/paraview-5.12/plugins/VofFlow/VofFlow.so'
output_path = '~/vofflow/output/fig_jet_bound.png'

label_jet = 15  # Use 20 for full dataset, 20 for "ds1", and 15 for "ds2"
label_drop = 8  # Use 12 for full dataset, 12 for "ds1", and 8 for "ds2"

# Load plugin
LoadPlugin(os.path.expanduser(plugin_path), ns=globals())

LoadPalette(paletteName='WhiteBackground')

renderView1 = GetActiveViewOrCreate('RenderView')

# Read smoothed bound data
reader = XMLPolyDataReader(registrationName='bounds', FileName=[os.path.expanduser(dataset)])
reader.PointArrayStatus = ['Labels', 'Normals']

# Cutout jet label
threshold1 = Threshold(registrationName='Threshold1', Input=reader)
threshold1.Scalars = ['POINTS', 'Labels']
threshold1.LowerThreshold = label_jet
threshold1.UpperThreshold = label_jet

# Mirror for esthetic reasons, that jets flows from left to right in image
reflect1 = Reflect(registrationName='Reflect1', Input=threshold1)
reflect1.Plane = 'Z Min'
reflect1.CopyInput = 0
reflect1.FlipAllInputArrays = 0

# Cut bounds in mirror plane to view "into" them
clip1 = Clip(registrationName='Clip1', Input=reflect1)
clip1.ClipType = 'Plane'
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# Cutout drop label
threshold2 = Threshold(registrationName='Threshold2', Input=reader)
threshold2.Scalars = ['POINTS', 'Labels']
threshold2.LowerThreshold = label_drop
threshold2.UpperThreshold = label_drop

# Mirror for esthetic reasons, that jets flows from left to right in image
reflect2 = Reflect(registrationName='Reflect2', Input=threshold2)
reflect2.Plane = 'Z Min'
reflect2.CopyInput = 0
reflect2.FlipAllInputArrays = 0

# Cut bounds in mirror plane to view "into" them
clip2 = Clip(registrationName='Clip2', Input=reflect2)
clip2.ClipType = 'Plane'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# Render
color_jet = [1.0, 0.4980392156862745, 0.054901960784313725]
color_drop = [0.12156862745098039, 0.4666666666666667, 0.7058823529411765]

clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.AmbientColor = color_jet
clip1Display.DiffuseColor = color_jet
clip1Display.Specular = 0.5
clip1Display.SpecularPower = 25.0
clip1Display.SelectNormalArray = 'Normals'

clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.AmbientColor = color_drop
clip2Display.DiffuseColor = color_drop
clip2Display.Specular = 0.5
clip2Display.SpecularPower = 25.0
clip2Display.SelectNormalArray = 'Normals'

# Camera setup
renderView1.ViewSize = [2400, 800]
renderView1.ResetActiveCameraToPositiveZ()
renderView1.ResetCamera(False)
renderView1.CameraPosition = [0.29, -0.024, 0.4]
renderView1.CameraFocalPoint = [0.29, -0.024, 0.0]
renderView1.CameraViewUp = [0.0, 1.0, 0.0]
renderView1.OrientationAxesVisibility = 0
renderView1.Update()

# Save screenshot
output_path = os.path.expanduser(output_path)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
SaveScreenshot(output_path, renderView1, ImageResolution=[2400, 800])
