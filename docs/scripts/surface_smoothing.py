import os
from paraview.simple import *

# Config
dataset = '~/vofflow/output/result/seeds.pvd'
if os.name == 'nt':
    plugin_path = '~/vofflow/install/bin/paraview-5.13/plugins/VofFlow/VofFlow.dll'
else:
    plugin_path = '~/vofflow/install/lib/paraview-5.13/plugins/VofFlow/VofFlow.so'
output_path = '~/vofflow/output/result'

# Load plugin
LoadPlugin(os.path.expanduser(plugin_path), ns=globals())

# Read seed points
seedspvd = PVDReader(registrationName='seeds.pvd', FileName=os.path.expanduser(dataset))

# Generate boundary
vofBoundary1 = VofBoundary(registrationName='VofBoundary1', InputData=seedspvd)
vofBoundary1.BoundaryMode = 'MappedNormals'

# Check output dir
output_path = os.path.expanduser(output_path)
os.makedirs(output_path, exist_ok=True)

# Save output
SaveData(os.path.join(output_path, 'bound_smooth.pvd'), proxy=vofBoundary1, CompressorType='ZLib', CompressionLevel='1')
