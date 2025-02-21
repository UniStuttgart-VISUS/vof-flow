import os
from paraview.simple import *

# Config
dataset = '~/vofflow/data/jet-collision-ds2/jet-collision-ds2.pvd'
if os.name == 'nt':
    plugin_path = '~/vofflow/install/bin/paraview-5.13/plugins/VofFlow/VofFlow.dll'
else:
    plugin_path = '~/vofflow/install/lib/paraview-5.13/plugins/VofFlow/VofFlow.so'
output_path = '~/vofflow/output/result'

init_timestep = 68  # The time steps used in the paper (fig. 14) are: 68, 82, 96, 110, 124, 138, 152.
ghost_cells = 10  # MPI ghost cells, use 36 for full dataset, 20 for "ds1", and 10 for "ds2".

# Load plugin
LoadPlugin(os.path.expanduser(plugin_path), ns=globals())

# Read dataset
reader = PVDReader(registrationName='dataset', FileName=os.path.expanduser(dataset))
reader.CellArrays = ['vof-function[-]', 'velocity[cm/s]', 'f3-function[-]', 'n_c_3ph[1]']
time_steps = reader.TimestepValues

# Configure VofTracking
vofTracking1 = VofTracking(registrationName='VofTracking1', InputGrid=reader)
vofTracking1.UseThreePhase = 1
vofTracking1.VoF = ['CELLS', 'vof-function[-]']
vofTracking1.VoF3 = ['CELLS', 'f3-function[-]']
vofTracking1.VoFNorm = ['CELLS', 'n_c_3ph[1]']
vofTracking1.Velocity = ['CELLS', 'velocity[cm/s]']
vofTracking1.UseComponents = 0
vofTracking1.ComponentsVoF = ['CELLS', 'None']
vofTracking1.ComponentsVoF3 = ['CELLS', 'None']
vofTracking1.UseTargetTimeStep = 1
vofTracking1.InitTimeStep = init_timestep
vofTracking1.TargetTimeStep = 152
vofTracking1.Refinement = 2
vofTracking1.NeighborCorrection = 1
vofTracking1.CellCorrection = 1
vofTracking1.PLICCorrection = 1
vofTracking1.IntegrationMethod = 'RK4'  # 'Euler'
vofTracking1.IntegrationSubSteps = 8
vofTracking1.Epsilon = 1e-05
vofTracking1.NumIterations = 15
vofTracking1.GhostCells = ghost_cells
vofTracking1.CutLabels = 0
vofTracking1.LabelCutType = 'Plane'
vofTracking1.BoundaryMethod = 'DiscreteFlyingEdges3D'  # 'DiscreteMarchingCubes'
vofTracking1.OutputDataType = 'vtkImageData'  # 'vtkRectilinearGrid'
vofTracking1.OutputState = 1
vofTracking1.OutputTimeMeasure = 1
vofTracking1.MirrorXMin = 0
vofTracking1.MirrorXMax = 0
vofTracking1.MirrorYMin = 0
vofTracking1.MirrorYMax = 0
vofTracking1.MirrorZMin = 1
vofTracking1.MirrorZMax = 0

# Check output dir
output_path = os.path.expanduser(output_path)
os.makedirs(output_path, exist_ok=True)

# Save output
SaveData(os.path.join(output_path, 'grid.pvd'), proxy=OutputPort(vofTracking1, 0), CompressorType='ZLib', CompressionLevel='1')
SaveData(os.path.join(output_path, 'seeds.pvd'), proxy=OutputPort(vofTracking1, 1), CompressorType='ZLib', CompressionLevel='1')
SaveData(os.path.join(output_path, 'part.pvd'), proxy=OutputPort(vofTracking1, 2), CompressorType='ZLib', CompressionLevel='1')
SaveData(os.path.join(output_path, 'bound.pvd'), proxy=OutputPort(vofTracking1, 3), CompressorType='ZLib', CompressionLevel='1')
