# This script must be run with pvpython.
# It opens a pre-computed activation_time.vtu file and creates a high-quality plot.

from paraview.simple import *
import os
import sys

# --- Configuration ---
IMAGE_RESOLUTION = [1200, 1000]

# --- Input Validation ---
if len(sys.argv) < 3:
    print("Usage: pvpython generate_single_plot.py <path_to_vtu> <output_png>")
    sys.exit(1)
vtu_file_path = sys.argv[1]
output_png_path = sys.argv[2]

if not os.path.exists(vtu_file_path):
    print(f"Error: Input file not found at {vtu_file_path}")
    sys.exit(1)

# --- 1. Load the pre-computed data ---
reader = OpenDataFile(vtu_file_path)

# --- 2. Visualization Pipeline ---
view = CreateView('RenderView')
view.ViewSize = IMAGE_RESOLUTION
view.Background = [0.32, 0.34, 0.43]

display = Show(reader, view)
ColorBy(display, ('POINTS', 'activationTime'))
lut = GetColorTransferFunction('activationTime')
lut.ApplyPreset('Turbo', True)

# --- THE DEFINITIVE FIX: Call the function with no arguments ---
# This automatically rescales to the actual data range of the 'activationTime' array.
lut.RescaleTransferFunctionToDataRange()

# Set final representation properties for a clean look
display.Representation = 'Surface'
display.ColorArrayName = ['POINTS', 'activationTime']
display.LookupTable = lut
display.Interpolation = 'Gouraud' 

scalar_bar = GetScalarBar(lut, view)
scalar_bar.Title = "activation time"
scalar_bar.ComponentTitle = ''

# --- 3. Frame the Shot and Save ---
ResetCamera()
view.CameraParallelScale *= 0.9

SaveScreenshot(output_png_path, view, ImageResolution=IMAGE_RESOLUTION)
print(f"Successfully saved final plot to {output_png_path}")