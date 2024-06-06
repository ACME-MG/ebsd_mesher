"""
 Title:         Mesher for 617 data
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from ebsd_mesher.interface import Interface

# Set up interface
itf = Interface()
itf.define_headers("x", "y", "grainId", "EulerMean_phi1", "EulerMean_Phi", "EulerMean_phi2")

# Read EBSD data
EBSD_FOLDER = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data/20240523 (ansto_617_s1)/hagb15/"
FILE_NAME = "ebsdExportColumnsTable_NoFill.csv"
itf.import_ebsd(f"{EBSD_FOLDER}/01_strain_0pct_on_stage_finalMapData20/{FILE_NAME}", 5)

# Plot the original EBSD map
itf.plot_ebsd(
    ebsd_path = "original",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Process the EBSD map
itf.decrease_resolution(20)
itf.remove_grains(2500)
itf.fill(15)
itf.clean(2)
itf.smooth(2)

# Plot the cleaned up EBSD map
itf.plot_ebsd(
    ebsd_path = "processed",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Mesh the EBSD map and plot
itf.mesh("~/cubit/psculpt.exe", z_length=200, z_voxels=3)
itf.plot_mesh(
    ipf      = "x",
    figure_x = 20
)
