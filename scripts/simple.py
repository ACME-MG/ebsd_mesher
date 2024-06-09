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

# Read EBSD map
DATA_FOLDER = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data"
itf.import_ebsd(f"{DATA_FOLDER}/2024-06-07 (ansto_617_s1_simple)/simple_converted.csv", 20)
itf.add_grips(20)

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "original",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Mesh the EBSD map and plot
itf.mesh("~/cubit/psculpt.exe", z_voxels=3)
itf.plot_mesh(
    ipf      = "x",
    figure_x = 20
)

# Fix the interfaces of the grips and scale the mesh
itf.fix_grip_interfaces(2300/4, 2300)
itf.scale_mesh(1571, "y")
itf.scale_mesh(323, "z")
itf.plot_mesh(
    ipf      = "x",
    figure_x = 20
)
