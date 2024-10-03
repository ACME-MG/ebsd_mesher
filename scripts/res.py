"""
 Title:         Mesher for 617 data
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from ebsd_mesher.interface import Interface

# Set up interface
itf = Interface()
itf.define_headers("x", "y", "grainId", "Euler_phi1", "Euler_Phi", "Euler_phi2")

# Constants
GRIP_SIZE = 7

# Read EBSD data
# itf.import_ebsd("/mnt/c/Users/z5208868/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/prior_with_stage/res5gs20/ebsdExportColumnsTableReduced_FillRegion.csv", 5)
# itf.import_ebsd("/mnt/c/Users/z5208868/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/prior_with_stage/res10gs10/ebsdExportColumnsTableReduced_FillRegion.csv", 10)
itf.import_ebsd("/mnt/c/Users/z5208868/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/prior_with_stage/res15gs10/ebsdExportColumnsTableReduced_FillRegion.csv", 15)
# itf.import_ebsd("/mnt/c/Users/z5208868/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/prior_with_stage/res20gs5/ebsdExportColumnsTableReduced_FillRegion.csv",  20)

# Plot raw EBSD map
itf.plot_ebsd(
    ebsd_path = "1_ebsd_original",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Remove voided regions
bounds = itf.get_bounds()
itf.redefine_domain(
    x_min = bounds["x_min"],
    x_max = bounds["x_max"],
    y_min = bounds["y_min"],
    y_max = bounds["y_max"]
)

# Set dimensions to nominal
x_size = bounds["x_max"]-bounds["x_min"]
y_size = bounds["y_max"]-bounds["y_min"]
nominal_length = 2200
nominal_width = 1560
itf.redefine_domain(
    x_min = (x_size-nominal_length)/2,
    x_max = (x_size+nominal_length)/2,
    y_min = (y_size-nominal_width)/2,
    y_max = (y_size+nominal_width)/2,
)
itf.export_bounds()

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "2_ebsd_rebound",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Process the EBSD map
itf.remove_grains(1000)
itf.decrease_resolution(3)
itf.fill(5)
# itf.smooth(1)
itf.add_grips(GRIP_SIZE)

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "3_ebsd_processed",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Mesh the EBSD map and plot
itf.mesh("~/cubit/psculpt.exe", z_elements=1)
itf.export_grains(degrees=False)
itf.export_elements(degrees=False)
itf.plot_mesh(
    mesh_path = "4_mesh_raw",
    ipf       = "x",
    figure_x  = 40
)

# Fix the interfaces of the grips and scale the mesh
itf.fix_grip_interfaces(nominal_length/3, nominal_length)
itf.scale_mesh(nominal_width, "y")
itf.scale_mesh(20, "z")
itf.plot_mesh(
    mesh_path = "5_mesh_fixed",
    ipf       = "x",
    figure_x  = 40
)
