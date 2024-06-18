"""
 Title:         Mesher for 617 data
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from ebsd_mesher.interface import Interface

# Set up interface
itf = Interface(f"617_s1")
itf.define_headers("x", "y", "grainId", "Euler_phi1", "Euler_Phi", "Euler_phi2")

# Read EBSD data
# EBSD_FOLDER = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data/2024-05-23 (ansto_617_s1)/hagb15/"
EBSD_FOLDER = "/mnt/c/Users/z5208868/OneDrive - UNSW/PhD/data/2024-05-23 (ansto_617_s1)/hagb15/"
FILE_PATH = "01_strain_0pct_on_stage_finalMapData20/ebsdExportColumnsTable_NoFill.csv"
itf.import_ebsd(
    ebsd_path = f"{EBSD_FOLDER}/{FILE_PATH}",
    step_size = 5,
    degrees   = True
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
nominal_length = 2170
nominal_width = 1700
itf.redefine_domain(
    x_min = (x_size-nominal_length)/2,
    x_max = (x_size+nominal_length)/2,
    y_min = (y_size-nominal_width)/2,
    y_max = (y_size+nominal_width)/2,
)
itf.export_bounds()

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "1_ebsd_original",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Process the EBSD map
itf.remove_grains(2500)
itf.decrease_resolution(4)
itf.fill(10)
itf.clean(3)
itf.smooth(3)
itf.add_grips(30)

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "2_ebsd_processed",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Mesh the EBSD map and plot
itf.mesh("~/cubit/psculpt.exe", z_elements=3)
itf.export_grains(degrees=False)
itf.export_elements(degrees=False)
itf.plot_mesh(
    mesh_path = "3_mesh_raw",
    ipf       = "x",
    figure_x  = 20
)

# Fix the interfaces of the grips and scale the mesh
itf.fix_grip_interfaces(nominal_length/3, nominal_length)
itf.scale_mesh(nominal_width, "y")
itf.scale_mesh(300, "z")
itf.plot_mesh(
    mesh_path = "4_mesh_fixed",
    ipf       = "x",
    figure_x  = 20
)
