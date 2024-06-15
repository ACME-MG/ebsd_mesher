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

# Read EBSD data
EBSD_FOLDER = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data/2024-05-23 (ansto_617_s1)/hagb15/"
FILE_PATH = "01_strain_0pct_on_stage_finalMapData20/ebsdExportColumnsTable_NoFill.csv"
itf.import_ebsd(
    ebsd_path = f"{EBSD_FOLDER}/{FILE_PATH}",
    step_size = 5,
    degrees   = True
)

# Redefine the EBSD domain
bounds = itf.get_bounds()
x_size = bounds["max_x"]-bounds["min_x"]
y_size = bounds["max_y"]-bounds["min_y"]
itf.export_bounds()
itf.redefine_domain(bounds["min_x"], bounds["max_x"], bounds["min_y"], bounds["max_y"])

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "1_ebsd_original",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Process the EBSD map
itf.remove_grains(10000)
itf.decrease_resolution(4)
itf.fill(20)
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
itf.mesh("~/cubit/psculpt.exe", z_elements=1)
itf.export_grains(degrees=False)
itf.export_elements(degrees=False)
itf.plot_mesh(
    mesh_path = "3_mesh_raw",
    ipf       = "x",
    figure_x  = 20
)

# Fix the interfaces of the grips and scale the mesh
itf.fix_grip_interfaces(x_size/3, x_size)
# itf.scale_mesh(y_size, "y")
# itf.scale_mesh(323, "z")
itf.plot_mesh(
    mesh_path = "4_mesh_fixed",
    ipf       = "x",
    figure_x  = 20
)
