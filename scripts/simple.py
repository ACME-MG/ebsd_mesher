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

# Read EBSD map
DATA_FOLDER = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data"
FILE_NAME = "2024-06-07 (ansto_617_s1_simple)/simple_converted.csv"
itf.import_ebsd(
    ebsd_path = f"{DATA_FOLDER}/{FILE_NAME}",
    step_size = 20,
    degrees   = True
)
itf.add_grips(30)

# Process and plot EBSD map
itf.plot_ebsd(
    ebsd_path = "original",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Mesh the EBSD map and plot
itf.mesh("~/cubit/psculpt.exe", z_elements=1)
itf.export_grains(degrees=True)
itf.export_elements(degrees=True)
itf.plot_mesh(
    mesh_path = "raw_mesh",
    ipf       = "x",
    figure_x  = 20
)

# Fix the interfaces of the grips and scale the mesh
itf.fix_grip_interfaces(2300/4, 2300)
itf.scale_mesh(1571, "y")
# itf.scale_mesh(323, "z")
itf.plot_mesh(
    mesh_path = "fixed_mesh",
    ipf       = "x",
    figure_x  = 20
)
