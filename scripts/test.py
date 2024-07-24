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
itf.import_ebsd(
    ebsd_path = f"data/example.csv",
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
itf.export_grains(degrees=False)
itf.export_elements(degrees=False)
itf.plot_mesh(
    mesh_path = "raw_mesh",
    ipf       = "x",
    figure_x  = 20
)
