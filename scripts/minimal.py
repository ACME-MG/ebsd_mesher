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
itf.import_ebsd(
    ebsd_path = f"data/minimal.csv",
    # ebsd_path = f"data/minimal_sachs.csv",
    # ebsd_path = f"data/minimal_taylor.csv",
    step_size = 1,
    degrees   = True
)
itf.increase_resolution(3)

# Add grips and plot EBSD map
itf.add_grips(1)
itf.plot_ebsd(
    ebsd_path = "1_ebsd",
    ipf       = "x",
    figure_x  = 20,
    grain_id  = {"fontsize": 10, "color": "black"},
    boundary  = True,
)

# Mesh the EBSD map and plot
itf.mesh("~/cubit/psculpt.exe", z_elements=1)
itf.fix_grip_interfaces(0.25, 2)
itf.export_grains(degrees=False)
itf.export_elements(degrees=False)
itf.plot_mesh(
    mesh_path = "2_mesh",
    ipf       = "x",
    figure_x  = 20
)
