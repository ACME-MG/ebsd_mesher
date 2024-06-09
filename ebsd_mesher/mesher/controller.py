"""
 Title:         Controller
 Description:   For conducting the meshing
 Author:        Janzen Choi

"""

# Libraries
import math
from copy import deepcopy
from ebsd_mesher.visualiser.ebsd_plotter import EBSDPlotter
from ebsd_mesher.visualiser.mesh_plotter import MeshPlotter
from ebsd_mesher.mapper.gridder import read_pixels, get_grain_ids, get_void_pixel_grid, VOID_PIXEL_ID
from ebsd_mesher.mapper.grain import Grain
from ebsd_mesher.mapper.improver import clean_pixel_grid, smoothen_edges, pad_edges, remove_small_grains
from ebsd_mesher.maths.orientation import deg_to_rad
from ebsd_mesher.mesher.exodus import map_spn_to_exo, renumber_grains, get_exodus_dimension, scale_exodus_mesh
from ebsd_mesher.mesher.exodus import straighten_interface, scale_gripped_microstructure
from ebsd_mesher.mesher.mesher import coarse_mesh
from ebsd_mesher.helper.io import get_file_path_exists, dict_to_csv
from ebsd_mesher.visualiser.plotter import save_plot

# Controller class
class Controller:

    def __init__(self, output_dir:str):
        """
        Initialises the constructor class

        Parameters:
        * `output_dir`: Path to the output directory
        """
        self.headers      = []
        self.pixel_grid   = None
        self.grain_map    = None
        self.step_size    = None
        self.input_path   = f"{output_dir}/input.i"
        self.spn_path     = f"{output_dir}/voxels.spn"
        self.exodus_path  = f"{output_dir}/mesh.e"
        self.map_path     = f"{output_dir}/grain_map.csv"
        self.ori_path     = f"{output_dir}/orientations.csv"
        self.summary_path = f"{output_dir}/summary.csv"
        self.has_meshed   = False
        self.grip_ids     = [None, None] # left, right
        self.spn_to_exo   = None

    def define_headers(self, x:str, y:str, grain_id:str, phi_1:str, Phi:str, phi_2:str) -> None:
        """
        Defines the necessary headers for the CSV files

        Parameters:
        * `x`:        Header for the x-coordinate
        * `y`:        Header for the y-coordinate
        * `grain_id`: Header for the grain ID
        * `phi_1`:    Header for the phi_1 values
        * `Phi`:      Header for the Phi values
        * `phi_2`:    Header for the phi_2 values
        """
        self.headers = [x, y, grain_id, phi_1, Phi, phi_2]

    def import_ebsd(self, ebsd_path:str, step_size:float) -> None:
        """
        Reads in an EBSD map

        Parameters:
        * `ebsd_path`: Path to the EBSD file as a CSV file
        * `step_size`: Step size between coordinates
        """
        pixel_grid, grain_map = read_pixels(ebsd_path, step_size, self.headers)
        self.pixel_grid = pixel_grid
        self.grain_map = grain_map
        self.step_size = step_size

    def redefine_domain(self, x_min:float, x_max:float, y_min:float, y_max:float) -> None:
        """
        Redefines the domain
        
        Parameters:
        * `x_min`: The lowest x value for the new domain
        * `x_max`: The highest x value for the new domain
        * `y_min`: The lowest y value for the new domain
        * `y_max`: The highest y value for the new domain
        """
        
        # Get boundaries
        x_min = round(x_min / self.step_size)
        x_max = round(x_max / self.step_size)
        y_min = round(y_min / self.step_size)
        y_max = round(y_max / self.step_size)

        # Get new and original lengths
        x_size = len(self.pixel_grid[0])
        y_size = len(self.pixel_grid)
        new_x_size = x_max - x_min
        new_y_size = y_max - y_min

        # Create new pixel grid and replace
        new_pixel_grid = get_void_pixel_grid(new_x_size, new_y_size)
        for row in range(y_size):
            for col in range(x_size):
                new_col, new_row = abs(col-x_min), abs(row-y_min)
                if new_col >= 0 and new_row >= 0 and new_col < new_x_size and new_row < new_y_size:
                    new_pixel_grid[new_row][new_col] = self.pixel_grid[row][col]
        self.pixel_grid = new_pixel_grid

    def decrease_resolution(self, factor:int) -> None:
        """
        Decreases the resolution of the voxellation
        
        Parameters:
        * `factor`: The factor of the resolution decrease
        """
        self.step_size *= factor
        new_x_size = math.ceil(len(self.pixel_grid[0]) / factor)
        new_y_size = math.ceil(len(self.pixel_grid) / factor)
        new_pixel_grid = get_void_pixel_grid(new_x_size, new_y_size)
        for row in range(new_y_size):
            for col in range(new_x_size):
                new_pixel_grid[row][col] = self.pixel_grid[row * factor][col * factor]
        self.pixel_grid = new_pixel_grid

    def clean(self, iterations:int=1) -> None:
        """
        Cleans the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the cleaning
        """
        for _ in range(iterations):
            self.pixel_grid = clean_pixel_grid(self.pixel_grid)

    def smooth(self, iterations:int=1) -> None:
        """
        Smoothens the edges of the grains in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the smoothing
        """
        for _ in range(iterations):
            self.pixel_grid = smoothen_edges(self.pixel_grid)
    
    def fill(self, iterations:int=1) -> None:
        """
        Fills in the voids in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the filling
        """
        for _ in range(iterations):
            self.pixel_grid = pad_edges(self.pixel_grid)
    
    def remove_grains(self, threshold:float) -> None:
        """
        Removes small grains in the EBSD map

        Parameters:
        * `threshold`: The threshold grain area to remove the grains
        """
        threshold /= self.step_size**2
        self.pixel_grid = remove_small_grains(self.pixel_grid, threshold)

    def add_grips(self, num_cells:int) -> None:
        """
        Adds grips to the left and right of the microstructure
        
        Parameters:
        * `num_cells`: The number of cells (in the x-direction) to use in the grips
        """
        
        # Create grain IDs for the grips
        grain_ids = get_grain_ids(self.pixel_grid)
        l_grip_id = max(grain_ids+[VOID_PIXEL_ID])+1
        r_grip_id = max(grain_ids+[VOID_PIXEL_ID])+2
        self.grip_ids[0] = l_grip_id
        self.grip_ids[1] = r_grip_id

        # Add grip IDs to pixel grid
        new_pixel_grid = deepcopy(self.pixel_grid)
        for row in range(len(new_pixel_grid)):
            new_pixel_grid[row] = [l_grip_id]*num_cells + new_pixel_grid[row] + [r_grip_id]*num_cells
        self.pixel_grid = new_pixel_grid

        # Add orientation to grips
        l_grip_grain = Grain(0, 0, 0, num_cells*len(new_pixel_grid))
        r_grip_grain = Grain(0, 0, 0, num_cells*len(new_pixel_grid))
        self.grain_map[l_grip_id] = l_grip_grain
        self.grain_map[r_grip_id] = r_grip_grain

    def plot_ebsd(self, ebsd_path:str, ipf:str="x", figure_x:float=10,
                  grain_id:bool=False, boundary:bool=False, id_list:list=None) -> None:
        """
        Plots the EBSD maps

        Parameters:
        * `ebsd_path`: The path to save the EBSD files
        * `ipf`:       The IPF colour ("x", "y", "z")
        * `figure_x`:  The initial horizontal size of the figures
        * `grain_id`:  Whether to include IDs in the EBSD maps;
                       define dictionary for custom settings
        * `boundary`:  Whether to include IDs in the EBSD maps;
                       define dictionary for custom settings
        * `id_list`:   The IDs of the grains to plot the IDs and boundaries;
                       if undefined, adds for all grains
        """

        # Define settings
        id_settings = grain_id if isinstance(grain_id, dict) else {"fontsize": 20, "color": "black"}
        boundary_settings = boundary if isinstance(boundary, dict) else {"linewidth": 1, "color": "black"}
        
        # Plot EBSD map
        plotter = EBSDPlotter(self.pixel_grid, self.grain_map, self.step_size, figure_x)
        plotter.plot_ebsd(ipf)

        # Add IDs and boundaries if specified
        id_list = get_grain_ids(self.pixel_grid) if id_list == None else id_list
        if (isinstance(grain_id, bool) and grain_id) or isinstance(grain_id, dict):
            plotter.plot_ids(id_list, settings=id_settings)
        if (isinstance(boundary, bool) and boundary) or isinstance(boundary, dict):
            plotter.plot_boundaries(id_list, settings=boundary_settings)

        # Save
        ebsd_path = get_file_path_exists(ebsd_path, "png")
        save_plot(ebsd_path)

    def mesh(self, psculpt_path:str, z_voxels:int, num_processors:int) -> None:
        """
        Generates a mesh based on an SPN file
        
        Parameters:
        * `psculpt_path`:   The path to PSculpt 
        * `z_voxels`:       The thickness of the mesh (in voxels) 
        * `num_processors`: The number of processors to use to create the mesh
        """

        # Generate the mesh and renumber the meshed grains (ascending and consecutive)
        print()
        self.has_meshed = True
        coarse_mesh(psculpt_path, z_voxels, num_processors, self.pixel_grid,
                    self.step_size, self.input_path, self.spn_path, self.exodus_path)
        print()
        renumber_grains(self.exodus_path)

        # Map the grains of the EBSD map to the mesh
        spn_size = (len(self.pixel_grid[0]), len(self.pixel_grid), z_voxels)
        self.spn_to_exo, confidence_list = map_spn_to_exo(self.exodus_path, self.spn_path, spn_size)

        # Save the mapping
        map_dict = {
            "ebsd_id":    list(self.spn_to_exo.keys()),
            "mesh_id":    list(self.spn_to_exo.values()),
            "confidence": confidence_list
        }
        dict_to_csv(map_dict, self.map_path)

        # Save orientations
        ori_dict = {"phi_1": [], "Phi": [], "phi_2": []}
        exo_to_spn = dict(zip(self.spn_to_exo.values(), self.spn_to_exo.keys()))
        for exo_id in exo_to_spn.keys():
            if exo_to_spn[exo_id] in [VOID_PIXEL_ID]:
                continue
            phi_1, Phi, phi_2 = self.grain_map[exo_to_spn[exo_id]].get_orientation()
            ori_dict["phi_1"].append(deg_to_rad(phi_1))
            ori_dict["Phi"].append(deg_to_rad(Phi))
            ori_dict["phi_2"].append(deg_to_rad(phi_2))
        dict_to_csv(ori_dict, self.ori_path, add_header=False)

    def fix_grip_interfaces(self, grip_length:float, micro_length:float) -> None:
        """
        Fixes the interface between the grip and the microstructure

        Parameters:
        * `grip_length`:  The desired length of the grip
        * `micro_length`: The desired length of the microstructure
        """
        l_grip_id = self.spn_to_exo[self.grip_ids[0]]
        r_grip_id = self.spn_to_exo[self.grip_ids[1]]
        straighten_interface(self.exodus_path, l_grip_id, left=True)
        straighten_interface(self.exodus_path, r_grip_id, left=False)
        scale_gripped_microstructure(self.exodus_path, l_grip_id, r_grip_id, grip_length, micro_length)

    def scale_mesh(self, length:float, direction:str="x") -> None:
        """
        Scales the mesh to a certain length along a certain direction

        Parameters:
        * `length`:    The length to scale the mesh to
        * `direction`: The direction to conduct the scaling
        """
        mesh_dimension = get_exodus_dimension(self.exodus_path, direction)
        factor = length / mesh_dimension
        scale_exodus_mesh(self.exodus_path, factor, direction)

    def plot_mesh(self, mesh_path:str, ipf:str="x", figure_x:float=10) -> None:
        """
        Plots the generated mesh;
        relies on the mesh to have already been generated

        Parameters:
        * `mesh_path`: Path to generate the plot
        * `ipf`:       The IPF colour ("x", "y", "z")
        * `figure_x`:  The initial horizontal size of the figures
        """
        mesh_plotter = MeshPlotter(self.exodus_path, self.pixel_grid, self.grain_map, self.step_size, figure_x)
        mesh_plotter.plot_mesh(self.spn_to_exo, ipf)
        mesh_path = get_file_path_exists(mesh_path, "png")
        save_plot(mesh_path)
