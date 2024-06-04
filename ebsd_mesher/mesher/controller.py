"""
 Title:         Controller
 Description:   For conducting the meshing
 Author:        Janzen Choi

"""

# Libraries
import math
from ebsd_mesher.visualiser.ebsd_plotter import EBSDPlotter
from ebsd_mesher.visualiser.mesh_plotter import MeshPlotter
from ebsd_mesher.mapper.gridder import read_pixels, get_grain_ids, get_void_pixel_grid
from ebsd_mesher.mapper.improver import clean_pixel_grid, smoothen_edges, pad_edges, remove_small_grains
from ebsd_mesher.maths.statistics import map_spn_to_exo
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
        self.headers     = []
        self.pixel_grid  = None
        self.grain_map   = None
        self.step_size   = None
        self.input_path  = f"{output_dir}/input.i"
        self.spn_path    = f"{output_dir}/voxels.spn"
        self.exodus_path = f"{output_dir}/mesh.e"
        self.map_path    = f"{output_dir}/grain_map.csv"
        self.thickness   = None
        self.spn_to_exo  = None

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

    def mesh(self, psculpt_path:str, thickness:int, num_processors:int) -> None:
        """
        Generates a mesh based on an SPN file
        
        Parameters:
        * `psculpt_path`:   The path to PSculpt 
        * `thickness`:      The thickness of the mesh (in voxels) 
        * `num_processors`: The number of processors to use to create the mesh
        """

        # Generate the mesh
        print()
        self.thickness = thickness
        coarse_mesh(psculpt_path, thickness, num_processors, self.pixel_grid,
                    self.step_size, self.input_path, self.spn_path, self.exodus_path)
        print()
        
        # Map the grains of the EBSD map to the mesh
        spn_size = (len(self.pixel_grid[0]), len(self.pixel_grid), self.thickness)
        self.spn_to_exo, confidence_list = map_spn_to_exo(self.exodus_path, self.spn_path, spn_size)

        # Save the mapping
        map_dict = {
            "ebsd_id":    list(self.spn_to_exo.keys()),
            "mesh_id":    list(self.spn_to_exo.values()),
            "confidence": confidence_list
        }
        dict_to_csv(map_dict, self.map_path)

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
