"""
 Title:         Controller
 Description:   For conducting the meshing
 Author:        Janzen Choi

"""

# Libraries
import math, numpy as np
from copy import deepcopy
from ebsd_mesher.visualiser.ebsd_plotter import EBSDPlotter
from ebsd_mesher.visualiser.mesh_plotter import MeshPlotter
from ebsd_mesher.mapper.gridder import read_elements, get_grain_ids, get_void_element, get_void_element_grid
from ebsd_mesher.mapper.gridder import get_areas, VOID_ELEMENT_ID
from ebsd_mesher.mapper.improver import clean_element_grid, smoothen_edges, pad_edges, remove_small_grains
from ebsd_mesher.mesher.exodus import map_spn_to_exo, renumber_grains, get_exodus_dimension, scale_exodus_mesh
from ebsd_mesher.mesher.exodus import get_element_info, straighten_interface, scale_gripped_microstructure
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
        self.headers       = []
        self.element_grid  = None
        self.step_size     = None
        self.output_dir    = output_dir
        self.input_path    = f"{output_dir}/input.i"
        self.spn_path      = f"{output_dir}/voxels.spn"
        self.exodus_path   = f"{output_dir}/mesh.e"
        self.mesh_elements = None
        self.grip_ids      = [None, None] # left, right
        self.spn_to_exo    = None

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

    def import_ebsd(self, ebsd_path:str, step_size:float, degrees:bool=True) -> None:
        """
        Reads in an EBSD map

        Parameters:
        * `ebsd_path`: Path to the EBSD file as a CSV file
        * `step_size`: Step size between coordinates
        * `degrees`:   Whether to the orientations are in degrees
        """
        self.element_grid = read_elements(ebsd_path, step_size, self.headers, degrees)
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
        x_size = len(self.element_grid[0])
        y_size = len(self.element_grid)
        new_x_size = x_max - x_min
        new_y_size = y_max - y_min

        # Create new element grid and replace
        new_element_grid = get_void_element_grid(new_x_size, new_y_size)
        for row in range(y_size):
            for col in range(x_size):
                new_col, new_row = abs(col-x_min), abs(row-y_min)
                if new_col >= 0 and new_row >= 0 and new_col < new_x_size and new_row < new_y_size:
                    new_element_grid[new_row][new_col] = deepcopy(self.element_grid[row][col])
        self.element_grid = new_element_grid

    def decrease_resolution(self, factor:int) -> None:
        """
        Decreases the resolution of the voxellation
        
        Parameters:
        * `factor`: The factor of the resolution decrease
        """
        self.step_size *= factor
        new_x_size = math.ceil(len(self.element_grid[0]) / factor)
        new_y_size = math.ceil(len(self.element_grid) / factor)
        new_element_grid = get_void_element_grid(new_x_size, new_y_size)
        for row in range(new_y_size):
            for col in range(new_x_size):
                new_element_grid[row][col] = deepcopy(self.element_grid[row*factor][col*factor])
        self.element_grid = new_element_grid

    def clean(self, iterations:int=1) -> None:
        """
        Cleans the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the cleaning
        """
        for _ in range(iterations):
            self.element_grid = clean_element_grid(self.element_grid)

    def smooth(self, iterations:int=1) -> None:
        """
        Smoothens the edges of the grains in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the smoothing
        """
        for _ in range(iterations):
            self.element_grid = smoothen_edges(self.element_grid)
    
    def fill(self, iterations:int=1) -> None:
        """
        Fills in the voids in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the filling
        """
        for _ in range(iterations):
            self.element_grid = pad_edges(self.element_grid)
    
    def remove_grains(self, threshold:float) -> None:
        """
        Removes small grains in the EBSD map

        Parameters:
        * `threshold`: The threshold grain area to remove the grains
        """
        threshold /= self.step_size**2
        self.element_grid = remove_small_grains(self.element_grid, threshold)

    def add_grips(self, num_cells:int) -> None:
        """
        Adds grips to the left and right of the microstructure
        
        Parameters:
        * `num_cells`: The number of cells (in the x-direction) to use in the grips
        """
        
        # Create grain IDs for the grips
        grain_ids = get_grain_ids(self.element_grid)
        l_grip_id = max(grain_ids+[VOID_ELEMENT_ID])+1
        r_grip_id = max(grain_ids+[VOID_ELEMENT_ID])+2
        self.grip_ids[0] = l_grip_id
        self.grip_ids[1] = r_grip_id

        # Add grip IDs to element grid
        new_element_grid = deepcopy(self.element_grid)
        for row in range(len(new_element_grid)):
            l_elements = [get_void_element(l_grip_id) for _ in range(num_cells)]
            r_elements = [get_void_element(r_grip_id) for _ in range(num_cells)]
            new_element_grid[row] = l_elements + new_element_grid[row] + r_elements
        self.element_grid = new_element_grid

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
        plotter = EBSDPlotter(self.element_grid, self.step_size, figure_x)
        plotter.plot_ebsd(ipf)

        # Add IDs and boundaries if specified
        id_list = get_grain_ids(self.element_grid) if id_list == None else id_list
        if (isinstance(grain_id, bool) and grain_id) or isinstance(grain_id, dict):
            plotter.plot_ids(id_list, settings=id_settings)
        if (isinstance(boundary, bool) and boundary) or isinstance(boundary, dict):
            plotter.plot_boundaries(id_list, settings=boundary_settings)

        # Save
        ebsd_path = get_file_path_exists(ebsd_path, "png")
        save_plot(ebsd_path)

    def mesh(self, psculpt_path:str, z_elements:int, num_processors:int) -> None:
        """
        Generates a mesh based on an SPN file
        
        Parameters:
        * `psculpt_path`:   The path to PSculpt 
        * `z_elements`:     The number of elements in the the mesh's thickness
        * `num_processors`: The number of processors to use to create the mesh
        """

        # Generate the mesh and renumber the meshed grains (ascending and consecutive)
        print()
        coarse_mesh(psculpt_path, z_elements, num_processors, self.element_grid,
                    self.step_size, self.input_path, self.spn_path, self.exodus_path)
        print()
        renumber_grains(self.exodus_path)

        # Map the grains of the EBSD map to the mesh
        spn_size = (len(self.element_grid[0]), len(self.element_grid), z_elements)
        self.spn_to_exo, confidence_list = map_spn_to_exo(self.exodus_path, self.spn_path, spn_size)

        # Save element information
        self.mesh_elements = get_element_info(self.exodus_path, self.element_grid, self.step_size)

        # Save the mapping
        map_dict = {
            "ebsd_id":    list(self.spn_to_exo.keys()),
            "mesh_id":    list(self.spn_to_exo.values()),
            "confidence": confidence_list
        }
        dict_to_csv(map_dict, f"{self.output_dir}/grain_map.csv")

    def export_grains(self, degrees:bool=True) -> None:
        """
        Exports the orientations and areas of each grain
        
        Parameters:
        * `degrees`: Whether to save the orientations as degrees
        """
        
        # Initialise
        exo_to_spn = dict(zip(self.spn_to_exo.values(), self.spn_to_exo.keys()))
        stats_dict = {"phi_1": [], "Phi": [], "phi_2": [], "area": []}
        area_dict = get_areas(self.element_grid)
        
        # Iterate through each grain
        for exo_id in exo_to_spn.keys():
            
            # Get grain ID and check if void
            spn_id = exo_to_spn[exo_id]
            if spn_id in [VOID_ELEMENT_ID]:
                continue

            # Get average orientation
            phi_1_list, Phi_list, phi_2_list = [], [], []
            for row in range(len(self.element_grid)):
                for col in range(len(self.element_grid[row])):
                    if spn_id == self.element_grid[row][col].get_grain_id():
                        phi_1, Phi, phi_2 = self.element_grid[row][col].get_orientation(degrees)
                        phi_1_list.append(phi_1)
                        Phi_list.append(Phi)
                        phi_2_list.append(phi_2)
            
            # Store statistics
            stats_dict["phi_1"].append(np.average(phi_1_list))
            stats_dict["Phi"].append(np.average(Phi_list))
            stats_dict["phi_2"].append(np.average(phi_2_list))
            stats_dict["area"].append(area_dict[spn_id])
        
        # Save statistics
        dict_to_csv(stats_dict, f"{self.output_dir}/grain_stats.csv", add_header=False)

    def export_elements(self, degrees:bool=True) -> None:
        """
        Exports the orientation and grain ID of each element
        
        Parameters:
        * `degrees`: Whether to save the orientations as degrees
        """
        
        # Initialise
        stats_dict = {"phi_1": [], "Phi": [], "phi_2": [], "grain_id": []}
        
        # Iterate through ordered elements
        for element in self.mesh_elements:
            phi_1, Phi, phi_2 = element.get_orientation(degrees)
            stats_dict["phi_1"].append(phi_1)
            stats_dict["Phi"].append(Phi)
            stats_dict["phi_2"].append(phi_2)
            stats_dict["grain_id"].append(element.get_grain_id())
        
        # Save statistics
        dict_to_csv(stats_dict, f"{self.output_dir}/element_stats.csv", add_header=False)

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
        mesh_plotter = MeshPlotter(self.exodus_path, self.mesh_elements, figure_x)
        mesh_plotter.plot_mesh(ipf)
        mesh_path = get_file_path_exists(mesh_path, "png")
        save_plot(mesh_path)
