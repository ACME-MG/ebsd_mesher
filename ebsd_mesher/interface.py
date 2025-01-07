"""
 Title:         Interface
 Description:   Interface for mapping grains between EBSD maps
 Author:        Janzen Choi

"""

# Libraries
import inspect, re, time
from ebsd_mesher.mesher.controller import Controller
from ebsd_mesher.helper.io import safe_mkdir

# Interface Class
class Interface:

    def __init__(self, title:str="", output_path:str="./results", verbose:bool=True, output_here:bool=False):
        """
        Class to interact with the optimisation code
        
        Parameters:
        * `title`:       Title of the output folder
        * `input_path`:  Path to the input folder
        * `output_path`: Path to the output folder
        * `verbose`:     If true, outputs messages for each function call
        * `output_here`: If true, just dumps the output in ths executing directory
        """

        # Initialise internal variables
        self.__print_index__ = 0
        self.__verbose__     = verbose

        # Initialise timing variables
        time_str = time.strftime("%A, %D, %H:%M:%S", time.localtime())
        self.__print__(f"\n  Starting on {time_str}\n", add_index=False)
        self.__start_time__ = time.time()
        time_stamp = time.strftime("%y%m%d%H%M%S", time.localtime(self.__start_time__))
        
        # Define input and output
        file_path = inspect.currentframe().f_back.f_code.co_filename
        file_name = file_path.split("/")[-1].replace(".py","")
        title = f"_{file_name}" if title == "" else f"_{title}"
        title = re.sub(r"[^a-zA-Z0-9_]", "", title.replace(" ", "_"))
        output_dir = "." if output_here else f"{output_path}/{time_stamp}{title}"
        self.__get_output__ = lambda x : f"{output_dir}/{x}"
        
        # Initialise controller
        self.__controller__  = Controller(output_dir)
        self.__controller__.define_headers("x", "y", "grain_id", "phi_1", "Phi", "phi_2")
        
        # Create directories
        if not output_here:
            safe_mkdir(output_path)
            safe_mkdir(output_dir)    

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
        self.__print__("Defining headers for EBSD files")
        self.__controller__.define_headers(x, y, grain_id, phi_1, Phi, phi_2)

    def import_ebsd(self, ebsd_path:str, step_size:float, degrees:bool=True) -> None:
        """
        Reads in an EBSD map

        Parameters:
        * `ebsd_path`: Path to the EBSD file as a CSV file
        * `step_size`: Step size between coordinates
        * `degrees`:   Whether to the orientations are in degrees
        """
        self.__print__(f"Adding the EBSD map")
        self.__controller__.import_ebsd(ebsd_path, step_size, degrees)

    def get_bounds(self) -> dict:
        """
        Returns the bounds of the imported EBSD map
        """
        self.__print__(f"Getting the bounds of the EBSD map")
        bounds = self.__controller__.get_bounds()
        return bounds
    
    def export_bounds(self) -> None:
        """
        Exports the bounds
        """
        self.__print__(f"Exporting the bounds of the EBSD map")
        self.__controller__.export_bounds()

    def redefine_domain(self, x_min:float, x_max:float, y_min:float, y_max:float) -> None:
        """
        Redefines the domain
        
        Parameters:
        * `x_min`: The lowest x value for the new domain
        * `x_max`: The highest x value for the new domain
        * `y_min`: The lowest y value for the new domain
        * `y_max`: The highest y value for the new domain
        """
        self.__print__("Redefining the domain")
        self.__check_ebsd__()
        self.__controller__.redefine_domain(x_min, x_max, y_min, y_max)

    def decrease_resolution(self, factor:int) -> None:
        """
        Decreases the resolution of the pixellated EBSD map
        
        Parameters:
        * `factor`: The factor of the resolution reduction
        """
        self.__print__(f"Decreasing the sample resolution by a factor of {factor}")
        self.__check_ebsd__()
        if factor < 0:
            raise ValueError("The 'factor' parameter needs to be greater than 0!")
        self.__controller__.decrease_resolution(factor)

    def increase_resolution(self, factor:int) -> None:
        """
        Increases the resolution of the pixellated EBSD map
                
        Parameters:
        * `factor`: The factor of the resolution increase
        """
        self.__print__(f"Increasing the sample resolution by a factor of {factor}")
        self.__check_ebsd__()
        if factor < 0:
            raise ValueError("The 'factor' parameter needs to be greater than 0!")
        self.__controller__.increase_resolution(factor)

    def clean(self, iterations:int) -> None:
        """
        Cleans the EBSD map; can mess up equivalent radius distribution

        Parameters:
        * `iterations`: The number of times to conduct the cleaning
        """
        self.__print__(f"Cleaning the EBSD map {iterations} time(s)")
        self.__check_ebsd__()
        self.__controller__.clean(iterations)

    def smooth(self, iterations:int) -> None:
        """
        Smoothens the edges of the grains in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the smoothing
        """
        self.__print__(f"Smoothing the EBSD map {iterations} time(s)")
        self.__check_ebsd__()
        self.__controller__.smooth(iterations)
    
    def fill(self, iterations:int) -> None:
        """
        Fills in the voids in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the filling
        """
        self.__print__(f"Filling the voids of the EBSD map {iterations} time(s)")
        self.__check_ebsd__()
        self.__controller__.fill(iterations)
    
    def remove_grains(self, threshold:float) -> None:
        """
        Removes small grains in the EBSD map

        Parameters:
        * `threshold`: The threshold grain area to remove the grains
        """
        self.__print__(f"Removing grains from the EBSD map with areas under {threshold}u2")
        self.__check_ebsd__()
        self.__controller__.remove_grains(threshold)

    def add_grips(self, num_elements:int) -> None:
        """
        Adds grips to the left and right of the microstructure
        
        Parameters:
        * `num_elements`: The number of elements (in the x-direction) to use in the grips
        """
        self.__print__(f"Adding grips to the left and right of the EBSD map")
        self.__check_ebsd__()
        self.__controller__.add_grips(num_elements)

    def plot_ebsd(self, ebsd_path:str="ebsd", ipf:str="x", figure_x:float=10,
                  grain_id:bool=False, boundary:bool=False, id_list:list=None) -> None:
        """
        Plots the EBSD maps

        Parameters:
        * `ebsd_path`: Path to save plot
        * `ipf`:       The IPF colour ("x", "y", "z")
        * `figure_x`:  The initial horizontal size of the figures
        * `grain_id`:  Whether to include IDs in the EBSD maps;
                       define dictionary for custom settings
        * `boundary`:  Whether to include IDs in the EBSD maps;
                       define dictionary for custom settings
        * `id_list`:   The IDs of the grains to plot the IDs and boundaries
                       if undefined, adds for all grains
        """
        self.__print__(f"Plotting the EBSD map")
        self.__check_ebsd__()
        ebsd_path = self.__get_output__(ebsd_path)
        self.__controller__.plot_ebsd(ebsd_path, ipf, figure_x, grain_id, boundary, id_list)

    def mesh(self, psculpt_path:str, z_elements:int=1, num_processors:int=1) -> None:
        """
        Generates a mesh based on an SPN file
        
        Parameters:
        * `psculpt_path`:   The path to PSculpt
        * `z_elements`:     The number of elements in the the mesh's thickness
        * `num_processors`: The number of processors to use to create the mesh
        """
        self.__print__(f"Generating a mesh of the EBSD map with a thickness of {z_elements} element(s)")
        self.__check_ebsd__()
        self.__controller__.mesh(psculpt_path, z_elements, num_processors)

    def export_grains(self, degrees:bool=True) -> None:
        """
        Exports the orientation and area of each grain
        
        Parameters:
        * `degrees`: Whether to save the orientations as degrees
        """
        self.__print__("Exporting the grain statistics from the mesh")
        self.__check_mesh__()
        self.__controller__.export_grains(degrees)

    def export_elements(self, degrees:bool=True) -> None:
        """
        Exports the orientation and grain ID of each element
        
        Parameters:
        * `degrees`: Whether to save the orientations as degrees
        """
        self.__print__("Exporting the element statistics from the mesh")
        self.__check_mesh__()
        self.__controller__.export_elements(degrees)

    def fix_grip_interfaces(self, grip_length:float, micro_length:float, straighten:bool=True) -> None:
        """
        Fixes the interfaces between the grips and the microstructure

        Parameters:
        * `grip_length`:  The desired length of the grip
        * `micro_length`: The desired length of the microstructure
        * `straighten`:   Whether to also straighten the grip interfaces
        """
        self.__print__("Fixing the interface between the grip and the microstructure")
        self.__check_grips__()
        self.__controller__.fix_grip_interfaces(grip_length, micro_length, straighten)

    def scale_mesh(self, length:float, direction:str) -> None:
        """
        Scales the mesh to a certain length along a certain direction

        Parameters:
        * `length`:    The length to scale the mesh to
        * `direction`: The direction to conduct the scaling ("x", "y", "z")
        """
        self.__print__(f"Scaling the mesh along the {direction}-axis to {length}")
        self.__check_mesh__()
        self.__controller__.scale_mesh(length, direction)

    def plot_mesh(self, mesh_path:str="mesh", ipf:str="x", figure_x:float=10,
                  directions:str="xy", positive:bool=True) -> None:
        """
        Plots the generated mesh;
        relies on the mesh to have already been generated

        Parameters:
        * `mesh_path`:  Path to generate the plot
        * `ipf`:        The IPF colour ("x", "y", "z")
        * `figure_x`:   The initial horizontal size of the figures
        * `directions`: The directions of to plot the mesh
        * `positive`:   Whether to plot the positive or negative face
        """
        self.__print__("Generating an image of the mesh")
        self.__check_mesh__()
        mesh_path = self.__get_output__(mesh_path)
        if len(directions) != 2 or not set(directions).issubset({"x", "y", "z"}):
            raise ValueError("The 'directions' field must contain two characters from 'x', 'y', and 'z'")
        self.__controller__.plot_mesh(mesh_path, ipf, figure_x, directions, positive)

    def __check_ebsd__(self) -> None:
        """
        Checks that the EBSD map has been imported in
        """
        if self.__controller__.element_grid == None:
            raise ValueError("The EBSD map has not been imported yet!")

    def __check_mesh__(self) -> None:
        """
        Checks that the EBSD mesh has been generated
        """
        if self.__controller__.mesh_grains == None:
            raise ValueError("The EBSD mesh has not been generated yet!")

    def __check_grips__(self) -> None:
        """
        Checks that the grips have been added
        """
        if None in self.__controller__.grip_ids:
            raise ValueError("The grips have not been added yet!")

    def __print__(self, message:str, add_index:bool=True, sub_index:bool=False) -> None:
        """
        Displays a message before running the command (for internal use only)
        
        Parameters:
        * `message`:   the message to be displayed
        * `add_index`: if true, adds a number at the start of the message
        * `sub_index`: if true, adds a number as a decimal
        """
        
        # Special printing cases
        if not self.__verbose__:
            return
        if not add_index:
            print(message)
            return
        
        # Prints with an index / subindex
        if sub_index:
            self.__print_subindex__ += 1
            print_index = f"     {self.__print_index__}.{self.__print_subindex__}"
        else:
            self.__print_index__ += 1
            self.__print_subindex__ = 0
            print_index = f"{self.__print_index__}"
        print(f"   {print_index})\t{message} ...")

    def __del__(self):
        """
        Prints out the final message (for internal use only)
        """
        time_str = time.strftime("%A, %D, %H:%M:%S", time.localtime())
        duration = round(time.time() - self.__start_time__)
        duration_h = duration // 3600
        duration_m = (duration - duration_h * 3600) // 60
        duration_s = duration - duration_h * 3600 - duration_m * 60
        duration_str_list = [
            f"{duration_h} hours" if duration_h > 0 else "",
            f"{duration_m} mins" if duration_m > 0 else "",
            f"{duration_s} seconds" if duration_s > 0 else ""
        ]
        duration_str = ", ".join([d for d in duration_str_list if d != ""])
        duration_str = f"in {duration_str}" if duration_str != "" else ""
        self.__print__(f"\n  Finished on {time_str} {duration_str}\n", add_index=False)
