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

    def import_ebsd(self, ebsd_path:str, step_size:float) -> None:
        """
        Reads in an EBSD map

        Parameters:
        * `ebsd_path`: Path to the EBSD file as a CSV file
        * `step_size`: Step size between coordinates
        """
        self.__print__(f"Adding the EBSD map")
        self.__controller__.import_ebsd(ebsd_path, step_size)

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
        Decreases the resolution of the pixellation;
        define the factor as a fraction to increase the resolution
        
        Parameters:
        * `factor`: The factor of the resolution reduction
        """
        self.__print__(f"Decreasing the sample resolution by a factor of {factor}")
        self.__check_ebsd__()
        if factor < 0:
            raise ValueError("The 'factor' parameter needs to be greater than 0!")
        self.__controller__.decrease_resolution(factor)

    def clean(self, iterations:int=1) -> None:
        """
        Cleans the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the cleaning
        """
        self.__print__(f"Cleaning the EBSD map {iterations} time(s)")
        self.__check_ebsd__()
        self.__controller__.clean(iterations)

    def smooth(self, iterations:int=1) -> None:
        """
        Smoothens the edges of the grains in the EBSD map

        Parameters:
        * `iterations`: The number of times to conduct the smoothing
        """
        self.__print__(f"Smoothing the EBSD map {iterations} time(s)")
        self.__check_ebsd__()
        self.__controller__.smooth(iterations)
    
    def fill(self, iterations:int=1) -> None:
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
        * `id_list`:   The IDs of the grains to plot the IDs and boundaries;
                       IDs are the ones of the first grain map;
                       if undefined, adds for all grains
        """
        self.__print__(f"Plotting the EBSD map")
        self.__check_ebsd__()
        ebsd_path = self.__get_output__(ebsd_path)
        self.__controller__.plot_ebsd(ebsd_path, ipf, figure_x, grain_id, boundary, id_list)

    def mesh(self, psculpt_path:str, z_length:float, z_voxels:int=1, num_processors:int=1) -> None:
        """
        Generates a mesh based on an SPN file
        
        Parameters:
        * `psculpt_path`:   The path to PSculpt
        * `z_length`:       The thickness of the mesh (in units)
        * `z_voxels`:       The thickness of the mesh (in voxels)
        * `num_processors`: The number of processors to use to create the mesh
        """
        self.__print__("Generating a mesh of the EBSD map")
        self.__check_ebsd__()
        self.__controller__.mesh(psculpt_path, z_length, z_voxels, num_processors)

    def plot_mesh(self, mesh_path:str="mesh", ipf:str="x", figure_x:float=10) -> None:
        """
        Plots the generated mesh;
        relies on the mesh to have already been generated

        Parameters:
        * `mesh_path`: Path to generate the plot
        * `ipf`:       The IPF colour ("x", "y", "z")
        * `figure_x`:  The initial horizontal size of the figures
        """
        self.__print__("Generating an image of the mesh")
        self.__check_mesh__()
        mesh_path = self.__get_output__(mesh_path)
        self.__controller__.plot_mesh(mesh_path, ipf, figure_x)

    def __check_ebsd__(self) -> None:
        """
        Checks that the EBSD map has been imported in
        """
        if self.__controller__.pixel_grid == None:
            raise ValueError("The EBSD map has not been imported yet!")

    def __check_mesh__(self) -> None:
        """
        Checks that the EBSD mesh has been generated
        """
        if not self.__controller__.has_meshed:
            raise ValueError("The EBSD mesh has not been generated yet!")

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
