"""
 Title:         EBSD Plotter
 Description:   Visualises a EBSD map using a plot
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
from ebsd_mesher.mapper.gridder import get_centroids
from ebsd_mesher.maths.ipf_cubic import euler_to_rgb
from ebsd_mesher.visualiser.plotter import get_coordinate, get_boundary

# EBSDPlotter class
class EBSDPlotter:
    
    def __init__(self, pixel_grid:list, grain_map:dict, step_size:float, figure_x:float=10):
        """
        Constructor for the plotter class
        
        Parameters:
        * `pixel_grid`: A grid of pixels
        * `grain_map`:  A mapping of the grains to the average orientations
        * `step_size`:  The step size of the EBSD map
        * `figure_x`:   Size of the horizontal axis of the figure (in inches)
        """
        
        # Initialise internal variables
        self.pixel_grid = pixel_grid
        self.grain_map = grain_map
        self.step_size = step_size
        self.figure_x = figure_x
        
        # Initialise plot
        x_max = len(pixel_grid[0])*self.step_size
        y_max = len(pixel_grid)*self.step_size
        self.figure, self.axis = plt.subplots(figsize=(figure_x, y_max/x_max*figure_x))
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        self.axis.invert_yaxis()
        
        # Define size of each square marker (55 magical number)
        self.square_size = 55*self.figure_x*self.step_size/x_max

    def plot_ebsd(self, ipf:str="x") -> None:
        """
        Plots the EBSD map using Matplotlib
        
        Parameters:
        * `pixel_grid`: A grid of pixels
        * `grain_map`:  A mapping of the grains to the average orientations
        * `step_size`:  The step size of the EBSD map
        * `ipf`:        The IPF direction
        """
        
        # Create colour map
        colour_map = {}
        for grain_id in self.grain_map.keys():
            orientation = self.grain_map[grain_id].get_orientation()
            colour = [rgb/255 for rgb in euler_to_rgb(*orientation, ipf=ipf)]
            colour_map[grain_id] = colour
        
        # Prepare pixel data
        x_list, y_list, colour_list = [], [], []
        for row in range(len(self.pixel_grid)):
            for col in range(len(self.pixel_grid[row])):
                if self.pixel_grid[row][col] in colour_map.keys():
                    x_list.append(get_coordinate(col, self.step_size))
                    y_list.append(get_coordinate(row, self.step_size))
                    colour_list.append(colour_map[self.pixel_grid[row][col]])

        # Plot
        plt.scatter(x_list, y_list, c=colour_list, s=self.square_size**2, marker="s")

    def plot_border(self) -> None:
        """
        Draws an outline at the outermost pixels;
        for checking that the whole microstructre is shown
        """
        x_list, y_list, colour_list = [], [], []
        for row in range(len(self.pixel_grid)):
            for col in range(len(self.pixel_grid[row])):
                if (col == 0 or col == len(self.pixel_grid[row])-1 or
                    row == 0 or row == len(self.pixel_grid)-1):
                        x_list.append(get_coordinate(col, self.step_size))
                        y_list.append(get_coordinate(row, self.step_size))
                        colour_list.append((0,0,0))
        plt.scatter(x_list, y_list, c=colour_list, s=self.square_size**2, marker="s")

    def plot_ids(self, id_list:list=None, settings:dict={}) -> None:
        """
        Writes the grain IDs at the centroids
        
        Parameters:
        * `id_list`: List of grain IDs to add centroids to
        """
        centroid_dict = get_centroids(self.pixel_grid)
        for grain_id in centroid_dict.keys():
            if id_list != None and not grain_id in id_list:
                continue
            x, y = centroid_dict[grain_id]
            x = get_coordinate(x, self.step_size)
            y = get_coordinate(y, self.step_size)
            plt.text(x, y, str(grain_id), ha="center", va="center", **settings)

    def plot_boundaries(self, id_list:list=None, settings:dict={}) -> None:
        """
        Plots the grain boundaries
        
        Parameters:
        * `id_list`: List of grain IDs to draw boundaries around
        """
        x_list, y_list = [], []
        for row in range(len(self.pixel_grid)):
            for col in range(len(self.pixel_grid[row])):
                if id_list == None or self.pixel_grid[row][col] in id_list:
                    x_boundary, y_boundary = get_boundary(row, col, self.pixel_grid, self.step_size)
                    x_list += x_boundary
                    y_list += y_boundary
        plt.plot(x_list, y_list, **settings)
