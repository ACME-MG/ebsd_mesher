"""
 Title:         Gridder
 Description:   Creates and manipulates pixel grids
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from copy import deepcopy
from ebsd_mesher.mapper.grain import Grain

# Special element IDs
VOID_PIXEL_ID  = 100000 # large number
NO_ORIENTATION = [0, 0, 0] # for void

def get_void_pixel_grid(x_cells:list, y_cells:list) -> list:
    """
    Creates a grid of void pixels
    
    Parameters:
    * `x_cells`:    The number of pixels on the horizontal axis
    * `y_cells`:    The number of pixels on the vertical axis
    * `init_value`: The initial value of the cell in the piel grid
    
    Returns a grid of void pixels
    """
    pixel_grid = []
    for _ in range(y_cells):
        pixel_list = []
        for _ in range(x_cells):
            pixel_list.append(VOID_PIXEL_ID)
        pixel_grid.append(pixel_list)
    return pixel_grid

def get_info(value_list:list, step_size:float) -> tuple:
    """
    Gets the range and step size from a list of values
    
    Parameters:
    * `value_list`: List of values
    * `step_size`:  The step size of each pixel
    
    Returns the number of values and minimum values
    """
    max_value = max(value_list)
    min_value = min(value_list)
    num_values = round((max_value - min_value) / step_size) + 1
    return num_values, min_value

def read_pixels(path:str, step_size:float, headers:list=None) -> tuple:
    """
    Converts a CSV file into a grid of pixels
    
    Parameters:
    * `path`:      The path to the CSV file
    * `step_size`: The step size of each pixel
    * `headers`:   The list of headers for the x coordinates, y coordinates,
                   grain IDs, phi_1, Phi, and phi_2 values
    
    Returns the pixel grid and grain map
    """

    # Open file and read headers
    file = open(path, "r")
    csv_headers = file.readline().replace("\n", "").split(",")
    rows = file.readlines()
    
    # Get column indexes
    if headers == None:
        headers = ["x", "y", "grain_id", "phi_1", "Phi", "phi_2"]
    index_x         = csv_headers.index(headers[0])
    index_y         = csv_headers.index(headers[1])
    index_grain_id  = csv_headers.index(headers[2])
    index_avg_phi_1 = csv_headers.index(headers[3])
    index_avg_Phi   = csv_headers.index(headers[4])
    index_avg_phi_2 = csv_headers.index(headers[5])

    # Get dimensions
    x_cells, x_min = get_info([float(row.split(",")[index_x]) for row in rows], step_size)
    y_cells, y_min = get_info([float(row.split(",")[index_y]) for row in rows], step_size)
    
    # Initialise pixel grid and grain map
    pixel_grid = get_void_pixel_grid(x_cells, y_cells)
    grain_map = {}

    # Read CSV and fill grid
    for row in rows:

        # Process data
        row_list = row.replace("\n", "").split(",")
        if "NaN" in row_list or "nan" in row_list:
            continue
        row_list = [float(val) for val in row_list]
        grain_id = round(row_list[index_grain_id])

        # Add to pixel grid
        x = round(float(row_list[index_x] - x_min) / step_size)
        y = round(float(row_list[index_y] - y_min) / step_size)
        pixel_grid[y][x] = grain_id

        # Add to grain map if not yet added
        if not grain_id in grain_map:
            new_grain = Grain(
                phi_1    = row_list[index_avg_phi_1],
                Phi      = row_list[index_avg_Phi],
                phi_2    = row_list[index_avg_phi_2],
                size     = 1,
            )
            grain_map[grain_id] = new_grain
        
        # Update grain map if already added
        else:
            grain_map[grain_id].increment_size()
    
    # Close file and return grid and map
    file.close()
    return pixel_grid, grain_map

def get_grain_ids(pixel_grid:list) -> list:
    """
    Gets the grain IDs of the pixel grid
    
    Parameters:
    * `pixel_grid`: A grid of pixels
    
    Returns the list of grain IDs
    """
    grain_ids = list(set([pixel for pixel_list in pixel_grid for pixel in pixel_list]))
    if VOID_PIXEL_ID in grain_ids:
        grain_ids.remove(VOID_PIXEL_ID)
    return grain_ids

def shift_pixel_grid(pixel_grid:list, shift_amount:int) -> list:
    """
    Shifts the grain IDs of the pixel grid
    
    Parameters:
    * `pixel_grid`:   A grid of pixels
    * `shift_amount`: Amount to shift the grain IDs
    
    Returns the shifted pixel grid
    """
    new_pixel_grid = deepcopy(pixel_grid)
    for row in range(len(pixel_grid)):
        for col in range(len(pixel_grid[0])):
            new_pixel_grid[row][col] = pixel_grid[row][col] + shift_amount
    return new_pixel_grid

def shift_grain_map(grain_map:dict, shift_amount:int) -> dict:
    """
    Shifts the grain IDs of the grain map
    
    Parameters:
    * `grain_map`:    A mapping of the grains to the average orientations
    * `shift_amount`: Amount to shift the grain IDs
    
    Returns the shifted grain map
    """
    new_grain_map = {}
    for grain_id in grain_map.keys():
        new_grain_map[grain_id+shift_amount] = grain_map[grain_id]
    return new_grain_map

def get_centroids(pixel_grid:list) -> dict:
    """
    Gets the centroids for all the grains
    
    Parameters:
    * `pixel_grid`: A grid of pixels
    
    Returns a dictionary of the grain IDs and their centroids
    """
    
    # Initialise dictionary to store pixel positions
    flattened = [pixel for pixel_list in pixel_grid for pixel in pixel_list]
    grain_ids = list(set(flattened))
    pixel_dict = {}
    for grain_id in grain_ids:
        pixel_dict[grain_id] = {"x": [], "y": []}
    
    # Add points to the pixel dictionary
    for row in range(len(pixel_grid)):
        for col in range(len(pixel_grid[row])):
            pixel_dict[pixel_grid[row][col]]["x"].append(col)
            pixel_dict[pixel_grid[row][col]]["y"].append(row)
    
    # Calculate centroids
    centroid_dict = {}
    for grain_id in grain_ids:
        x_mean = np.average(pixel_dict[grain_id]["x"])
        y_mean = np.average(pixel_dict[grain_id]["y"])
        centroid_dict[grain_id] = (x_mean, y_mean)
    return centroid_dict
