"""
 Title:         Gridder
 Description:   Creates and manipulates element grids
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from ebsd_mesher.mapper.element import Element

# Special element IDs
VOID_ORIENTATION = (0,0,0)
VOID_ELEMENT_ID = 100000 # large number

def get_void_element(grain_id:int=VOID_ELEMENT_ID) -> Element:
    """
    Returns a void element
    
    Parameters:
    * `grain_id`: The ID of the element
    """
    return Element(*VOID_ORIENTATION, grain_id=grain_id)

def get_void_element_grid(x_cells:list, y_cells:list) -> list:
    """
    Creates a grid of void elements
    
    Parameters:
    * `x_cells`:    The number of elements on the horizontal axis
    * `y_cells`:    The number of elements on the vertical axis
    * `init_value`: The initial value of the cell in the element grid
    
    Returns a grid of void elements
    """
    element_grid = []
    for _ in range(y_cells):
        element_list = []
        for _ in range(x_cells):
            element = get_void_element()
            element_list.append(element)
        element_grid.append(element_list)
    return element_grid

def get_info(value_list:list, step_size:float) -> tuple:
    """
    Gets the range and step size from a list of values
    
    Parameters:
    * `value_list`: List of values
    * `step_size`:  The step size of each element
    
    Returns the number of values and minimum values
    """
    max_value = max(value_list)
    min_value = min(value_list)
    num_values = round((max_value - min_value) / step_size) + 1
    return num_values, min_value

def read_elements(path:str, step_size:float, headers:list=None, degrees:bool=True) -> tuple:
    """
    Converts a CSV file into a grid of elements
    
    Parameters:
    * `path`:      The path to the CSV file
    * `step_size`: The step size of each element
    * `headers`:   The list of headers for the x coordinates, y coordinates,
                   grain IDs, phi_1, Phi, and phi_2 values
    * `degrees`:   Whether to store the orientations as degrees
    
    Returns the element grid
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
    
    # Read CSV and fill grid
    element_grid = get_void_element_grid(x_cells, y_cells)
    for row in rows:

        # Process data
        row_list = row.replace("\n", "").split(",")
        if "NaN" in row_list or "nan" in row_list:
            continue
        row_list = [float(val) for val in row_list]
        grain_id = round(row_list[index_grain_id])

        # Add to element grid
        x = round(float(row_list[index_x] - x_min) / step_size)
        y = round(float(row_list[index_y] - y_min) / step_size)
        element_grid[y][x] = Element(
            phi_1    = row_list[index_avg_phi_1],
            Phi      = row_list[index_avg_Phi],
            phi_2    = row_list[index_avg_phi_2],
            grain_id = grain_id,
            degrees  = degrees,
        )
    
    # Close file and return element grid
    file.close()
    return element_grid

def get_grain_ids(element_grid:list) -> list:
    """
    Gets the grain IDs of the element grid
    
    Parameters:
    * `element_grid`: A grid of elements
    
    Returns the list of grain IDs
    """
    grain_ids = list(set([element.get_grain_id() for element_list in element_grid for element in element_list]))
    if VOID_ELEMENT_ID in grain_ids:
        grain_ids.remove(VOID_ELEMENT_ID)
    return grain_ids

def get_areas(element_grid:list) -> dict:
    """
    Gets the areas for all the grains
    
    Parameters:
    * `element_grid`: A grid of elements
    
    Returns a dictionary that maps the grain IDs to their areas
    """
    grain_ids = get_grain_ids(element_grid)
    initial_areas = [0] * len(grain_ids)
    area_dict = dict(zip(grain_ids, initial_areas))
    for row in range(len(element_grid)):
        for col in range(len(element_grid[row])):
            grain_id = element_grid[row][col].get_grain_id()
            if grain_id == VOID_ELEMENT_ID:
                continue
            area_dict[grain_id] += 1
    return area_dict

def get_grain_positions(element_grid:list) -> dict:
    """
    Gets the x and y positions for all the grains
    
    Parameters:
    * `element_grid`: A grid of elements
    
    Returns a dictionary that maps the grain IDs
    to their x and y positions
    """
    
    # Initialise dictionary to store element positions
    grain_ids = get_grain_ids(element_grid)
    position_dict = {}
    for grain_id in grain_ids:
        position_dict[grain_id] = {"x": [], "y": []}
    
    # Add points to the element dictionary
    for row in range(len(element_grid)):
        for col in range(len(element_grid[row])):
            grain_id = element_grid[row][col].get_grain_id()
            if grain_id in [VOID_ELEMENT_ID]:
                continue
            position_dict[grain_id]["x"].append(col)
            position_dict[grain_id]["y"].append(row)
    
    # Returns the position mapping
    return position_dict

def get_centroids(element_grid:list) -> dict:
    """
    Gets the centroids for all the grains
    
    Parameters:
    * `element_grid`: A grid of elements
    
    Returns a dictionary that maps the grain IDs to their centroids
    """

    # Initialise
    grain_ids = get_grain_ids(element_grid)
    position_dict = get_grain_positions(element_grid)
    centroid_dict = {}
    
    # Calculate centroids for each grain
    for grain_id in grain_ids:
        x_mean = np.average(position_dict[grain_id]["x"])
        y_mean = np.average(position_dict[grain_id]["y"])
        centroid_dict[grain_id] = (x_mean, y_mean)
    
    # Return the centroids
    return centroid_dict
