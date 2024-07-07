"""
 Title:         Improver
 Description:   Manipulates the values in the element grid
 Author:        Janzen Choi

"""

# Libraries
from random import randint
from copy import deepcopy
from ebsd_mesher.mapper.gridder import VOID_ELEMENT_ID, get_void_element, get_void_element_grid, get_areas, get_grain_ids
from ebsd_mesher.mapper.neighbour import get_neighbours

def clean_element_grid(element_grid:list) -> list:
    """
    Cleans the element grid by replacing stray void / live elements
    
    Parameters:
    * `element_grid`: The uncleaned grid of elements
    
    Returns the cleaned element grid
    """

    # Dimensions of the element grid
    x_size = len(element_grid[0])
    y_size = len(element_grid)
    
    # Iterate through each element
    for row in range(y_size):
        for col in range(x_size):
            element = element_grid[row][col]

            # Evaluate neighbouring elements
            neighbours = get_neighbours(col, row, x_size, y_size)
            void_neighbours = [n for n in neighbours if element_grid[n[1]][n[0]].get_grain_id() == VOID_ELEMENT_ID]
            num_void = len(void_neighbours)
            
            # If half (or less) of the neighbours are void, then fill a void element
            if element.get_grain_id() == VOID_ELEMENT_ID and num_void <= len(neighbours) / 2:
                copy = neighbours[randint(0, len(neighbours) - 1)]
                element_grid[row][col] = deepcopy(element_grid[copy[1]][copy[0]])

            # If more than half of the neighbours are void, then remove a live element
            if element.get_grain_id() != VOID_ELEMENT_ID and num_void > len(neighbours) / 2:
                element_grid[row][col] = get_void_element()

    # Return cleaned element grid
    return element_grid

def smoothen_edges(element_grid:list) -> list:
    """
    Smoothen the edges of grains by merging elements
    
    Parameters:
    * `element_grid`: The unsmoothed grid of elements
    
    Returns the smoothed element grid
    """

    # Dimensions of the element grid
    x_size = len(element_grid[0])
    y_size = len(element_grid)
    
    # Iterate through each element
    for row in range(y_size):
        for col in range(x_size):

            # Evaluate neighbouring elements
            neighbours = get_neighbours(col, row, x_size, y_size)
            foreign_neighbours = []
            for n in neighbours:
                if element_grid[n[1]][n[0]].get_grain_id() != element_grid[row][col].get_grain_id():
                    foreign_neighbours.append(n)

            # If there are more than 2 foreign neighbours, get absorbed
            if len(foreign_neighbours) > 2:
                foreign = foreign_neighbours[randint(0, len(foreign_neighbours) - 1)]
                element_grid[row][col] = deepcopy(element_grid[foreign[1]][foreign[0]])

    # Return cleaned element grid
    return element_grid

def pad_edges(element_grid:list) -> list:
    """
    Pads the element grid by replicating unvoided elements
    
    Parameters:
    * `element_grid`: The unpadded grid of elements
    
    Returns the padded element grid
    """
    
    # Dimensions of the element grid
    x_size = len(element_grid[0])
    y_size = len(element_grid)
    
    # Replicate it
    padded_element_grid = get_void_element_grid(x_size, y_size)

    # Iterate through each element
    for row in range(y_size):
        for col in range(x_size):
            element = element_grid[row][col]

            # If live, copy and skip
            if element.get_grain_id() != VOID_ELEMENT_ID:
                padded_element_grid[row][col] = deepcopy(element)
                continue

            # Get live neighbouring elements
            neighbours = get_neighbours(col, row, x_size, y_size)
            live_neighbours = [n for n in neighbours if element_grid[n[1]][n[0]].get_grain_id() != VOID_ELEMENT_ID]

            # If there is a live neighbour, then fill this void element
            if len(live_neighbours) > 0:
                copy = live_neighbours[randint(0, len(live_neighbours) - 1)]
                padded_element_grid[row][col] = deepcopy(element_grid[copy[1]][copy[0]])
    
    # Return padded element grid
    return padded_element_grid

def remove_small_grains(element_grid:list, threshold:int) -> list:
    """
    Removes small grains
    
    Parameters:
    * `element_grid`: The unremoved grid of elements
    * `threshold`:  The grain size threshold to start removing grains
    
    Returns the element grid without the small grains
    """
    
    # Get grain IDs that do not pass the threshold
    grain_ids = get_grain_ids(element_grid)
    area_dict = get_areas(element_grid)
    to_remove = [grain_id for grain_id in grain_ids if area_dict[grain_id] <= threshold]
    
    # Remove small grains
    for row in range(len(element_grid)):
        for col in range(len(element_grid[row])):
            if element_grid[row][col].get_grain_id() in to_remove:
                element_grid[row][col] = get_void_element()
    
    # Print removed grains and return the new element grid
    print(f"\n    Removed {len(to_remove)}/{len(grain_ids)} grain(s)\n")
    return element_grid
