"""
 Title:         Neighbour
 Description:   For determining neighbours of grains
 Author:        Janzen Choi

"""

def get_neighbour_dict(pixel_grid:list) -> dict:
    """
    Gets neighbouring grain information;
    does not repeat (if a->b then b-/>a)

    Parameters:
    * `pixel_grid`: A grid of pixels

    Returns the dictionary of grain neighbours
    """

    # Initialise
    x_size = len(pixel_grid[0])
    y_size = len(pixel_grid)
    neighbour_dict = {}

    # Iterate through grid of pixels
    for row in range(y_size):
        for col in range(x_size):

            # Initialise grain information if uninitialised
            grain_id = pixel_grid[row][col]
            if not grain_id in neighbour_dict.keys():
                neighbour_dict[grain_id] = []

            # Identify neighbouring grains
            neighbours = get_neighbours(col, row, x_size, y_size)
            for n in neighbours:
                n_id = pixel_grid[n[1]][n[0]]
                if (grain_id != n_id and not n_id in neighbour_dict[grain_id] and
                    (not n_id in neighbour_dict.keys() or not grain_id in neighbour_dict[n_id])):
                    neighbour_dict[grain_id].append(n_id)
    
    # Return dictionary
    return neighbour_dict

def get_neighbours(x:float, y:float, x_size:int, y_size:int) -> list:
    """
    Gets the neighbouring indices of a pixel
    
    Parameters:
    * `x`:      The x coordinate
    * `y`:      The y coordinate
    * `x_size`: The maximum x value 
    * `y_size`: The maximum y value
    
    Returns a list of the neighbouring coordinates 
    """
    neighbours = [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]
    neighbours = [
        neighbour for neighbour in neighbours
        if neighbour[0] >= 0 and neighbour[0] < x_size
        and neighbour[1] >= 0 and neighbour[1] < y_size
    ]
    return neighbours

def get_common_neighbours(pixel_grid:list, x:float, y:float, x_size:int, y_size:int) -> list:
    """
    Gets the neighbouring indices of a pixel that have the
    same values in the pixel grid
    
    Parameters:
    * `pixel_grid`: A list of pixels
    * `x`:          The x coordinate
    * `y`:          The y coordinate
    * `x_size`:     The maximum x value 
    * `y_size`:     The maximum y value
    
    Returns a list of the common neighbouring coordinates 
    """
    current_value = pixel_grid[y][x]
    neighbours = get_neighbours(x, y, x_size, y_size)
    common_neighbours = []
    for neighbour in neighbours:
        neighbour_value = pixel_grid[neighbour[1]][neighbour[0]]
        if neighbour_value == current_value:
            common_neighbours.append(neighbour)
    return common_neighbours

def get_all_neighbours(x_list:list, y_list:list, x_size:int, y_size:int):
    """
    Gets the neighbouring indices of a group of pixels
    
    Parameters:
    * `x_list`: The list of x coordinates
    * `y_list`: The list of y coordinates
    * `x_size`: The maximum x value 
    * `y_size`: The maximum y value
    
    Returns a list of all the neighbouring coordinates 
    """
    
    # Gets all the neighbours
    all_neighbours = []
    for i in range(len(x_list)):
        neighbours = get_neighbours(x_list[i], y_list[i], x_size, y_size)
        all_neighbours += neighbours
    
    # Remove duplicates and neighbours in the group
    all_neighbours = list(set(all_neighbours))
    group = [(x_list[i],y_list[i]) for i in range(len(x_list))]
    all_neighbours = [neighbour for neighbour in all_neighbours if not neighbour in group]

    # Return
    return all_neighbours
