"""
 Title:         Plotter
 Description:   Helper functions for the plotter
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import numpy as np

def get_boundary(row:int, col:int, pixel_grid:list, step_size:float) -> tuple:
    """
    Gets coordinates for drawing the boundaries

    Parameters:
    * `row`:        The row of the pixel
    * `col`:        The column of the pixel
    * `pixel_grid`: The grid of pixels
    * `step_size`:  The step size

    Returns the x and y lists
    """

    # Initialise
    x_list, y_list = [], []
    x = get_coordinate(col, step_size)
    y = get_coordinate(row, step_size)
    
    # Check to add boundary on the right
    if col+1 >= len(pixel_grid[0]) or pixel_grid[row][col] != pixel_grid[row][col+1]:
        x_list += [x + step_size/2]*2 + [np.NaN]
        y_list += [y - step_size/2, y + step_size/2] + [np.NaN]

    # Check to add boundary on the left
    if col-1 < 0 or pixel_grid[row][col] != pixel_grid[row][col-1]:
        x_list += [x - step_size/2]*2 + [np.NaN]
        y_list += [y - step_size/2, y + step_size/2] + [np.NaN]

    # Check to add boundary on the top
    if row+1 >= len(pixel_grid) or pixel_grid[row][col] != pixel_grid[row+1][col]:
        x_list += [x - step_size/2, x + step_size/2] + [np.NaN]
        y_list += [y + step_size/2]*2 + [np.NaN]

    # Check to add boundary on the bottom
    if row-1 < 0 or pixel_grid[row][col] != pixel_grid[row-1][col]:
        x_list += [x - step_size/2, x + step_size/2] + [np.NaN]
        y_list += [y - step_size/2]*2 + [np.NaN]

    # Return the coordinates
    return x_list, y_list

def get_coordinate(index:int, step_size:float) -> float:
    """
    Converts an index into a coordinate value

    Parameters:
    * `index`:     The index
    * `step_size`: The step size

    Returns the coordinate value
    """
    return index*step_size + step_size/2

def get_positions(grain_id:int, pixel_grid:list) -> tuple:
    """
    Gets the positions of a grain in a pixel grid

    Parameters:
    * `grain_id`:   The ID of the grain
    * `pixel_grid`: The grid of pixels

    Returns the column and row positions
    """
    col_list, row_list = [], []
    for row in range(len(pixel_grid)):
        for col in range(len(pixel_grid[row])):
            if grain_id == pixel_grid[row][col]:
                col_list.append(col)
                row_list.append(row)
    return col_list, row_list

def save_plot(file_path:str) -> None:
    """
    Saves the plot and clears the figure

    Parameters:
    * `file_path`: The path to save the plot
    """
    plt.savefig(file_path)
    plt.cla()
    plt.clf()
    plt.close()
