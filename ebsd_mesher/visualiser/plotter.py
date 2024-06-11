"""
 Title:         Plotter
 Description:   Helper functions for the plotter
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import numpy as np

def get_boundary(row:int, col:int, element_grid:list, step_size:float) -> tuple:
    """
    Gets coordinates for drawing the boundaries

    Parameters:
    * `row`:          The row of the element
    * `col`:          The column of the element
    * `element_grid`: The grid of elements
    * `step_size`:    The step size

    Returns the x and y lists
    """

    # Initialise
    x_list, y_list = [], []
    x = get_coordinate(col, step_size)
    y = get_coordinate(row, step_size)
    
    # Check to add boundary on the right
    if col+1 >= len(element_grid[0]) or element_grid[row][col].get_grain_id() != element_grid[row][col+1].get_grain_id():
        x_list += [x + step_size/2]*2 + [np.NaN]
        y_list += [y - step_size/2, y + step_size/2] + [np.NaN]

    # Check to add boundary on the left
    if col-1 < 0 or element_grid[row][col].get_grain_id() != element_grid[row][col-1].get_grain_id():
        x_list += [x - step_size/2]*2 + [np.NaN]
        y_list += [y - step_size/2, y + step_size/2] + [np.NaN]

    # Check to add boundary on the top
    if row+1 >= len(element_grid) or element_grid[row][col].get_grain_id() != element_grid[row+1][col].get_grain_id():
        x_list += [x - step_size/2, x + step_size/2] + [np.NaN]
        y_list += [y + step_size/2]*2 + [np.NaN]

    # Check to add boundary on the bottom
    if row-1 < 0 or element_grid[row][col].get_grain_id() != element_grid[row-1][col].get_grain_id():
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

def get_positions(grain_id:int, element_grid:list) -> tuple:
    """
    Gets the positions of a grain in a element grid

    Parameters:
    * `grain_id`:     The ID of the grain
    * `element_grid`: The grid of elements

    Returns the column and row positions
    """
    col_list, row_list = [], []
    for row in range(len(element_grid)):
        for col in range(len(element_grid[row])):
            if grain_id == element_grid[row][col].get_grain_id():
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
