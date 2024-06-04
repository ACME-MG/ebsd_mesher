"""
 Title:         Coincidence Site Lattice (CSL)
 Description:   For generating euler angles that satisfy CSL criteria
 Referennces:   http://pajarito.materials.cmu.edu/lectures/L14-CSL_Theory_GBE-17Mar16.pdf
                https://github.com/heprom/pymicro/blob/master/pymicro/crystal/microstructure.py
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from ebsd_mesher.maths.orientation import random_euler, deg_to_rad, euler_to_matrix, matrix_to_euler, get_matrix_product

# Dictionary of CSLs
CSL_DICT = {
    "3":    {"mori": 60.00, "euler": [45, 70.53, 45]},
    "5":    {"mori": 36.86, "euler": [0, 90, 36.86]},
    "7":    {"mori": 38.21, "euler": [26.56, 73.4, 63.44]},
    "9":    {"mori": 38.94, "euler": [26.56, 83.62, 26.56]},
    "11":   {"mori": 50.47, "euler": [33.68, 79.53, 33.68]},
    "13a":  {"mori": 22.62, "euler": [0, 90, 22.62]},
    "13b":  {"mori": 27.79, "euler": [18.43, 76.66, 71.57]},
    "15":   {"mori": 48.19, "euler": [19.65, 82.33, 42.27]},
    "17a":  {"mori": 28.07, "euler": [0, 90, 28.07]},
    "17b":  {"mori": 61.90, "euler": [45, 86.63, 45]},
    "19a":  {"mori": 26.53, "euler": [18.44, 89.68, 18.44]},
    "19b":  {"mori": 46.80, "euler": [33.69, 71.59, 56.31]},
    "21a":  {"mori": 21.78, "euler": [14.03, 79.02, 75.97]},
    "21b":  {"mori": 44.41, "euler": [22.83, 79.02, 50.91]},
    "23":   {"mori": 40.45, "euler": [15.25, 82.51, 52.13]},
    "25a":  {"mori": 16.26, "euler": [0, 90, 16.26]},
    "25b":  {"mori": 51.68, "euler": [36.87, 90, 53.13]},
    "27a":  {"mori": 31.59, "euler": [21.8, 85.75, 21.8]},
    "27b":  {"mori": 35.43, "euler": [15.07, 85.75, 31.33]},
    "29a":  {"mori": 43.60, "euler": [0, 90, 43.6]},
    "29b":  {"mori": 46.40, "euler": [33.69, 84.06, 56.31]},
    "31a":  {"mori": 17.90, "euler": [11.31, 80.72, 78.69]},
    "31b":  {"mori": 52.20, "euler": [27.41, 78.84, 43.66]},
    # "33a":  {"mori": 20.10, "euler": [12.34, 83.04, 58.73]}, # issue
    # "33b":  {"mori": 33.60, "euler": [37.51, 76.84, 37.51]}, # issue
    "33c":  {"mori": 59.00, "euler": [38.66, 75.97, 38.66]},
    "35a":  {"mori": 34.00, "euler": [16.86, 80.13, 60.46]},
    "35b":  {"mori": 43.20, "euler": [30.96, 88.36, 59.04]},
}

def get_symmetry_matrices(type:str="cubic") -> list:
    """
    Returns the symmetry matrices

    Parameters:
    * `type`: The crystal structure type

    Returns a list of the symmetry matrices given the type
    """
    if type == "cubic":
        return get_cubic_symmetry_matrices()
    elif type == "hexagonal":
        return get_hexagonal_symmetry_matrices()
    elif type == "tetrahedral":
        return get_tetrahedral_symmetry_matrices()

def get_cubic_symmetry_matrices() -> list:
    """
    Returns a list of cubic symmetry matrices
    """
    return [
        [[1,0,0],  [0,1,0],  [0,0,1]],
        [[0,0,1],  [1,0,0],  [0,1,0]],
        [[0,1,0],  [0,0,1],  [1,0,0]],
        [[0,-1,0], [0,0,1],  [-1,0,0]],
        [[0,-1,0], [0,0,-1], [1,0,0]],
        [[0,1,0],  [0,0,-1], [-1,0,0]],
        [[0,0,-1], [1,0,0],  [0,-1,0]],
        [[0,0,-1], [-1,0,0], [0,1,0]],
        [[0,0,1],  [-1,0,0], [0,-1,0]],
        [[-1,0,0], [0,1,0],  [0,0,-1]],
        [[-1,0,0], [0,-1,0], [0,0,1]],
        [[1,0,0],  [0,-1,0], [0,0,-1]],
        [[0,0,-1], [0,-1,0], [-1,0,0]],
        [[0,0,1],  [0,-1,0], [1,0,0]],
        [[0,0,1],  [0,1,0],  [-1,0,0]],
        [[0,0,-1], [0,1,0],  [1,0,0]],
        [[-1,0,0], [0,0,-1], [0,-1,0]],
        [[1,0,0],  [0,0,-1], [0,1,0]],
        [[1,0,0],  [0,0,1],  [0,-1,0]],
        [[-1,0,0], [0,0,1],  [0,1,0]],
        [[0,-1,0], [-1,0,0], [0,0,-1]],
        [[0,1,0],  [-1,0,0], [0,0,1]],
        [[0,1,0],  [1,0,0],  [0,0,-1]],
        [[0,-1,0], [1,0,0],  [0,0,1]],
    ]

def get_hexagonal_symmetry_matrices() -> list:
    """
    Returns the hexagonal symmetry matrices
    """
    a = (3 ** 0.5) / 2
    return [
        [[1,0,0],     [0,1,0],     [0,0,1]],
        [[-0.5,a,0],  [-a,-0.5,0], [0,0,1]],
        [[-0.5,-a,0], [a,-0.5,0],  [0,0,1]],
        [[0.5,a,0],   [-a,0.5,0],  [0,0,1]],
        [[-1,0,0],    [0,-1,0],    [0,0,1]],
        [[0.5,-a,0],  [a,0.5,0],   [0,0,1]],
        [[-0.5,-a,0], [-a,0.5,0],  [0,0,-1]],
        [[1,0,0],     [0,-1,0],    [0,0,-1]],
        [[-0.5,a,0],  [a,0.5,0],   [0,0,-1]],
        [[0.5,a,0],   [a,-0.5,0],  [0,0,-1]],
        [[-1,0,0],    [0,1,0],     [0,0,-1]],
        [[0.5,-a,0],  [-a,-0.5,0], [0,0,-1]],
    ]

def get_tetrahedral_symmetry_matrices() -> list:
    """
    Returns the tetrahedral symmetry matrices
    """
    return [
        [[1,0,0],  [0,1,0],  [0,0,1]],
        [[-1,0,0], [0,1,0],  [0,0,-1]],
        [[1,0,0],  [0,-1,0], [0,0,-1]],
        [[-1,0,0], [0,-1,0], [0,0,1]],
        [[0,1,0],  [-1,0,0], [0,0,1]],
        [[0,-1,0], [1,0,0],  [0,0,1]],
        [[0,1,0],  [1,0,0],  [0,0,-1]],
        [[0,-1,0], [-1,0,0], [0,0,-1]],
    ]

def get_csl_euler_angles(csl_sigma:str, euler_1:list=None) -> list:
    """
    Generates two sets of euler angles that conform to CSL3

    Parameters:
    * `csl_sigma`: The sigma value of the CSL
    * `euler_1`:   The optional first set of euler angles

    Returns a list of two euler angles
    """

    # Generate a set of random euler angles if none specified
    euler_1 = euler_1 if euler_1 != None else random_euler()
    
    # Specify rotational offset
    euler_offset = deg_to_rad(CSL_DICT[csl_sigma]["euler"])

    # Determine second set of euler angles
    matrix_1 = euler_to_matrix(euler_1)
    matrix_offset = euler_to_matrix(euler_offset)
    matrix_2 = get_matrix_product(matrix_offset, matrix_1)
    euler_2 = matrix_to_euler(matrix_2)

    # Return
    return [euler_1, euler_2]

def get_misorientations(euler_1:list, euler_2:list, type:str) -> list:
    """
    Determines the misorientations of two sets of euler angles (rads)

    Parameters:
    * `euler_1`: The first euler angle
    * `euler_2`: The second euler angle
    * `type`:    The crystal structure type
    
    Returns a list of the misorientation angles from the symmetry matrices
    """

    # Get orientation and symmetry matrices
    orientation_1 = np.array(euler_to_matrix(euler_1))
    orientation_2 = np.array(euler_to_matrix(euler_2))
    symmetries = np.array(get_symmetry_matrices(type))

    # Iterate through symmetry matrices
    misorientation_list = []
    for symmetry_1 in symmetries:
        operator_1 = np.dot(symmetry_1, orientation_1)
        for symmetry_2 in symmetries:
            operator_2 = np.dot(symmetry_2, orientation_2)
            delta = np.dot(operator_2, operator_1.T)
            cw = 0.5 * (delta.trace() - 1)
            if cw > 1. and cw - 1. < 10 * np.finfo("float32").eps:
                cw = 1.
            misorientation = np.arccos(cw)
            misorientation_list.append(misorientation)
    return misorientation_list

def get_disorientation(euler_1:list, euler_2:list, type:str) -> float:
    """
    Determines the minimal misorientation of two sets of euler angles (rads)

    Parameters:
    * `euler_1`: The first euler angle (rads)
    * `euler_2`: The second euler angle (rads)
    * `type`:    The crystal structure type
    
    Returns the disorientation angle
    """
    misorientations = get_misorientations(euler_1, euler_2, type)
    misorientations += get_misorientations(euler_2, euler_1, type)
    return min(misorientations)

# Testing
# from orientation import rad_to_deg
# euler_pairs = [
#     [[-89,23,38], [-89,-38,23]], # 49.03
#     [[97,-19,-4], [-97,-19,-4]], # 135.09
# ]
# for euler_pair in euler_pairs:
#     euler_rad_1 = deg_to_rad(euler_pair[0])
#     euler_rad_2 = deg_to_rad(euler_pair[1])
#     print(rad_to_deg(get_disorientation(euler_rad_1, euler_rad_2, "cubic")))
