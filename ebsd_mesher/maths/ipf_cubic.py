"""
 Title:         IPF Cubic
 Description:   For converting orientations from euler-bunge into RGB form for IPF
 Resource:      Adapted from https://mooseframework.inl.gov/docs/doxygen/modules/Euler2RGB_8h.html
                which is generalised for all/most structures
 Author:        Janzen Choi

"""

# Libraries
import math
from ebsd_mesher.maths.orientation import deg_to_rad, euler_to_matrix
from scipy.optimize import minimize

def get_cubic_symmetry_matrices():
    """
    Returns the cubic symmetry matrices
    """
    return [
        [[1,0,0], [0,1,0], [0,0,1]],
        [[0,0,1], [1,0,0], [0,1,0]],
        [[0,1,0], [0,0,1], [1,0,0]],
        [[0,-1,0], [0,0,-1], [1,0,0]],
        [[0,-1,0], [0,0,-1], [1,0,0]],
        [[0,1,0], [0,0,-1], [-1,0,0]],
        [[0,0,-1], [1,0,0], [0,-1,0]],
        [[0,0,-1], [-1,0,0], [0,1,0]],
        [[0,0,1], [-1,0,0], [0,-1,0]],
        [[-1,0,0], [0,1,0], [0,0,-1]],
        [[-1,0,0], [0,-1,0], [0,0,1]],
        [[1,0,0], [0,-1,0], [0,0,-1]],
        [[0,0,-1], [0,-1,0], [-1,0,0]],
        [[0,0,1], [0,-1,0], [1,0,0]],
        [[0,0,1], [0,1,0], [-1,0,0]],
        [[0,0,-1], [0,1,0], [1,0,0]],
        [[-1,0,0], [0,0,-1], [0,-1,0]],
        [[1,0,0], [0,0,-1], [0,1,0]],
        [[1,0,0], [0,0,1], [0,-1,0]],
        [[-1,0,0], [0,0,1], [0,1,0]],
        [[0,-1,0], [-1,0,0], [0,0,-1]],
        [[0,1,0], [-1,0,0], [0,0,-1]], # [[0,1,0], [-1,0,0], [0,0,1]]
        [[0,1,0], [1,0,0], [0,0,-1]],
        [[0,-1,0], [1,0,0], [0,0,1]],
    ]

def euler_to_rgb(phi_1:float, Phi:float, phi_2:float, ipf="x") -> tuple:
    """
    Converts orientation from euler-bunge to RGB (for cubic only)
    
    Parameters:
    * `phi_1`: The average phi_1 orientation of the grain
    * `Phi`:   The average Phi orientation of the grain
    * `phi_2`: The average phi_2 orientation of the grain
    * `ipf`:   The IPF colouring scheme
    
    Returns the RGB values
    """
    
    # Get IPF direction list
    if ipf == "x":
        ipf_list = [1,0,0]
    elif ipf == "y":
        ipf_list = [0,1,0]
    elif ipf == "z":
        ipf_list = [0,0,1]

    # Convert orientation
    phi_1 = deg_to_rad(phi_1)
    Phi   = deg_to_rad(Phi)
    phi_2 = deg_to_rad(phi_2)

    # Define auxiliary variables
    eta_min = deg_to_rad(0)
    eta_max = deg_to_rad(45)
    chi_min = deg_to_rad(0)
    chi_max = math.acos(1/math.sqrt(2+math.pow(math.tan(eta_max), 2)))

    # Assign black for out of RGB domain
    if phi_1 > 2*math.pi or Phi > math.pi or phi_2 > 2*math.pi:
        return 255, 255, 255
    
    # Get matrices
    orientation_matrix = euler_to_matrix([phi_1, Phi, phi_2])
    symmetry_matrices = get_cubic_symmetry_matrices()

    # Sort euler angles into SST
    for i in range(len(symmetry_matrices)):

        # Calculate temporary matrix
        temp_matrix = [[0,0,0],[0,0,0],[0,0,0]]
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    temp_matrix[j][k] += symmetry_matrices[i][j][l] * orientation_matrix[l][k]
        
        # Get multiple orientation matrix
        hkl = [0,0,0]
        for j in range(3):
            for k in range(3):
                hkl[j] += temp_matrix[j][k] * ipf_list[k]
        
        # Convert to spherical coordinates
        eta = abs(math.atan2(hkl[1], hkl[0]))
        chi = math.acos(abs(hkl[2]))

        # Check if eta and chi values are within SST, and keep searching otherwise
        if eta >= eta_min and eta < eta_max and chi >= chi_min and chi < chi_max:
            break

    # Calculate auxiliary variables
    chi_max_2 = math.acos(1/math.sqrt(2+math.pow(math.tan(eta),2)))
    eta_diff  = abs((eta-eta_min)/(eta_max-eta_min))

    # Calculate RGB colours
    red   = math.sqrt(abs(1-chi/chi_max_2))
    green = math.sqrt((1-eta_diff)*(chi/chi_max_2))
    blue  = math.sqrt(eta_diff*(chi/chi_max_2))

    # Normalise, round, and return
    max_rgb = max(red, green, blue)
    red     = round(red/max_rgb*255)
    green   = round(green/max_rgb*255)
    blue    = round(blue/max_rgb*255)
    return red, green, blue

def rgb_to_euler(red:int, green:int, blue:int, ipf:str="x"):
    """
    Converts RGB to euler-bunge orientation through optimisation;
    this function is mostly reliable, but can take some compute time
    
    Parameters:
    * `red`:   The amount of red
    * `blue`:  The amount of blue
    * `green`: The amount of green
    * `ipf`:   The IPF colouring scheme
    
    Returns the euler-bunge orientations
    """

    # Define the objective function
    def obj_func(x):
        phi_1, Phi, phi_2 = x
        red_, green_, blue_ = euler_to_rgb(phi_1, Phi, phi_2, ipf)
        return math.pow(red-red_, 2) + math.pow(green-green_, 2) + math.pow(blue-blue_, 2)
    
    # Get the euler angle
    euler = minimize(
        fun=obj_func,
        x0=(100,100,100),
        method="Powell",
    ).x

    # Wrap and return
    euler = [e % 360 for e in euler]
    return euler[0], euler[1], euler[2]
