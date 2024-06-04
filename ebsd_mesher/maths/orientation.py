"""
 Title:         Orientation
 Description:   Functions related to orientations
 References:    https://www.researchgate.net/publication/324088567_Computing_Euler_angles_with_Bunge_convention_from_rotation_matrix
 Author:        Janzen Choi

"""

# Libraries
import numpy as np, math, random

def get_matrix_product(matrix_1:list, matrix_2:list) -> list:
    """
    Performs a 3x3 matrix multiplication

    Parameters:
    * `matrix_1`: The first matrix
    * `matrix_2`: The second matrix

    Returns the final matrix
    """
    result = [[0,0,0], [0,0,0], [0,0,0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                result[i][j] += matrix_1[i][k] * matrix_2[k][j]
    return result

def get_inverted(matrix:list) -> list:
    """
    Inverts a matrix

    Parameters:
    * `matrix`: The matrix to be inverted

    Returns the inverted matrix as a list of lists
    """
    matrix = np.array(matrix)
    inverted = [list(i) for i in np.linalg.inv(matrix)]
    return inverted

def euler_to_matrix(euler:list) -> list:
    """
    Determines the orientation matrix of a set of euler-bunge angles (rads)
    
    Parameters:
    * `euler`: The euler angle in euler-bunge form

    Returns the 3x3 orientation matrix as a list of lists
    """
    om_11 = math.cos(euler[0])*math.cos(euler[2]) - math.sin(euler[0])*math.sin(euler[2])*math.cos(euler[1])
    om_12 = math.sin(euler[0])*math.cos(euler[2]) + math.cos(euler[0])*math.sin(euler[2])*math.cos(euler[1])
    om_13 = math.sin(euler[2])*math.sin(euler[1])
    om_21 = -math.cos(euler[0])*math.sin(euler[2]) - math.sin(euler[0])*math.cos(euler[2])*math.cos(euler[1])
    om_22 = -math.sin(euler[0])*math.sin(euler[2]) + math.cos(euler[0])*math.cos(euler[2])*math.cos(euler[1])
    om_23 = math.cos(euler[2])*math.sin(euler[1])
    om_31 = math.sin(euler[0])*math.sin(euler[1])
    om_32 = -math.cos(euler[0])*math.sin(euler[1])
    om_33 = math.cos(euler[1])
    om = [[om_11, om_12, om_13],
          [om_21, om_22, om_23],
          [om_31, om_32, om_33]]
    return om

def matrix_to_euler(matrix:list) -> list:
    """
    Determines the euler-bunge angles based on an orientation matrix (rads)

    Parameters:
    * `matrix`: The orientation matrix

    Returns the euler angles
    """
    Phi = math.acos(matrix[2][2])
    if Phi == 0:
        phi_1 = math.atan2(-matrix[1][0], matrix[0][0])
        phi_2 = 0
    elif Phi == math.pi:
        phi_1 = math.atan2(matrix[1][0], matrix[0][0])
        phi_2 = 0
    else:
        phi_1 = math.atan2(matrix[2][0], -matrix[2][1])
        phi_2 = math.atan2(matrix[0][2], matrix[1][2])
    phi_1 = phi_1 + 2*math.pi if phi_1 < 0 else phi_1
    phi_2 = phi_2 + 2*math.pi if phi_2 < 0 else phi_2
    return [phi_1, Phi, phi_2]

def rad_to_deg(radians:float) -> float:
    """
    Converts radians to degrees

    Parameters:
    * `radians`: The radians to be converted

    Returns the converted degrees
    """
    if isinstance(radians, list):
        return [rad_to_deg(r) for r in radians]
    return radians * 180 / math.pi

def deg_to_rad(degrees:float) -> float:
    """
    Converts degrees to radians

    Parameters:
    * `degrees`: The degrees to be converted

    Returns the converted radians
    """
    if isinstance(degrees, list):
        return [deg_to_rad(d) for d in degrees]
    return degrees * math.pi / 180

def random_euler() -> list:
    """
    Generates a set of (uniformly) random euler-bunge angles;
    """
    phi_1 = random.random() * 360.
    Phi = 180. * math.acos(2 * random.random() - 1) / np.pi
    phi_2 = random.random() * 360.
    return [phi_1, Phi, phi_2]

def random_quat():
    """
    Generates a (uniformly) random quaternion
    """
    u = [random.uniform(0, 1) for _ in range(3)]
    x = math.sqrt(1 - u[0]) * math.sin(2 * math.pi * u[1])
    y = math.sqrt(1 - u[0]) * math.cos(2 * math.pi * u[1])
    z = math.sqrt(u[0]) * math.sin(2 * math.pi * u[2])
    w = math.sqrt(u[0]) * math.cos(2 * math.pi * u[2])
    return [x, y, z, w]

def euler_to_quat(phi_1:float, Phi:float, phi_2:float) -> float:
    """
    Converts a set of euler-bunge angles into a quaternion (rads)

    Parameters:
    `phi_1`: The first euler angle
    `Phi`:   The second euler angle
    `phi_2`: The third euler angle

    Returns the quaternion as a list
    """
    cy = math.cos(phi_2 * 0.5)
    sy = math.sin(phi_2 * 0.5)
    cp = math.cos(Phi * 0.5)
    sp = math.sin(Phi * 0.5)
    cr = math.cos(phi_1 * 0.5)
    sr = math.sin(phi_1 * 0.5)
    x = sr * cp * cy - cr * sp * sy
    y = cr * sp * cy + sr * cp * sy
    z = cr * cp * sy - sr * sp * cy
    w = cr * cp * cy + sr * sp * sy
    return [x, y, z, w]

def quat_to_euler(x:float, y:float, z:float, w:float) -> list:
    """
    Converts a quaternion into a set of euler-bunge angles (rads)

    Parameters:
    * `x`: The first quaternion value
    * `y`: The second quaternion value
    * `z`: The third quaternion value
    * `w`: The fourth quaternion value

    Returns the euler angle as a list
    """
    phi_1 = math.atan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y))
    Phi   = math.asin(max([min([2 * (w * y - z * x), 1]), -1]))
    phi_2 = math.atan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z))
    return [phi_1, Phi, phi_2]
