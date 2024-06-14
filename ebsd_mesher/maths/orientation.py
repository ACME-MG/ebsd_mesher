"""
 Title:         Orientation
 Description:   Functions related to orientations
 References:    [1] https://www.researchgate.net/publication/324088567_Computing_Euler_angles_with_Bunge_convention_from_rotation_matrix
                [2] https://stackoverflow.com/questions/12374087/average-of-multiple-quaternions
 Author:        Janzen Choi

"""

# Libraries
import numpy as np, math, random
from scipy.spatial.transform import Rotation

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
    Determines the orientation matrix of a set of euler-bunge angles (rads);
    from Ref. [1]
    
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
    Determines the euler-bunge angles based on an orientation matrix (rads);
    from Ref. [1]

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

def euler_to_quat(euler:list) -> list:
    """
    Converts a set of euler-bunge angles (rads) into a quaternion

    Parameters:
    `euler`: The euler angle (rads)

    Returns the quaternion as a list
    """
    euler_array = np.array(euler)
    rotation = Rotation.from_euler("zxz", euler_array, degrees=False)
    quat = rotation.as_quat()
    return list(quat)

def quat_to_euler(quat:list) -> list:
    """
    Converts a quaternion into a set of euler-bunge angles (rads)

    Parameters:
    * `quaternion`: The quaternion

    Returns the euler angle as a list
    """
    quat_array = np.array(quat)
    rotation = Rotation.from_quat(quat_array)
    euler = rotation.as_euler("zxz", degrees=False)
    euler = [e+2*math.pi if e<0 else e for e in euler]
    return euler

def get_average_quat(quat_list:list) -> list:
    """
    Gets the average quaternion;
    from Ref. [2]

    Parameters:
    * `quat_list`: List of quaternions

    Returns the averaged quaternion
    """
    quat_array = np.array(quat_list)
    weights = np.array([1 for _ in range(len(quat_array))])
    average = np.linalg.eigh(np.einsum('ij,ik,i->...jk', quat_array, quat_array, weights))[1][:, -1]
    return list(average)

def get_average_euler(euler_list:list, degrees:bool=True) -> list:
    """
    Gets the average euler angle

    Parameters:
    * `euler_list`: List of euler-bunge angles (rads)
    * `degrees`:    Whether the euler angles are in degrees
    
    Returns the averaged euler angles
    """
    euler_list = [list(euler) for euler in euler_list]
    if degrees:
        euler_list = deg_to_rad(euler_list)
    quat_list = [euler_to_quat(euler) for euler in euler_list]
    average_quat = get_average_quat(quat_list)
    average_euler = quat_to_euler(average_quat)
    if degrees:
        average_euler = rad_to_deg(average_euler)
    return average_euler
