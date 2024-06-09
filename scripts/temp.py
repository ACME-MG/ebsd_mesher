import netCDF4 as nc
import pyvista as pv
import numpy as np
import shutil

def get_points(grain:pv.core.pointset.UnstructuredGrid) -> list:
    """
    Gets the points within a grain

    Parameters:
    * `grain`: The grain

    Returns the list of points
    """
    return [list(p) for p in grain.points]

def get_all_points(mesh:pv.core.composite.MultiBlock, exclude:list=None) -> list:
    """
    Gets the coordinates of all the points

    Parameters:
    * `mesh`:    The multiblock of grains
    * `exclude`: The list of grain IDs to avoid
    * `x_lower`: The lower bound on the horizontal axis to extract the points
    * `y_upper`: The upper bound on the horizontal axis to extract the points
    
    Returns the list of foreign points
    """
    exclude = [] if exclude == None else exclude
    grain_id_list = [i+1 for i in range(len(mesh)) if not i+1 in exclude]
    point_list = []
    for grain_id in grain_id_list:
        grain = mesh[grain_id-1]
        point_list += get_points(grain)
    return point_list

def get_border_points(curr_points:list, foreign_points:list) -> list:
    """
    Gets the coordinates of the points neighbouring a certain grain

    Parameters:
    * `curr_points`:    List of points in the grain
    * `foreign_points`: List of foreign points near the grain
    
    Returns the list of bordering points
    """
    border_point_list = []
    for point in curr_points:
        if point in foreign_points and not point in border_point_list:
            border_point_list.append(point)
    return border_point_list

def get_x_lists(exodus_path:str, grip_id:int) -> tuple:
    """
    Gets useful lists of x coordinates
    
    Parameters:
    * `exodus_path`: The path to the exodus mesh
    * `grip_id`:     The ID corresponding to the grip

    Returns the x coordinates in the grip, the other grains,
    and the grip interface
    """

    # Gets the points from the grip, other grains, and grip interface
    pv_mesh       = pv.read(exodus_path)[0]
    grip_points   = get_points(pv_mesh[grip_id-1])
    other_points  = get_all_points(pv_mesh, [grip_id])
    border_points = get_border_points(grip_points, other_points)

    # Extracts x coordinates from the points
    grip_x_list   = list(set([point[0] for point in grip_points]))
    border_x_list = list(set([point[0] for point in border_points]))

    # Returns extracted x coordinates
    return grip_x_list, border_x_list

def straighten_interface(exodus_path:str, grip_id:int, left:bool) -> None:
    """
    Straightens the interface between the grip and the microstructure
    
    Parameters:
    * `exodus_path`: The path to the exodus mesh
    * `grip_id`:     The ID corresponding to the grip
    * `direction`:   Whether the grip is on the left
    """

    # Get x coordinates and calculate bounds
    grip_x_list, border_x_list = get_x_lists(exodus_path, grip_id)
    max_x   = max(border_x_list)
    range_x = abs(max_x-min(border_x_list))
    
    # Fix interface
    nc_mesh = nc.Dataset(exodus_path, mode="a")
    for i in range(len(nc_mesh.variables["coordx"])):
        x_coord = nc_mesh.variables["coordx"][i]
        if x_coord in border_x_list:
            nc_mesh.variables["coordx"][i] = max_x
        elif left and not x_coord in grip_x_list:
            nc_mesh.variables["coordx"][i] += range_x
        elif not left and x_coord in grip_x_list:
            nc_mesh.variables["coordx"][i] += range_x
    nc_mesh.close()

def scale_mesh(exodus_path:str, l_grip_id:int, r_grip_id:int,
               grip_length:float, micro_length:float) -> None:
    """
    Scales the mesh to desired dimensions
    
    Parameters:
    * `exodus_path`:  The path to the exodus mesh
    * `l_grip_id`:    The ID of the left grip
    * `r_grip_id`:    The ID of the right grip
    * `grip_length`:  The desired length of the grip
    * `micro_length`: The desired length of the microstructure
    """

    # Gets the x coordinates in the grips
    pv_mesh       = pv.read(exodus_path)[0]
    l_grip_x_list = list(set([point[0] for point in get_points(pv_mesh[l_grip_id-1])]))
    r_grip_x_list = list(set([point[0] for point in get_points(pv_mesh[r_grip_id-1])]))
    
    # Gets the dimensions of the unscaled grips and microstructure
    l_grip_max = max(l_grip_x_list)
    r_grip_min = min(r_grip_x_list) # max of microstructure
    r_grip_max = max(r_grip_x_list)
    l_grip_mesh_length = l_grip_max
    micro_mesh_length  = r_grip_min - l_grip_max
    r_grip_mesh_length = r_grip_max - r_grip_min

    # Define scaling functions
    scale_l_grip = lambda x : x*grip_length/l_grip_mesh_length
    scale_micro  = lambda x : (x-l_grip_max)*micro_length/micro_mesh_length+grip_length
    scale_r_grip = lambda x : (x-r_grip_min)*grip_length/r_grip_mesh_length+grip_length+micro_length

    # Conduct scaling
    nc_mesh = nc.Dataset(exodus_path, mode="a")
    for i in range(len(nc_mesh.variables["coordx"])):
        x_coord = nc_mesh.variables["coordx"][i]
        if x_coord in l_grip_x_list:
            nc_mesh.variables["coordx"][i] = scale_l_grip(x_coord)
        elif x_coord in r_grip_x_list:
            nc_mesh.variables["coordx"][i] = scale_r_grip(x_coord)
        else:
            nc_mesh.variables["coordx"][i] = scale_micro(x_coord)
    nc_mesh.close()

old_path = "./results/240607210036_simple_premesh/mesh.e"
new_path = "./mesh.e"
shutil.copy(old_path, new_path)

# 72, 73
l_grip_id = 72
r_grip_id = 73
grip_length = 600
micro_length = 2500

straighten_interface(new_path, l_grip_id, left=True)
straighten_interface(new_path, r_grip_id, left=False)
scale_mesh(new_path, l_grip_id, r_grip_id, grip_length, micro_length)
