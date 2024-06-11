"""
 Title:         Exodus
 Description:   For manipulating exodus files
 Author:        Janzen Choi

"""

# Libraries
import math
import pyvista as pv
import netCDF4 as nc
from ebsd_mesher.mapper.gridder import get_grain_positions

def map_spn_to_exo(exodus_path:str, spn_path:str, spn_size:tuple) -> tuple:
    """
    Gets the grain IDs of exodus grains from the SPN file
    
    Parameters:
    * `exodus_path`: The path to the exodus file
    * `spn_path`:    The path to the SPN file
    * `spn_size`:    The size of the SPN file as a tuple (x, y, z)
    
    Returns a mapping of the SPN to exodus IDs and the confidence of the mapping;
    the mapping is in the form of a dictionary with keys "exo_id", "confidence"
    """

    # Reads the contents of the exodus file
    exo_grains = pv.read(exodus_path)[0]
    bounds = exo_grains.bounds
    exo_bounds = [{"min": bounds[2*i], "max": bounds[2*i+1], "range": bounds[2*i+1] - bounds[2*i]} for i in range(3)]

    # Read the contents of the SPN file
    with open(spn_path, "r") as fh:
        voxel_string = " ".join(fh.readlines())
        voxel_list = [int(voxel) for voxel in voxel_string.split(" ") if voxel != ""]

    # Initialise dictionaries
    spn_to_exo = {}
    confidence_list = []

    # Iterate through the exodus grains
    for i in range(exo_grains.n_blocks):
        
        # Get grain elements
        exo_grain = exo_grains[i]
        elements = exo_grain.cell_centers().points
        elements = [list(element) for element in elements]

        # Get the grain ids based on element coordinates
        id_list = []
        for element in elements:
            if math.nan in element:
                continue
            pos = [math.floor((element[j] - exo_bounds[j]["min"]) / exo_bounds[j]["range"] * spn_size[j]) for j in range(3)]
            grain_id_index = pos[0] * spn_size[1] * spn_size[2] + pos[1] * spn_size[2] + pos[2]
            grain_id = voxel_list[grain_id_index]
            id_list.append(grain_id)
        
        # Add exodus grain id
        mode = max(set(id_list), key=id_list.count)
        freq = id_list.count(mode)
        total = len(id_list)
        confidence = round(freq / total * 100, 2)

        # Update dictionaries
        exo_id = int(str(exo_grains.get_block_name(i)).split(" ")[-1])
        spn_to_exo[mode] = exo_id
        confidence_list.append(confidence)

    # Return
    return spn_to_exo, confidence_list

def get_element_info(exodus_path:str, element_grid:list, spn_to_exo:dict, step_size:float) -> list:
    """
    Gets a list of elements; the elements are ordered by the grains
    then elements within the grains
    
    Parameters:
    * `exodus_path`:  The path to the exodus file
    * `element_grid`: The grid of elements
    * `spn_to_exo`:   Mapping from spn id to exodus id
    * `step_size`:    The size of each element

    Returns an ordered list element objects
    """
    
    # Print message (it takes a while)
    print("\n=================================================")
    print("Mapping EBSD elements to mesh elements...")
    print("=================================================\n")

    # Initialise
    exo_to_spn = dict(zip(spn_to_exo.values(), spn_to_exo.keys()))
    position_dict = get_grain_positions(element_grid)
    get_distance = lambda a, b : math.sqrt(math.pow(a[0]-b[0],2) + math.pow(a[1]-b[1],2))
    element_list = []

    # Read grains and iterate through them
    mesh = pv.read(exodus_path)[0]
    for i in range(mesh.n_blocks):

        # Get grain and grain information
        grain = mesh[i]
        exo_id = int(str(mesh.get_block_name(i)).split(" ")[-1])
        grain_id = exo_to_spn[exo_id]

        # Get centroids of the elements in the grain
        elements = grain.cell_centers().points
        centroid_list = [list(element) for element in elements]

        # Iterate through centroids
        for centroid in centroid_list:

            # Get the positions of all elements in the grain
            positions = position_dict[grain_id]
            positions = [(positions["x"][i], positions["y"][i]) for i in range(len(positions["x"]))]
            positions_scaled = [(position[0]*step_size, position[1]*step_size) for position in positions]

            # Find element closest to centroid
            distances = [get_distance(position, centroid) for position in positions_scaled]
            min_index = distances.index(min(distances))
            opt_position = positions[min_index]
            element = element_grid[opt_position[1]][opt_position[0]]
            element_list.append(element)
    
    # Return the list of element objects
    return element_list

def renumber_grains(exodus_path:str) -> None:
    """
    Renumbers the grain IDs in an exodus file so that they
    start from 1 and end at num_grains

    Parameters:
    * `exodus_path`: Path to the exodus file
    """
    ds = nc.Dataset(exodus_path, mode="r+")
    block_ids = ds.variables["eb_prop1"]
    block_ids[:] = range(1, len(block_ids)+1)
    ds.close()

def get_exodus_dimension(exodus_path:str, direction:str="z") -> float:
    """
    Gets the dimension of the mesh in a certaind irection

    Parameters:
    * `exodus_path`: Path to the exodus file
    * `direction`:   The direction to scale ("x", "y", "z")

    Returns the dimension
    """
    ds = nc.Dataset(exodus_path, mode="r")
    dimension = max(ds.variables[f"coord{direction}"][:])
    ds.close()
    return dimension

def scale_exodus_mesh(exodus_path:str, factor:float, direction:str="z") -> None:
    """
    Scales the coordinates of the cell in an exodus file

    Parameters:
    * `exodus_path`: Path to the exodus file
    * `factor`:      The amount to scale
    * `direction`:   The direction to scale ("x", "y", "z")
    """
    ds = nc.Dataset(exodus_path, mode="r+")
    coords = ds.variables[f"coord{direction}"]
    coords[:] = coords[:] * factor
    ds.close()

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
    grain_id_list = [int(str(mesh.get_block_name(i)).split(" ")[-1])
                     for i in range(len(mesh)) if not i+1 in exclude]
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

def scale_gripped_microstructure(exodus_path:str, l_grip_id:int, r_grip_id:int,
                                 grip_length:float, micro_length:float) -> None:
    """
    Scales the grips and microstructure to desired dimensions
    
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
