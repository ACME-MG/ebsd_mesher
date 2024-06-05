"""
 Title:         Exodus
 Description:   For manipulating exodus files
 Author:        Janzen Choi

"""

# Libraries
import math, pyvista as pv
import netCDF4 as nc

def map_spn_to_exo(exo_path:str, spn_path:str, spn_size:tuple) -> tuple:
    """
    Gets the grain IDs of exodus grains from the SPN file
    
    Parameters:
    * `exo_path`: The path to the exodus file
    * `spn_path`: The path to the SPN file
    * `spn_size`: The size of the SPN file as a tuple (x, y, z)
    
    Returns a mapping of the SPN to exodus IDs and the confidence of the mapping;
    the mapping is in the form of a dictionary with keys "exo_id", "confidence"
    """

    # Reads the contents of the exodus file
    exo_grains = pv.read(exo_path)[0]
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
