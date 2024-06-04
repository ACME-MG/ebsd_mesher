"""
 Title:         Exporter
 Description:   Exports stuff
 Author:        Janzen Choi

"""

# Libraries
import math, pyvista as pv

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
        spn_to_exo[mode] = i+1
        confidence_list.append(confidence)

    # Return
    return spn_to_exo, confidence_list
