"""
 Title:         Mesher
 Description:   For converting the element grid into a mesh;
                Does not use adaptive meshing because it's broken
 Author:        Janzen Choi

"""

# Libraries
import subprocess, os

# Input file
INPUT_FILE_CONTENT = """
BEGIN SCULPT
    
    # Dimensions
    nelx = {x_cells}
    nely = {y_cells}
    nelz = {z_cells}
    scale = {step_size}
    
    # Fixed mesh improvement
    smooth = 3
    defeature = 1
    pillow_curves = true
    pillow_boundaries = true
    micro_shave = true
    
    # Variable mesh improvement
    opt_threshold = 0.7
    pillow_curve_layers = 3
    pillow_curve_thresh = 0.3

    # Solver
    laplacian_iters = 5
    max_opt_iters = 50
    
    # Output
    input_spn = {spn_file}
    exodus_file = {exodus_file}

END SCULPT
"""

def coarse_mesh(psculpt_path:str, thickness:int, num_processors:int, element_grid:list,
                step_size:float, input_path:str, spn_path:str, exodus_path:str) -> None:
    """
    Generates a mesh based on an SPN file
    
    Parameters:
    * `psculpt_path`:   The path to PSculpt 
    * `thickness`:      The thickness of the mesh (in voxels)
    * `num_processors`: The number of processors to use to create the mesh
    * `element_grid`:   The grid of elements
    * `step_size`:      The size of the elements
    * `input_path`:     The path to the input file
    * `spn_path`:       The path to the spn file
    * `exodus_path`:    The path to the mesh file
    """

    # Write SPN file (x = gauge, y = height, z = thickness)
    with open(spn_path, "w+") as fh:
        for i in range(len(element_grid[0])):  # x
            for j in range(len(element_grid)): # y
                for _ in range(thickness):     # z
                    fh.write(f"{element_grid[j][i].get_grain_id()} ")

    # Create input file
    with open(input_path, "w+", newline="") as fh:
        fh.write(INPUT_FILE_CONTENT.format(
            x_cells     = len(element_grid[0]),
            y_cells     = len(element_grid),
            z_cells     = thickness,
            step_size   = step_size,
            spn_file    = spn_path,
            exodus_file = exodus_path,
        ))

    # Run mesh command
    command = f"mpiexec -n {num_processors} {psculpt_path} -j {num_processors} -i {input_path}"
    subprocess.run([command], shell=True, check=True)
    os.rename(f"{exodus_path}.e.1.0", f"{exodus_path}")
