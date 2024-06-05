"""
 Title:         Mesh Plotter
 Description:   Visualises a EBSD mesh using a plot
 Author:        Janzen Choi

"""

# Libraries
import pyvista as pv
import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from ebsd_mesher.maths.ipf_cubic import euler_to_rgb

# Mesh Plotter class
class MeshPlotter:

    def __init__(self, exodus_path:str, pixel_grid:list, grain_map:dict, step_size:float, figure_x:float):
        """
        Constructor for the mesh plotter

        Parameters:
        * `exodus_path`: The path to the mesh file
        * `pixel_grid`:  The grid of pixels
        * `grain_map`:   The dictionary of grain orientations
        * `step_size`:   The step size
        * `figure_x`:    The initial horizontal size of the figures
        """

        # Initialise internal variables
        self.exodus_path = exodus_path
        self.grain_map = grain_map
        
        # Initialise plot
        x_max = len(pixel_grid[0])*step_size
        y_max = len(pixel_grid)*step_size
        self.figure, self.axis = plt.subplots(figsize=(figure_x, y_max/x_max*figure_x))
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        self.axis.invert_yaxis()

    def plot_mesh(self, spn_to_exo:dict, ipf:str="x") -> None:
        """
        Creates a plot of the mesh

        Parameters:
        * `spn_to_exo`: Dictionary that maps the EBSD map to the mesh
        * `ipf`:        The IPF direction to plot the mesh
        """

        # Map mesh_id to ebsd_id
        exo_to_spn = dict(zip(spn_to_exo.values(), spn_to_exo.keys()))

        # Read mesh and iterate through grains
        mesh_info = pv.read(self.exodus_path)[0]
        for i, grain in enumerate(mesh_info):
            
            # Get IPF colour
            ebsd_id = exo_to_spn[i+1]
            euler = self.grain_map[ebsd_id].get_orientation()
            ipf_colour = euler_to_rgb(*euler, ipf=ipf)
            ipf_colour = [ipf/255 for ipf in ipf_colour]

            # Plot cells for the grain
            for cell_id in range(grain.n_cells):
                cell_coordinates = get_cell_coordinates(grain, cell_id)
                if cell_coordinates == []:
                    continue # only consider surface cells
                cell_coordinates = order_vertices(cell_coordinates)
                polygon = patches.Polygon(cell_coordinates, closed=True, fill=True, facecolor=ipf_colour, edgecolor="black")
                self.axis.add_patch(polygon)

def get_cell_coordinates(grain:pv.core.pointset.UnstructuredGrid, cell_id:int) -> list:
    """
    Gets the coordinates of a cell

    Parameters:
    * `grain`:   The grain object
    * `cell_id`: The cell ID

    Returns the list of coordinates for the surface corners of the cell
    """
    point_ids = grain.get_cell(cell_id).point_ids
    corners = [list(point) for point in grain.points[point_ids]]
    coordinates = [corner[:2] for corner in corners if corner[2] == 0]
    return coordinates

def on_segment(p:list, q:list, r:list): 
    """
    Check if point q lies on the line segment pr;
    only for colinear points

    Parameters:
    * `p`: The first point
    * `q`: The second point
    * `r`: The three point
    """
    return (
        (q[0] <= max(p[0], r[0])) and
        (q[0] >= min(p[0], r[0])) and
        (q[1] <= max(p[1], r[1])) and
        (q[1] >= min(p[1], r[1]))
    )
  
def orientation(p:list, q:list, r:list): 
    """
    Gets the orientations of three points

    Parameters:
    * `p`: The first point
    * `q`: The second point
    * `r`: The three point
    """
    val = (float(q[1] - p[1]) * (r[0] - q[0])) - (float(q[0] - p[0]) * (r[1] - q[1])) 
    if (val > 0): 
        return 1
    elif (val < 0): 
        return 2
    else: 
        return 0
  
def do_intersect(p_1:list, q_1:list, p_2:list, q_2:list) -> bool: 
    """
    Checks whether two segments intersect

    Parameters:
    * `p_1`: The first point of the first segment
    * `q_1`: The second point of the first segment
    * `p_2`: The first point of the second segment
    * `q_2`: The second point of the second segment

    Returns whether an intersection occurs or not
    """
    o_1 = orientation(p_1, q_1, p_2) 
    o_2 = orientation(p_1, q_1, q_2) 
    o_3 = orientation(p_2, q_2, p_1) 
    o_4 = orientation(p_2, q_2, q_1) 
    return (
        ((o_1 != o_2) and (o_3 != o_4)) or
        ((o_1 == 0) and on_segment(p_1, p_2, q_1)) or
        ((o_2 == 0) and on_segment(p_1, q_2, q_1)) or
        ((o_3 == 0) and on_segment(p_2, p_1, q_2)) or
        ((o_4 == 0) and on_segment(p_2, q_1, q_2))
    )

def order_vertices(unordered_vertices:list) -> list:
    """
    Orders a list of four vertices so that
    they form a quadrilateral

    Parameters:
    * `unordered_vertices`: List of unordered vertices

    Returns the ordered lists of vertices
    """

    # Iterate through all permutations of vertices    
    for permutation in itertools.permutations(unordered_vertices[1:]):
        vertices = [unordered_vertices[0]] + list(permutation)
        intersect = do_intersect(vertices[0], vertices[2], vertices[1], vertices[3])
        if intersect:
            return vertices
