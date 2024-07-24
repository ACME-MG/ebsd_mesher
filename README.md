# EBSD 1:1 Mesher
For creating a mesh that projects the microstructure from an EBSD map. Note that the mesh does not consider changing grain geometries along the z-axis. The following README file was last updated on the 24th of July, 2024.

# CSV format for EBSD data

The code relies on the EBSD data to be discretised and stored in a CSV file. This CSV file must have columns for the x-coordinate, y-coordinate, grain ID, phi_1, Phi, and phi_2 values for each of the discretised cells of the EBSD map. The phi_1, Phi, and phi_2 values should be in degrees for the for the code to work properly. The following is an example of how the CSV file should be formatted.

|  x   |  y   | grain_id | phi_1 |  Phi  | phi_2 |
|:----:|:----:|:--------:|:-----:|:-----:|:-----:|
| 5.0  | 5.5  |    1     | 60.1  | 88.2  | 54.2  |
| 5.0  | 6.0  |    1     | 57.3  | 91.0  | 53.7  |
| 5.0  | 6.5  |    1     | 57.2  | 89.5  | 54.1  |
| 5.0  | 7.0  |    2     | 87.4  | 29.4  | 94.1  |
| 5.0  | 7.5  |    2     | 87.5  | 29.3  | 93.9  |
| ...  | ...  |   ...    |  ...  |  ...  |  ...  |

# Importing `ebsd_mesher`

The user can add the `ebsd_mesher` package to the system path so that they can import the code anywhere on their machine. This involves simply adding the following line to the `~/.bashrc` file, as follows. Note that `<path_to_ebsd_mesher>` refers to the absolute path to the `ebsd_mesher` directory.
```
export PYTHONPATH=$PYTHONPATH:<path_to_ebsd_mesher>
```

Alternatively, the user can import `ebsd_mesher` dynamically depending on the relative path between the main `ebsd_mesher` directory and their script. For instance, if the user wants to run a script in the `ebsd_mesher/scripts/` directory, they can use the following.
```
import sys; sys.path += [".."]
import ebsd_mesher
```

# Interacting with `ebsd_mesher`

The following section provides a tutorial on using the `ebsd_mesher` code.

## Interface Class

The `Interface` class allows the user to interact with the `ebsd_mesher` code. To access the `Interface` class, the user must first import the `ebsd_mesher` package and initialise an `Interface` object. The following is an example of doing so from the `ebsd_mesher/scripts/` directory.
```py
import sys; sys.path += [".."]
from ebsd_mesher import Interface
itf = Interface()
```

The `Interface` class contains several optional arguments. These include:
* `title`: This optional argument appends a title in front of the directory in which the output files will be placed in. The default value for this argument is an empty string.
* `output_path`: This optional argument defines the relative path to the output directory, which tells the script where to place the output files. The default value for this arguemnt is `"./results"`.
* `verbose`: This optional argument tells the script whether to display any information about the actions of the `Interface` class in the terminal. The default value for this argument is `True`, meaning that the information will be displayed in the terminal.
* `output_here`: This optional argument tells the script whether to just place the output files in the same directory as the script. The default value for this is `False`. Note that when the user sets the argument to `True`, the `title` and `output_path` values will not have any effect.

The implementation of the `Interface` class can be accessed via `ebsd_mesher/ebsd_mesher/interface.py`. The next sections provide descriptions of the available functions, their available arguments, and how to use them. Note that additional descriptions of the `Interface` functions can also be accessed by hovering your cursor over the functions. However, this functionality is only supported by some IDEs (e.g., Visual Studio Code).

## Defining headers of EBSD file (`define_headers`)

The `define_headers` function defines the headers for the EBSD CSV files. If the function isn't called, then default values for the headers (i.e., "x", "y", "grain_id", "phi_1", "Phi", "phi_2") will be used.

Parameters:
* `x`:        Header for the x-coordinate.
* `y`:        Header for the y-coordinate.
* `grain_id`: Header for the grain ID.
* `phi_1`:    Header for the phi_1 values.
* `Phi`:      Header for the Phi values.
* `phi_2`:    Header for the phi_2 values.

## Adding EBSD map and map grains (`import_ebsd`)

The `import_ebsd` function adds an EBSD map and conducts the mapping.

Parameters:
* `ebsd_path`: Path to the EBSD file as a CSV file.
* `step_size`: Step size between coordinates.
* `degrees`:   Whether the orientations in the file are in degrees or radians. The default value is `True`, which imports the orientations as degrees.

## Getting the bounds of the EBSD map (`get_bounds`)

The `get_bounds` function returns the bounds of the imported EBSD map. The returned bounds is in the form of a dictionary, with the following format:
```
{"x_min": x_min, "x_max": x_max, "y_min": y_min, "y_max": y_max}
```

Note that the bounds ignore the 'void' elements, where no information (e.g., coordinates, orientations) have been provided in the imported CSV file for the EBSD data.

## Exporting the bounds (`export_bounds`)

The `export_bounds` exports the bounds to a file named, `bounds.csv`.

## Redefining the domain (`redefine_domain`)

The `redefine_domain` function redefines the domain of the EBSD map. Particularly useful for when the user wants to mesh a certain portion of the EBSD map.

Parameters:
* `x_min`: The minimum x value of the domain.
* `x_max`: The maximum x value of the domain.
* `y_min`: The minimum y value of the domain.
* `y_max`: The maximum y value of the domain.

## Decreasing the resolution of the EBSD map (`decrease_resolution`)

The `decrease_resolution` decreases the resolution of the imported EBSD map. Note that the function essentially changes the step size of the imported EBSD data.

Parameters:
* `factor`: The factor to reduce the resolution by. Using a factor > 1 will decrease the resolution while using a factor < 1 will increase the resolution. 

## Cleaning the EBSD map (`clean`)

The `clean` function cleans the elements of the EBSD map to improve the mesh quality.

Parameters:
* `iterations`: The number of times to conduct the cleaning. Note that a small number of iterations is recommended as to not drastically change the grain geometry of the microstructure.

## Smoothing the EBSD map (`smooth`)

The `smooth` function smooths the elements of the EBSD map to improve the mesh quality.

Parameters:
* `iterations`: The number of times to conduct the smoothing. Note that a small number of iterations is recommended as to not drastically change the grain geometry of the microstructure.

## Filling the voids in the EBSD map (`fill`)

The `fill` function fills in the 'void' elements of the EBSD map. This function is useful when incomplete information (e.g., coordinates, orientations) has been provided in the imported CSV file for the EBSD data.

Parameters:
* `iterations`: The number of times to conduct the cleaning.

## Removing small grains (`remove_grains`)

The `remove_grains` removes small grains in the EBSD map.

Parameters:
* `threshold`: The threshold grain area to remove the grains.

## Adding 'grips' to the mesh (`add_grips`)

The `add_grips` adds elements to the left and right of the microstructure.

Parameters:
* `num_elements`: The number of elements (in the x-direction) to use in the grips.

## Plotting the EBSD maps (`plot_ebsd`)

The `plot_ebsd` function plots the EBSD maps.

Parameters:
* `ebsd_path`:        The path to save the plotted EBSD map. The default value is `"ebsd"`, which saves the plot to the defined results path with the name, `ebsd.png`.
* `ipf`:              The IPF colour ("x", "y", "z"). The default value is "x".
* `figure_x`:         The initial horizontal size of the figures. The default value is 10.
* `include_id`:       Whether to include IDs in the EBSD maps; the user can define a dictionary for custom plotting settings. The default value is False.
* `include_boundary`: Whether to include IDs in the EBSD maps; the user can define a dictionary for custom plotting settings. The default value is False.
* `id_list`:          The IDs of the grains to plot the grain IDs and boundaries. The IDs should be defined based on the first grain map. If some of the defined IDs are not mappable, they will be ignored. The default value is None, which adds for all mappable grains.

## Generating the mesh (`mesh`)

The `mesh` function generates a mesh from the imported EBSD map.

Parameters:
* `psculpt_path`:   The path to PSculpt.
* `z_elements`:     The number of elements in the the mesh's thickness (z-axis). The default value is 1, which creates a mesh with a thickness of 1 element.
* `num_processors`: The number of processors to use to create the mesh. The default value is 1, which tells the script to generate the mesh using only 1 processor.

## Exporting the grain orientations (`export_grains`)

The `export_grains` exports the orientation and area of each grain.

Parameters:
* `degrees`: Whether to save the orientations as degrees or radians. The default value is `True`, which saves the orientations as degrees.

## Exporting the element orientations (`export_elements`)

The `export_elements` exports the orientation and grain ID of each element.

Parameters:
* `degrees`: Whether to save the orientations as degrees or radians. The default value is `True`, which saves the orientations as degrees.

## Fixing grip interfaces (`fix_grip_interfaces`)

The `fix_grip_interfaces` function fices the interfaces between the grips and the microstructure.

Parameters:
* `grip_length`:  The desired length of the grip.
* `micro_length`: The desired length of the microstructure.

## Scaling the mesh (`scale_mesh`)

The `scale_mesh` function scales the mesh to a certain length along a certain direction.

Parameters:
* `length`:    The length to scale the mesh to.
* `direction`: The direction to conduct the scaling ("x", "y", "z").

## Plotting the mesh (`plot_mesh`)

The `plot_mesh` function plots the generated mesh. Note that the mesh must have already been generated for the function to work.

Parameters:
* `mesh_path`: Path to generate the plot. The default value is `"mesh"`, which saves the plot as `mesh.png` in the defined results path.
* `ipf`:       The IPF colour ("x", "y", "z"). The default value is "x".
* `figure_x`:  The initial horizontal size of the figures. The default value is 10.

