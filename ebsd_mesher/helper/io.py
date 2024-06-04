"""
 Title:         Input / Output
 Description:   Reading / writing functions
 Author:        Janzen Choi

"""

# Libraries
import math, os
import pandas as pd
from ebsd_mesher.helper.general import round_sf

def get_file_path_writable(file_path:str, extension:str):
    """
    Appends a number after a path if it is not writable

    Parameters:
    * `file_path`: Path to file without the extension
    * `extension`: The extension for the file
    """
    new_file_path = f"{file_path}.{extension}"
    if not os.path.exists(new_file_path):
        return new_file_path
    index = 1
    while True:
        try:
            with open(new_file_path, 'a'):
                return new_file_path
        except IOError:
            new_file_path = f"{file_path} ({index}).{extension}"
            index += 1

def get_file_path_exists(file_path:str, extension:str):
    """
    Appends a number after a path if it exists

    Parameters:
    * `file_path`: Path to file without the extension
    * `extension`: The extension for the file
    """
    new_file_path = f"{file_path}.{extension}"
    index = 1
    while os.path.exists(new_file_path):
        new_file_path = f"{file_path} ({index}).{extension}"
        index += 1
    return new_file_path

def csv_to_dict(csv_path:str, delimeter:str=",") -> dict:
    """
    Converts a CSV file into a dictionary
    
    Parameters:
    * `csv_path`:  The path to the CSV file
    * `delimeter`: The separating character
    
    Returns the dictionary
    """

    # Read all data from CSV (assume that file is not too big)
    csv_fh = open(csv_path, "r", encoding="utf-8-sig")
    csv_lines = csv_fh.readlines()
    csv_fh.close()

    # Initialisation for conversion
    csv_dict = {}
    headers = csv_lines[0].replace("\n", "").split(delimeter)
    csv_lines = csv_lines[1:]
    for header in headers:
        csv_dict[header] = []

    # Start conversion to dict
    for csv_line in csv_lines:
        csv_line_list = csv_line.replace("\n", "").split(delimeter)
        for i in range(len(headers)):
            value = csv_line_list[i]
            if value == "":
                continue
            try:
                value = float(value)
            except:
                pass
            csv_dict[headers[i]].append(value)
    
    # Convert single item lists to items and things multi-item lists
    for header in headers:
        if len(csv_dict[header]) == 1:
            csv_dict[header] = csv_dict[header][0]
    
    # Return
    return csv_dict

def dict_to_csv(data_dict:dict, csv_path:str, add_header:bool=True) -> None:
    """
    Converts a dictionary to a CSV file
    
    Parameters:
    * `data_dict`: The dictionary to be converted
    * `csv_path`:  The path that the CSV file will be written to
    * `header`:    Whether to include the header or not
    """
    
    # Extract headers and turn all values into lists
    headers = data_dict.keys()
    for header in headers:
        if not isinstance(data_dict[header], list):
            data_dict[header] = [data_dict[header]]
    
    # Open CSV file and write headers
    csv_fh = open(csv_path, "w+")
    if add_header:
        csv_fh.write(",".join(headers) + "\n")
    
    # Write data and close
    max_list_size = max([len(data_dict[header]) for header in headers])
    for i in range(max_list_size):
        row_list = [str(data_dict[header][i]) if i < len(data_dict[header]) else "" for header in headers]
        row_str = ",".join(row_list)
        csv_fh.write(row_str + "\n")
    csv_fh.close()

def read_excel(excel_path:str, sheet:str, column:int) -> list:
    """
    Reads an excel file

    Parameters:
    * `excel_path`: The path to the excel file
    * `sheet`:      The name of the sheet to read from
    * `column`:     The column index

    Returns a list of values corresponding to that column
    """
    data_frame = pd.read_excel(excel_path, sheet_name=sheet)
    data_list = list(data_frame.iloc[:,column])
    data_list = list(filter(lambda x: not math.isnan(x), data_list))
    data_list = [round_sf(data, 8) for data in data_list]
    return data_list

def safe_mkdir(dir_path:str) -> None:
    """
    For safely making a directory

    Parameters:
    * `dir_path`: The path to the directory
    """
    try:
        os.mkdir(dir_path)
    except FileExistsError:
        pass
