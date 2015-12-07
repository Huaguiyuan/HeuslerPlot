import os
import re
from pathlib import Path

def FindEs(dir_path, subdir_path=None, bands_name=None):
    '''Look for subdirectories of dir_path which have names of the form
    e[0-9][0-9] (such as 'e18', 'e29', etc.). For each such subdirectory, 
    search it for subdirectories with the appropriate band files (as given
    by FindBands).

    Return a dictionary with keys given by system names, as specified by the
    names of the subdirectories of e##. The corresponding values are the
    same as those returned by FindBands, giving another dictionary with keys
    "outcar_path", "oszicar_path", "eigenval_path", and "kpoints_path" giving
    the appropriate file paths.
    '''
    result = {}
    # Enumerate the subdirectories of dir_path.
    p = Path(dir_path)
    sub_paths = [x for x in p.iterdir() if x.is_dir()]
    # Iterate over the subdirectories if dir_path.
    # Look for those with the form "e[0-9][0-9]".
    pattern = re.compile("e[0-9][0-9]")
    for sub_path in sub_paths:
        # Reduce sub_path: /path/to/bands/e18 --> e18.
        sub_name = os.path.basename(str(sub_path))

        if pattern.fullmatch(sub_name) != None:
            # sub_path has the correct form; search its subdirectories
            # and add their data.
            this_e_result = FindBands(str(sub_path), subdir_path, bands_name)
            for k, v in this_e_result.items():
                result[k] = v

    return result

def FindBands(dir_path, additional_subdir_path=None, bands_name=None):
    '''For each subdirectory sub_path contained in dir_path, look for
    sub_path/OUTCAR, sub_path/OSZICAR, sub_path/bands_name/EIGENVAL, and
    sub_path/bands_name/KPOINTS.

    If additional_subdir_path != None, add it to sub_path.

    Return a dictionary with keys given by system names, as specified by the
    names of the subdirectories of dir_path. The corresponding values give
    another dictionary with keys "outcar_path", "oszicar_path", "eigenval_path",
    and "kpoints_path" giving the appropriate file paths.
    '''
    if bands_name == None:
        bands_name = "BANDS"

    result = {}
    # Enumerate the subdirectories of dir_path.
    p = Path(dir_path)
    sub_paths = [x for x in p.iterdir() if x.is_dir()]
    # Iterate over the subdirectories if dir_path.
    for sub_path in sub_paths:
        sub_name = os.path.basename(str(sub_path))
        if additional_subdir_path != None:
            sub_path = sub_path / additional_subdir_path
        # Assemble paths to required data files in this subdirectory.
        subdir_result = {}
        subdir_result["outcar_path"] = str(sub_path / "OUTCAR")
        subdir_result["oszicar_path"] = str(sub_path / "OSZICAR")
        subdir_result["eigenval_path"] = str(sub_path / bands_name / "EIGENVAL")
        subdir_result["kpoints_path"] = str(sub_path / bands_name / "KPOINTS")
        # Check if all of these files actually exist.
        subdir_ok = True
        for data_file in subdir_result.values():
            if not os.path.exists(data_file):
                subdir_ok = False
        # If all required files exist, add this sub_path to the result.
        if subdir_ok:
            result[sub_name] = subdir_result

    return result

def FindCollected(dir_path, structure_from_filename = False):
    '''Look for files contained in dir_path with the form:
        FILETYPE_structuretype_compound
    where FILETYPE is KPOINTS, EIGENVAL, OSZICAR, and OUTCAR.

    Return a dictionary with keys given by system names, as specified by the
    names of the subdirectories of e##. The corresponding values are the
    same as those returned by FindBands, giving another dictionary with keys
    "outcar_path", "oszicar_path", "eigenval_path", and "kpoints_path" giving
    the appropriate file paths.
    '''
    filetype_keys = {"KPOINTS": "kpoints_path", "EIGENVAL": "eigenval_path",
            "OSZICAR": "oszicar_path", "OUTCAR": "outcar_path"}

    result = {}
    # Enumerate the files in dir_path.
    p = Path(dir_path)
    sub_paths = [x for x in p.iterdir() if not x.is_dir()]
    # Iterate over the files in dir_path.
    for sub_path in sub_paths:
        basename = os.path.basename(str(sub_path))
        name_split = basename.split('_')
        filetype = name_split[0]
        key = filetype_keys[filetype]

        structure = name_split[1]
        compound = name_split[-1]

        if compound not in result:
            result[compound] = {}
        result[compound][key] = str(sub_path)

        if structure_from_filename:
            result[compound]["structure"] = structure

    return result
