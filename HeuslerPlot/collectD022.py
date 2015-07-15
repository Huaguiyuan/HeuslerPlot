import argparse
import os
import shutil

def _copy_scf(scf_Es_path, collect_path, csv_data):
    '''Copy SCF data from scf_Es_path/e#/system/(normalmoment or bigmoment)/
    to collect_path/e#/system/.
    '''
    moment_dirs = {'NORMAL':'normalmoment', "BIG":'bigmoment'}
    for system_name, system_attrs in csv_data.items():
        moment_type = system_attrs[0]
        enum = system_attrs[1]
        moment_dirname = moment_dirs[moment_type]
        estr = "e{}".format(str(enum))

        scf_src_path = os.path.join(scf_Es_path, estr, system_name, moment_dirname)
        scf_dest_path = os.path.join(collect_path, estr, system_name)
        os.makedirs(scf_dest_path, exist_ok=True)

        filenames = ["OUTCAR", "OSZICAR"]
        _copy_files(filenames, scf_src_path, scf_dest_path)

def _copy_bands(bands_path, collect_path, csv_data):
    '''Copy bands data from bands_path/system/ to collect_path/e#/system/BANDS/.
    '''
    for system_name, system_attrs in csv_data.items():
        enum = system_attrs[1]
        estr = "e{}".format(str(enum))

        bands_src_path = os.path.join(bands_path, system_name)
        bands_dest_path = os.path.join(collect_path, estr, system_name, "BANDS")
        os.makedirs(bands_dest_path, exist_ok=True)

        filenames = ["EIGENVAL", "KPOINTS"]
        _copy_files(filenames, bands_src_path, bands_dest_path)

def _copy_files(filenames, src_dir_path, dest_dir_path):
    for name in filenames:
        file_src_path = os.path.join(src_dir_path, name)
        file_dest_path = os.path.join(dest_dir_path, name)
        shutil.copy(file_src_path, file_dest_path)

def _parse_D022_csv(csv_path):
    '''Return a dictionary with keys given by system names ("Fe2MnIn", etc.)
    and values of the form (moment_type, enum) with moment_type = 'NORMAL' or
    'BIG' and enum = number of valence electrons.
    '''
    fp = open(csv_path, 'r')
    lines = fp.readlines()
    fp.close()

    result = {}
    for line in lines:
        lsp = line.split(',')
        system_name = lsp[0]
        moment_type = lsp[1]
        enum = int(lsp[2])
        result[system_name] = (moment_type, enum)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect D022 scf and bands data into format usable by plot.py.")
    parser.add_argument('D022_scf_Es_path', type=str, help="Path to directory with D022 SCF files (e#/sys/ = SCF).")
    parser.add_argument('D022_bands_path', type=str, help="Path to directory with D022 band files (sys/ = bands).")
    parser.add_argument('collect_path', type=str, help="Path to directory where plotting data will be assembled.")
    parser.add_argument('D022_csv_path', type=str, help="Path to CSV file describing D022s (need to know 'NORMAL' or 'BIG' moment).")
    args = parser.parse_args()

    csv_data = _parse_D022_csv(args.D022_csv_path)
    _copy_scf(args.D022_scf_Es_path, args.collect_path, csv_data)
    _copy_bands(args.D022_bands_path, args.collect_path, csv_data)
