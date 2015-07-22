import argparse
import os
from pathlib import Path
import subprocess

def _in_dir_with_suffix(dir_path, suffix):
    '''Return a list of paths to all files in dir_path with filenames
    ending in the specified suffix.
    '''
    p = Path(dir_path)
    sub_paths = p.iterdir()
    suffix_paths = []
    for sub_path in sub_paths:
        if sub_path.is_dir():
            continue
        sub_path_str = str(sub_path)
        prefix_len = len(sub_path_str) - len(suffix)
        sub_path_suffix = sub_path_str[prefix_len:]
        if sub_path_suffix == suffix:
            suffix_paths.append(str(sub_path))
    return suffix_paths

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compress all *.png files in a directory.")
    parser.add_argument('--dir_path', type=str, default=None,
            help="Path to directory with data to compress (use current directory if not specified).")
    parser.add_argument('--compressed_name', type=str, default="",
            help="Add this to filename after compression instead of replacing original file.")
    parser.add_argument('--colors', type=int, default=4,
            help="Number of colors to use in output image (more colors = better antialiasing).")
    args = parser.parse_args()

    dir_path = args.dir_path
    if dir_path == None:
        dir_path = os.getcwd()

    png_paths = _in_dir_with_suffix(dir_path, ".png")
    for path in png_paths:
        prefix_len = len(path) - len(".png")
        out_path_prefix = path[:prefix_len]
        out_path = "{}{}.png".format(out_path_prefix, args.compressed_name)

        subprocess.call(["pngquant", str(args.colors), path, "--ext", "-compress.png", "--force"])
        subprocess.call(["mv", "{}-compress.png".format(out_path_prefix), out_path])
