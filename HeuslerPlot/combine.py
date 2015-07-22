from argparse import ArgumentParser
import sys
import os
import shutil

def _compositions_from_csv(csv_path):
    fp = open(csv_path, 'r')
    lines = fp.readlines()
    fp.close()

    compositions = []
    for i, line in enumerate(lines):
        # Ignore first line (header)
        if i == 0:
            continue
        # Composition = first column
        comp = line.strip().split(',')[0]
        # First column may be empty.
        if comp == "":
            continue
        # Non-empty; add to list.
        compositions.append(comp)

    return compositions

if __name__ == "__main__":
    parser = ArgumentParser(description="Take L21 and D022 plots as specified in csv files and combine in one directory.")
    parser.add_argument('--L21_plots_dir', type=str, default=None,
            help="Directory containing L21 plots")
    parser.add_argument('--D022_plots_dir', type=str, default=None,
            help="Directory containing D022 plots")
    parser.add_argument('--L21_csv_path', type=str, default=None,
            help="CSV file describing L21 systems")
    parser.add_argument('--D022_csv_path', type=str, default=None,
            help="CSV file describing D022 systems")
    parser.add_argument('--combined_dir', type=str, default=None,
            help="Directory to put combined L21/D022 plots")
    parser.add_argument('--img_root_dir', type=str, default=None,
            help="If specified, ignore other arguments and determine paths as img_root_dir/(L21, D022, D0_22 Full.csv, FullHeusler.csv, combined)")
    args = parser.parse_args()

    L21_plots_dir = args.L21_plots_dir
    D022_plots_dir = args.D022_plots_dir
    L21_csv_path = args.L21_csv_path
    D022_csv_path = args.D022_csv_path
    combined_dir = args.combined_dir
    if args.img_root_dir != None:
        L21_plots_dir = os.path.join(args.img_root_dir, "L21")
        D022_plots_dir = os.path.join(args.img_root_dir, "D022")
        L21_csv_path = os.path.join(args.img_root_dir, "FullHeusler.csv")
        D022_csv_path = os.path.join(args.img_root_dir, "D0_22 Full.csv")
        combined_dir = os.path.join(args.img_root_dir, "combined")
    elif None in [L21_plots_dir, D022_plots_dir, L21_csv_path, D022_csv_path, combined_dir]:
        print("Error: if img_root_dir not specified, must specify all other arguments.")
        sys.exit(2)

    L21_compositions = _compositions_from_csv(L21_csv_path)
    D022_compositions = _compositions_from_csv(D022_csv_path)

    # Take all D022 systems.
    for comp in D022_compositions:
        image_name = comp + "_D022.png"
        image_path = os.path.join(D022_plots_dir, image_name)
        new_image_path = os.path.join(combined_dir, image_name)
        shutil.copy(image_path, new_image_path)

    # Take L21s which are not listed as D022.
    for comp in L21_compositions:
        if comp in D022_compositions:
            continue
        # If we get here, system is not D022.
        image_name = comp + "_L21.png"
        image_path = os.path.join(L21_plots_dir, image_name)
        new_image_path = os.path.join(combined_dir, image_name)
        shutil.copy(image_path, new_image_path)
