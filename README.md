# Installation

Instructions for installation on Debian derivatives (Debian/Ubuntu/Mint/etc.).

Requires setuptools, numpy, and matplotlib:

    sudo apt-get install python3-setuptools python3-dev python3-numpy python3-matplotlib python3-tk zlib1g zlib1g-dev

compress.py requires libpng >= 1.6:

    cd ~
    wget ftp://ftp.simplesystems.org/pub/libpng/png/src/libpng16/libpng-1.6.21.tar.gz
    tar -xvzf libpng-1.6.17.tar.gz
    cd libpng-1.6.17.tar.gz
    ./configure
    make check
    sudo make install

Add /usr/local/lib to `LD_LIBRARY_PATH` -- add the following to ~/.bashrc and then restart bash:

    if [ -z "$LD_LIBRARY_PATH" ]; then
        export LD_LIBRARY_PATH=/usr/local/lib
    else
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
    fi

compress.py also requires [pngquant](https://pngquant.org/):

    cd ~
    git clone https://github.com/pornel/pngquant.git
    cd pngquant
    ./configure
    sudo make install

Go back to HeuslerPlot directory and install using setup.py:

    python3 setup.py install --user

To have changes to the source reflected immediately:

    python3 setup.py develop --user

# Usage: band plotting

To plot bands in a directory at "dirpath" with a structure dirpath/e29/Co2MnSi/BANDS (where there may be many e## and e##/systemname subdirectories) use (where `structure_type` is set appropriately):

    python3 plot.py --searchEs --logo --structure_type "L21" "dirpath"

If the bands subdirectory has a different name (such as dirpath/e29/Co2MnSi/bands instead of .../BANDS) but the directory structure is otherwise the same as above (this is the structure for half-Heusler cubic compounds); also exclude compounds containing `bands_t` subdirectory indicating the compound is tetragonal:

    python3 plot.py --searchEs --logo --structure_type "C1b" --bands_name "bands" --exclude_if_present "bands_t" "dirpath"

If there is an additional subdirectory containing tetragonal scf and bands runs, e.g. `dirpath/e29/Co2MnSi/bands_t` gives the tetragonal SCF and `dirpath/e29/Co2MnSi/bands_t/bands` gives the tetragonal bands (this is the structure for half-Heusler tetragonal compounds):

    python3 plot.py --searchEs --logo --structure_type "tetragonal" --bands_subdir_path "bands_t" --bands_name "bands" "dirpath"

To plot bands where files have been collected in a directory with filenames of the form `EIGENVAL_Xa_Fe2CoAs` (from bands calculation), `OUTCAR_Xa_Fe2CoAs` (from SCF), `OSZICAR_Xa_Fe2CoAs` (from SCF) (inverse Heuslers have been collected like this):

    python3 plot.py --logo --collected_files --structure_from_filename "dirpath"

Handled special case FeCrSn with unique directory structure (later adjustment to half-metal) by copying necessary files following the above `--colected_files` format.

# Usage: compression

To compress all plots in the HeuslerPlot subdirectory, run:

    python3 compress.py

To compress all plots in a different directory with path `dir_path`, run:

    python3 compress.py "dir_path"

PNG files are already compressed losslessly with DEFLATE. The additional compression employed here is lossy and removes colors. A high compression ratio (on the order of 6 to 1) with minimal impact on the quality of the plots generated by plot.py can be obtained by shrinking the number of colors to 4 (which is used by default by compress.py).
