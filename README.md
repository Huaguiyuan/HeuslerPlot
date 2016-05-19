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

    sudo python3 setup.py install

To have changes to the source reflected immediately:

    sudo python3 setup.py develop

# Usage

To plot bands in a directory at "dirpath" with a structure dirpath/e29/Co2MnSi/BANDS (where there may be many e## and e##/systemname subdirectories) use (where `structure_type` is set appropriately):

    python3 plot.py --searchEs --logo --structure_type "L21" "dirpath"

If the bands subdirectory has a different name (such as dirpath/e29/Co2MnSi/bands instead of .../BANDS) but the directory structure is otherwise the same as above (this is the structure for half-Heusler cubic compounds); also exclude compounds containing `bands_t` subdirectory indicating the compound is tetragonal:

    python3 plot.py --searchEs --logo --structure_type "C1b" --bands_name "bands" --exclude_if_present "bands_t" "dirpath"

If there is an additional subdirectory containing tetragonal scf and bands runs, e.g. `dirpath/e29/Co2MnSi/bands_t` gives the tetragonal SCF and `dirpath/e29/Co2MnSi/bands_t/bands` gives the tetragonal bands (this is the structure for half-Heusler tetragonal compounds):

    python3 plot.py --searchEs --logo --structure_type "tetragonal" --bands_subdir_path "bands_t" --bands_name "bands" "dirpath"

To plot bands where files have been collected in a directory with filenames of the form `EIGENVAL_Xa_Fe2CoAs` (from bands calculation), `OUTCAR_Xa_Fe2CoAs` (from SCF), `OSZICAR_Xa_Fe2CoAs` (from SCF) (inverse Heuslers have been collected like this):

    python3 plot.py --logo --collected_files --structure_from_filename "dirpath"
