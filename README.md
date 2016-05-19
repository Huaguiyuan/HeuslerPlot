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

To plot bands in a directory at "dirpath" with a structure dirpath/e29/Co2MnSi/BANDS (where there may be many e## and e##/systemname subdirectories) use:

    python3 plot.py --searchEs --logo "dirpath"
