# Installation

Instructions for installation on Debian derivatives (Debian/Ubuntu/Mint/etc.).

Requires numpy and matplotlib:

    sudo apt-get install python3-numpy python3-matplotlib python3-tk

Get setuptools and install using setup.py:

    sudo apt-get install python3-setuptools
    sudo python3 setup.py install

To have changes to the source reflected immediately:

    sudo python3 setup.py develop

# Usage

To plot bands in a directory at "dirpath" with a structure dirpath/e29/Co2MnSi/BANDS (where there may be many e## and e##/systemname subdirectories) use:

    python3 plot.py --searchEs --logo "dirpath"
