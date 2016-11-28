Porosity-tools
For more details see: Kulasinski K. and Guyer R, Quantification of Nanopore Networks: Application to Amorphous Polymers, J Phys Chem C, 2016

Python script to quantify the pore network of a Molecular Dynamics system.
Copyright Karol Kulasinski 2016

Required python libraries:
* Python 2.x
* Numpy
* Scipy
* Matplotlib
* Getopt

Installation:
Copy the python files to your local directory.

Usage:

python potosity-tools.py [OPTIONS]

OPTIONS:
     -f <input.gro>  # input structure file, in .gro format
     
        <grid.npy>   # or previously calculated grid in .npy format
        
     -o <out>        # generic name for output files
     
     -d <output>     # output directory
     
     -r <float>      # resolution of the grid elements in nm, default: 0.01 nm
     
     -m <float>      # radius of the probe molecule in nm, default is 0.14 nm (H2O radius). If radius is 0 Van der Waals surface is probed.
                       
     -g <gro|npy|no> # save the grid to .gro file, numpy .npy file, or don't (default).
                       
     -p <y/n>        # if your structure is fully periodic, default is yes.
     
     -c <y/n>        # calculate chord and diameter distribution
     
     -t <x|y|z>      # calculate tortuosity in given direction
     
Example usage:

python porosity-tools.py -f conf.gro -r 0.05 -d foodir -o cellulose -g npy

python porosity-tools.py -f foodir/cellulose_grid.npy -d foodir -o cellulose -r 0.05 -c y

python porosity-tools.py -f foodir/cellulose_grid.npy -d foodir -o cellulose -r 0.05 -t x
