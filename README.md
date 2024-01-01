# LLSM5DTools

Tools for the processing of petabyte-scale 5D live images or images from large specimens from lattice light-sheet microscopy (LLSM) and other light sheet microscopies (also other imaging modalities). It is featured by fast image readers and writers (for Tiff and Zarr), combined image deskew/rotation, instantly converged Richardson-Lucy (RL) deconvolution and scalable Zarr-based stitching. It also contains some other useful tools, including 3D point detection and point tracking (based on Aguet, Upadhyayula et al., 2016), cropping, resampling, max projection, PSF analysis and visualization and so on.

## Usage

The tools have been tested with MATLAB R2022b-R2023a. Here are the steps to use the software:
1. Download the source code of the software with git repository or the zip file. If downloading the zip file, unzip the file to a directory; 
2. Launch MATLAB, and go to the root directory of the software and add the software to path with `setup.m` in the command window;
```
    setup
```
3. Create a script or function to setup the workflows for your image processing tasks by calling related functions. You may follow the examples in the [demos](https://github.com/abcucberkeley/LLSM5DTools/tree/dev/demos). The documentation of the parameters can refer to the [GUI wiki page](https://github.com/abcucberkeley/LLSM_Processing_GUI/wiki) or the parameter list in the related functions.

## GUI
The software also has an easy-to-use Graphical User Interface (GUI) without writing any code, which can be find here: [LLSM_Processing_GUI](https://github.com/abcucberkeley/LLSM_Processing_GUI). The GUI has support for Windows, MacOS and Linux (Ubuntu). For instructions on installation and use of the LLSM_Processing_GUI, vist the [GUI wiki page](https://github.com/abcucberkeley/LLSM_Processing_GUI/wiki). 

## Reference:
Please cite our paper if you found the software is useful for your research: 
```
Ruan X., Mueller, M., Liu G., Görlitz F., Fu T., Milkie D., Lillvis J., Killilea A., Betzig, E. and Upadhyayula S. (2023) Image processing tools for petabyte-scale light sheet microscopy data. bioRxiv. 
```
