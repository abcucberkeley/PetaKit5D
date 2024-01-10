# LLSM5DTools

Tools for efficient and scalable processing of petabyte-scale 5D live images or large specimen images from lattice light-sheet microscopy (LLSM) and other light sheet microscopies (also other imaging modalities). It is featured by fast image readers and writers (for Tiff and Zarr), combined image deskew/rotation, instantly converged Richardson-Lucy (RL) deconvolution and scalable Zarr-based stitching. It also contains some other useful tools, including 3D point detection and point tracking (based on Aguet, Upadhyayula et al., 2016), cropping, resampling, max projection, PSF analysis and visualization and more.

## Usage

The tools have been tested with MATLAB R2022b-R2023a for Linux (Ubuntu), Windows, and MacOS. Here are the steps to use the software:
1. Get the source code by either cloning the GitHub repository or downloading the ZIP file. If downloading the zip file, unzip the file to a directory.
2. Launch MATLAB, and navigate to the software's root directory and add the software to path with `setup.m` in the command window.
```
    setup
```
3. Create a script or function to setup the workflows for your image processing tasks by calling related functions. You may follow the examples in the [demos](https://github.com/abcucberkeley/LLSM5DTools/tree/master/demos). The documentation of the parameters can refer to the [GUI wiki page](https://github.com/abcucberkeley/LLSM_Processing_GUI/wiki) or the parameter list in the related functions.


## Demos
The main demos for the paper:
- `demo_generic_computing_framework.m`: demo to illustrate how to use generic computing framework for user-defined functions.
- `demo_fast_tiff_zarr_readers_writers.m`: demo to illustrate Cpp-Tiff and Cpp-zarr readers and writers, and comparison with conventional readers and writers. 
- `demo_geometric_transformation.m`: demo to illustrate how to run deskew/rotation, and comparision between separated and combined deskew/rotation. 
- `demo_RL_deconvolution.m`: demo to illustrate how to run deconvolution, and comparision between conventional RL, WB (Guo et al. 2020), and OMW methods. 
- `demo_zarr_stitching.m`: demo to illustrate how to run stitching in both skewed space and DSR space. 
- `demo_large_scale_processing.m`: demo to illustrate how to setup large-scale processing for stitching, deconvlution, and deskew/rotation. 


## GUI
The software has an easy-to-use Graphical User Interface (GUI) without writing any code, which can be find here: [LLSM_Processing_GUI](https://github.com/abcucberkeley/LLSM_Processing_GUI). The GUI has support for Windows, MacOS and Linux (Ubuntu). For instructions on installation and use of the LLSM_Processing_GUI, vist the [GUI wiki page](https://github.com/abcucberkeley/LLSM_Processing_GUI/wiki). 

## Fast Tiff/Zarr readers and writers
We also created independent repositories for Tiff/Zarr readers and writers for users who only need those functions: [Cpp-Tiff](https://github.com/abcucberkeley/cpp-tiff) and [Cpp-zarr](https://github.com/abcucberkeley/cpp-zarr). 

Based on these readers and writers, we also developed a Fiji Plugin for faster reading and writing of Tiff and Zarr files within Fiji, which can be accessed from [Parallel_Fiji_Visualizer](https://github.com/abcucberkeley/Parallel_Fiji_Visualizer).


## Reference:
Please cite our paper if you found the software is useful for your research: 

`Ruan X., Mueller, M., Liu G., Görlitz F., Fu T., Milkie D., Lillvis J., Killilea A., Betzig, E. and Upadhyayula S. (2024) Image processing tools for petabyte-scale light sheet microscopy data. bioRxiv. bioRxiv 2023.12.31.573734; doi: https://doi.org/10.1101/2023.12.31.573734`
