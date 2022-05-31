# Microfluidics-Abstraction-Levels
Source code for the paper _Simulation of Pressure-Driven and Channel-Based Microfluidics on Different Abstract Levels: A Case Study_

The source code contains three main elements:
* Source code for the 1D approach
* Source code for the Finite Volume Method (_FVM_), using [OpenFOAM](https://openfoam.org/download/8-ubuntu/).
* Source code for the Lattice-Boltzmann Method (_LBM_), using [Palabos](https://palabos.unige.ch/).

### 1D
For the 1D approach, there are two python scripts in this repository, namely for the Non-Newtonian and Fluid-Mixing use cases. The scripts can be run using python3.

### FVM
To run the 2D and 3D use cases, using the Finite Volume Method, please be sure to have a working installation of OpenFOAM v8 on your system. The cases can be run using standard OpenFOAM solvers.

### LBM
To run the 2D and 3D use cases, using the Lattice-Boltzmann Method, please be sure to have a working installation of Palabos on your system. Additionally, in the CMake files, be sure to correctly refer to the location on your system where Palabos is installed. The cases must first be compiled in the `build/` folder. In Linux the makefile can be created by using the command `cmake ..` in the `build/` folder. The code can be compiled using the `make` command, after which an executable will be created, which will run the cases using the parameters from `param.xml`. For more information on how to run code using Palabos (also on Windows), please have a look at the Palabos [User Guide](https://palabos.unige.ch/files/9515/6509/3036/Palabos_UserGuide.pdf).
