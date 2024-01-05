## Code for J.M. Osborne  "An adaptive numerical method for multi-cellular simulations of organ development and disease" 

This project contains the code necesary for running the simulations presented in "An adaptive numerical method for multi-cellular simulations of organ development and disease"

...

Before looking at this, you may wish to look at some of the basic user tutorials for Chaste https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.

## Getting the code and installing dependencies 

Before running these examples you will need to install Chaste's dependencies (https://chaste.github.io/getting-started/) and the source code for the latest release (release 2021.1) https://github.com/Chaste/Chaste.
The easiest way to do this is using an Ubuntu machine (or an Ubuntu virtual machine) as discussed on https://chaste.github.io/docs/installguides/ubuntu-package/. 
Note that Chaste is only fully supported on Linux/Unix systems, so users of Windows or MacOS we reccomend using Docker (https://github.com/Chaste/chaste-docker).

Once the chaste dependencies and installed and the source code is downloaded go to the `projects` folder and use the command 
`git clone https://github.com/jmosborne/CellBasedNumericalMethods.git`  
to download the code for this project.

Now the project should be installed, and everything should compile and run correctly. 
You can now run the tests or simulations, or create your own test suites.

## Documentation
There are two folders - `src` and `test`.
 1. The `src` folder contains the classes necessary to run the simulation. These define the additional numerical methods not in the core chaste code.
  * `AdamsMoultonNumericalMethod.hpp (cpp)` - AM2 implicit method
  * `BackwardEulerNumericalMethod.hpp (cpp)` - BE implicit method
  * `MidpointNumericalMethod.hpp (cpp)` - Midpoint (RK2) explicit method
  * `RK3NumericalMethod.hpp (cpp)` - RK3 method
  * `RK4NumericalMethod.hpp (cpp)` - RK3 method
  * `ForTests - extra code used to build the exemplar simulations, boundary conditions, cell cycle models etc.
 
 2. The `test` folder contains:
  * `TestNumerics.hpp` - this file can be run to generate the simulation results used to make all 1D Compression 2D Proliferation 2D Monloayer and 3D Spheroid results  in the paper.
  * `TestNumericsCrypt3D.hpp` - this file can be run to generate the simulation results used to make all the 3D Organ results  in the paper.
  * `run_1d_sweeps.sh` - Script to run multiple 1D Compression or. 1D Proliferation simulations.
  * `run_2d_3d_sweeps.sh` - Script to run multiple 2D Monolayer or 3D Spheroid simulations.
  * `run_3d_crypt_sweeps.sh` - Script to run multiple 3D Organ simulations.
  
## Running tests
You can then run tests and simulations with (note this assumes the file structure used in the Chaste Docker),
```
cd /lib
cmake ./src
make TestNumerics

```

Note that this will only compile the test. The following commands will run some of the parameters sweeps detailed in the paper:

```
cd projects/CellBasedNumericalMethods/test/
sh run_script.sh
```

**NB**: the paper was developed with the release 2021.1, it will not work with with release 2019.1 or under.

For further information on using Chaste, see the extensive guide material (https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides).
You may also wish to look at some of the basic user tutorials (https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials/). Note these links will be updated to github website soon.
