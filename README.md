# AnasenSim

AnasenSim is a Monte Carlo nuclear reaction kinematics simulation with the ANASEN detector geometry, designed around the active target configuration. AnasenSim uses the CAtima library to calculate the energy loss through the gas volume and detector elements. The geometry of the ANASEN detector componets is adjustable to different experimental configurations.

## Building

To download AnasenSim use `git clone --recursive https://github.com/gwm17/AnasenSim.git`. To build AnasenSim on Mac or Linux, enter the AnasenSim repository and run the following commands:

- `mkdir build`
- `cd build`
- `cmake .. && make`

## Simulation configurations

To specify the reaction of interest to AnasenSim, a lightweight text input file is used. An example of the format is given with the repository (input.txt). In general the input requires the specification of the target gas, the reaction chain, and a location to which data will be written. For the reaction specification, AnasenSim by default can calculate Reactions of up to 3 steps (one primary reaction and subsequent decays). Other configurations will require modification of the kinematics simulation.

To run the simulation use the following command structure: `./bin/AnasenSim <your_input_file>`

## Plotting

AnasenSim comes with a pre-packaged generic plotter (Plotter). This tool will take a simulation file and generate kinematics plots for the nuclei. It is very generic, so typically one would want to either tweak it to fit a specific use case, or design a custom plotter from scratch. Note that AnasenSim data is written using a ROOT dictionary, so a new plotter will need to link against the dictionary (found in lib).

## Requirements

- ROOT version which is compatible with CMake (tested on 6.25.06)
- CMake version  > 3.12
