## RheoCell: Multi-Phase Field model for cell monolayers

Simulation code used for (preprint arXiv:2508.06461) that models cell monolayers as active droplets using the multi-phase field model.

## Building

We use cmake to build the executable files. The input is the following from the main folder:

```
mkdir build
cd build
cmake ..
make
```

We use the Eigen3 package for matrix operations (https://libeigen.gitlab.io/). This package can be used as header only, but should be supplied to cmake.

## Running

The code is run from the command line and an input file must always be given as the first argument:

`./rheocell input`

An executable that generates different initial conditions specified by the input file is also included:

`./confgenerator input NNxNN`

where NN should be replaced by integers and will specify the the length and width of the simulation box, respectively.

## Examples

Some example input files are included in the `examples` directory. The `scripts` directory includes multiple python scripts that are used to analyse the results. In particular, the `read_conf.py` allows the user to visualize a configuration by supplying the topology file (see `examples` directory) and the configuration file extracted from a run (e.g. `last_conf.dat`).
