[![Build Status](https://travis-ci.org/PhaseTracer/PhaseTracer.svg?branch=master)](https://travis-ci.org/PhaseTracer/PhaseTracer)

# PhaseTracer

`PhaseTracer` is a C++ software package for mapping out cosmological phases, and potential transitions between them, for Standard Model extensions with any number of scalar fields.

## Citing

If you use `PhaseTracer`, please cite the accompanying manual

    @article{PhaseTracer,
      author         = "Athron, Peter and Balazs, Csaba and Fowlie, Andrew and Zhang, Yang",
      title          = "{PhaseTracer: tracing cosmological phases and calculating transition properties}",
      year           = "2020",
      eprint         = "2003.02859",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "CoEPP-MN-20-3",
    }
    

## Quickstart

### Requirements

The following libraries are required to build `PhaseTracer`:

* [nlopt](http://ab-initio.mit.edu/wiki/index.php/NLopt/)  is used for locating the extrema of potentials.
* [Eigen](https://eigen.tuxfamily.org) is an excellent linear algebra library
* The following [BOOST](http://www.boost.org/) libraries are required:

  * `Boost.Filesystem`
  * `Boost.Log`

* Our `EffectivePotential` library is built automatically in our build process, but itself requires

  * [alglib](http://www.alglib.net/) for spline interpolation

On Ubuntu, the dependencies can be installed by

    sudo apt install libalglib-dev libnlopt-dev libeigen3-dev libboost-filesystem-dev libboost-log-dev

On Ubuntu 20.04.1 LTS, nlopt for C++ needs to be installed by
    sudo apt install libnlopt-cxx-dev
    
On Mac, they can be installed by
    
    brew install alglib nlopt eigen boost

In case [alglib](http://www.alglib.net/) is not found in package manager, try [alglib-cmake](https://github.com/S-Dafarra/alglib-cmake).

### Building

To build the shared library and the examples:

    mkdir build
    cd build
    cmake ..
    make

### Running

If the build was succesful, run the examples with:

    cd ..
    ./bin/run_1D_test_model
    ./bin/run_2D_test_model
    ./bin/scan_Z2_scalar_singlet_model

### Building and running examples with `FlexibleSUSY` and `BSMPT`

To build the example `THDMIISNMSSMBCsimple` with `FlexibleSUSY`:

    cmake -D BUILD_WITH_FS=ON ..
    make

Then run the example with:

    cd ..
    ./bin/run_THDMIISNMSSMBCsimple

FlexibleSUSY has additional dependencies and will report errors if
these are not present. See the FlexibleSUSY documentation for details
and/or follow the suggestions from the cmake output.

To build the examples with BSMPT:

    cmake -D BUILD_WITH_BSMPT=ON ..
    make

Then run the examples with:

    cd ..
    ./bin/run_R2HDM
    ./bin/run_C2HDM
    ./bin/run_N2HDM

Please note that the BSMPT examples in PhaseTacer are just for checking that PhaseTacer and BSMPT can give consistent results.  Unsuccessful compilation of BSMPT will not affect other examples and BSMPT is not neccessary for PhaseTracer users unless they wish to use potentials from BSMPT.

