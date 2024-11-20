<h1 align="center">
PhaseTracer
</h1>

<div align="center">
<i>Trace cosmological phases, phase transitions, and gravitational waves in scalar-field theories</i>
</div>
<br>
<div align="center">
<img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/PhaseTracer/PhaseTracer/cmake-single-platform.yml">
<img alt="GitHub License" src="https://img.shields.io/github/license/PhaseTracer/PhaseTracer">
<img alt="Static Badge" src="https://img.shields.io/badge/arXiv-2003.02859-blue?link=https%3A%2F%2Farxiv.org%2Fabs%2F2003.02859">
</div>
<br>

**PhaseTracer** is a C++14 software package for tracing cosmological phases, finding potential phase transitions, computing the bounce action, and plotting the gravitational wave spectrum for Standard Model extensions with any number of scalar fields.

## Dependencies

You need a C++14 compliant compiler and our dependencies. The dependencies can be installed by

*Ubuntu/Debian*

    sudo apt install libalglib-dev llibnlopt-cxx-dev libeigen3-dev libboost-filesystem-dev libboost-log-dev libgsl-dev
    
*Fedora*

    sudo dnf install alglib-devel nlopt-devel eigen3-devel boost-devel gsl-devel
    
*Mac*

    brew install alglib nlopt eigen boost gsl

If alglib is not found, see https://github.com/S-Dafarra/alglib-cmake

## Building

To build the shared library and the examples:

    git clone https://github.com/PhaseTracer/PhaseTracer
    cd PhaseTracer
    mkdir build
    cd build
    cmake ..
    make

## Running

If the build was succesful, run the examples and tests with:

    cd ..
    ./bin/run_1D_test_model
    ./bin/run_2D_test_model
    ./bin/scan_Z2_scalar_singlet_model
    ./bin/unit_tests

## FlexibleSUSY
<details>
<summary>Click me</summary>

To build the example `THDMIISNMSSMBCsimple` with FlexibleSUSY:

    cmake -D BUILD_WITH_FS=ON ..
    make

Then run the example with:

    cd ..
    ./bin/run_THDMIISNMSSMBCsimple

FlexibleSUSY has additional dependencies and will report errors if
these are not present. See the FlexibleSUSY documentation for details
and/or follow the suggestions from the cmake output.
</details>

## BSMPT
<details>
<summary>Click me</summary>
To build the examples with BSMPT:

    cmake -D BUILD_WITH_BSMPT=ON ..
    make

Then run the examples with:

    cd ..
    ./bin/run_R2HDM
    ./bin/run_C2HDM
    ./bin/run_N2HDM

Please note that the BSMPT examples in PhaseTacer are just for checking that PhaseTacer and BSMPT can give consistent results.  Unsuccessful compilation of BSMPT will not affect other examples and BSMPT is not neccessary for PhaseTracer users unless they wish to use potentials from BSMPT.
</details>
    
## Citing

If you use PhaseTracer, please cite the accompanying manual

    @article{Athron:2020sbe,
        author = "Athron, Peter and Bal\'azs, Csaba and Fowlie, Andrew and Zhang, Yang",
        title = "{PhaseTracer: tracing cosmological phases and calculating transition properties}",
        eprint = "2003.02859",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        reportNumber = "CoEPP-MN-20-3",
        doi = "10.1140/epjc/s10052-020-8035-2",
        journal = "Eur. Phys. J. C",
        volume = "80",
        number = "6",
        pages = "567",
        year = "2020"
    }


