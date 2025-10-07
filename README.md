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
<a href="https://arxiv.org/abs/2003.02859"><img alt="Static Badge" src="https://img.shields.io/badge/arXiv-2003.02859-blue"></a>
<a href="https://arxiv.org/abs/2412.04881"><img alt="Static Badge" src="https://img.shields.io/badge/arXiv-2412.04881-blue"></a>
</div>
<br>

**PhaseTracer** is a C++14 software package for tracing cosmological phases, finding potential phase transitions, computing the bounce action, and plotting the gravitational wave spectrum for Standard Model extensions with any number of scalar fields.

## Dependencies

You need a C++14 compliant compiler and our dependencies. The dependencies can be installed by

*Ubuntu/Debian*

    sudo apt install libalglib-dev libnlopt-cxx-dev libeigen3-dev libboost-filesystem-dev libboost-log-dev libgsl-dev
    
*Fedora*

    sudo dnf install alglib-devel nlopt-devel eigen3-devel boost-devel gsl-devel
    
*Mac*

    brew install nlopt eigen boost gsl alglib

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
    
If you want to see debugging information or obtain plots of the phases and potential for the first two examples above you can add the -d flag, i.e.

    ./bin/run_1D_test_model -d 
    ./bin/run_2D_test_model -d


## BubbleProfiler
<details>
<summary>Click me</summary>

To use `BubbleProfiler` for calculation of bounce action:

    cmake -D BUILD_WITH_BP=ON ..
    make

Then run the example with:

    cd ..
    ./bin/run_BP_2d
    ./bin/run_BP_scale 1 0.6 200

or in other examples by setting
    
    PhaseTracer::ActionCalculator ac(model);
    ac.set_action_calculator(PhaseTracer::ActionMethod::BubbleProfiler);

</details>


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
 
    @article{Athron:2024xrh,
        author = "Athron, Peter and Balazs, Csaba and Fowlie, Andrew and Morris, Lachlan and Searle, William and Xiao, Yang and Zhang, Yang",
        title = "{PhaseTracer2: from the effective potential to gravitational waves}",
        eprint = "2412.04881",
        archivePrefix = "arXiv",
        primaryClass = "astro-ph.CO",
        doi = "10.1140/epjc/s10052-025-14258-y",
        journal = "Eur. Phys. J. C",
        volume = "85",
        number = "5",
        pages = "559",
        year = "2025"
    }

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
