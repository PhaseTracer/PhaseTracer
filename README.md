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

**PhaseTracer** is a C++17 software package for tracing cosmological phases, finding potential phase transitions, computing the bounce action, and plotting the gravitational wave spectrum for Standard Model extensions with any number of scalar fields. The full pipeline can be run from C++ or from Python.

## Quick start

    git clone https://github.com/PhaseTracer/PhaseTracer
    cd PhaseTracer
    ./scripts/install_dependencies.sh            # or --local (no sudo); add --python for Python

    # C++
    cmake -B build && cmake --build build -j 4
    ./bin/run_1D_test_model

    # Python
    python3 -m venv ~/.venvs/pt && source ~/.venvs/pt/bin/activate
    pip install .

## Dependencies

You need a C++17 compliant compiler and our dependencies (Boost, ALGLIB, Eigen3, NLopt and GSL). The easiest way to get them is the install script, which has two modes:

    ./scripts/install_dependencies.sh            # system packages via apt, dnf or brew (needs sudo)
    ./scripts/install_dependencies.sh --local    # build pinned versions into ./.deps (no sudo)

The local mode builds Boost, NLopt and GSL as static libraries and ALGLIB as a shared one. CMake
picks up `./.deps` automatically; keep the directory, as the build loads `libalglib.so` from there.
- `--download-only`: fetch the sources for a machine without internet access.
- `--python`: also install what the Python interface needs.
- `--help`: list all options.

Alternatively, the dependencies can be installed by hand:

*Ubuntu/Debian*

    sudo apt install libalglib-dev libnlopt-cxx-dev libeigen3-dev libboost-filesystem-dev libboost-log-dev libgsl-dev
    
*Fedora*

    sudo dnf install alglib-devel nlopt-devel eigen3-devel boost-devel gsl-devel
    
*Mac*

    brew install nlopt eigen boost gsl alglib

If alglib is not found, see https://github.com/S-Dafarra/alglib-cmake

## Installation (C++)

To build the shared libraries and the examples:

    git clone https://github.com/PhaseTracer/PhaseTracer
    cd PhaseTracer
    ./scripts/install_dependencies.sh
    cmake -B build
    cmake --build build -j 4

This is equivalent to `mkdir build && cd build && cmake .. && make -j 4`. The libraries are written
to `lib/` and the programs to `bin/`. Raise `-j` if you have more cores and memory (each job can
need about 1 GB).

**Useful CMake options** (`cmake -B build -D<option>=<value>`):
- `-DCMAKE_BUILD_TYPE=Release`: optimised build, recommended for timing runs and scans.
- `-DBUILD_WITH_HG=OFF`: build without HydroGrav, the SoundShell GW backend. It is on by default and is cloned at configure time.
- `-DBUILD_PYTHON=ON`: also build the Python module, for development (see [python/README.md](python/README.md)).
- `-DPT_DEPS_PREFIX=<dir>`: dependencies installed with `./scripts/install_dependencies.sh --local --prefix <dir>`.

### Running

If the build was succesful, run the examples and tests with:

    ./bin/run_1D_test_model
    ./bin/run_2D_test_model
    ./bin/scan_Z2_scalar_singlet_model
    ./bin/run_thermo_finder
    ./bin/run_hydrograv
    ./bin/run_runner
    ./bin/unit_tests
    
If you want to see debugging information or obtain plots of the phases and potential for the first two examples above you can add the -d flag, i.e.

    ./bin/run_1D_test_model -d 
    ./bin/run_2D_test_model -d

## Installation (Python)

***Note the Python installation will ultimately be handled using PyPI once a stable release build is ready. These instructions apply in the meantime.***

The `phasetracer` Python package runs the same pipeline from Python, with the models shipped with
PhaseTracer or with potentials you write in Python. After installing the dependencies (as above,
with `--python`), install it into a virtual environment:

    python3 -m venv ~/.venvs/pt
    source ~/.venvs/pt/bin/activate
    pip install -U pip
    pip install .
    python -c "import phasetracer"

pip runs CMake itself and builds only the libraries and the module, which takes a few minutes.
After that, `import phasetracer` works from any directory or notebook that uses this environment.

- **SoundShell GW backend:** `pip install . -C cmake.define.BUILD_WITH_HG=ON`.
- **Update:** `git pull && pip install .`.
- **Uninstall:** `pip uninstall phasetracer`.
- **Jupyter:** run `pip install ipykernel && python -m ipykernel install --user --name pt`.

### Running

A minimal run:

    import phasetracer as pt
    from phasetracer.models import TwoDimModel

    config = pt.Config()
    config.pipeline.to_print = False

    runner = pt.Runner(TwoDimModel(), config)

    status = runner.run()
    print(status)
    
    for tps in runner.get_thermal_parameters():
        print(tps.TC, tps.percolation.temperature, tps.percolation.alpha)
    print(runner.get_spectra()[0].SNR)

See [python/README.md](python/README.md) for the full interface, including how to define a model
in Python, and `python/examples/` for complete scripts.

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
        month = "12",
        year = "2024"
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


