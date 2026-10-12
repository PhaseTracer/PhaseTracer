# PhaseTracer Python interface

The `phasetracer` Python package runs the full PhaseTracer pipeline (phases, transitions,
thermal parameters and the gravitational wave spectrum) from Python. You can use the models
shipped with PhaseTracer, or define your own potential in Python.

```python
import phasetracer as pt
from phasetracer.models import TwoDimModel

config = pt.Config()
config.phase_finder.seed = 1
config.pipeline.to_print = False

runner = pt.Runner(TwoDimModel(), config)
status = runner.run()
print(status)                                   # [Success]

for tps in runner.get_thermal_parameters():
    p = tps.percolation
    print(tps.TC, p.temperature, p.alpha, p.betaH)

spectrum = runner.get_spectra()[0]
spectrum.frequency, spectrum.total_amplitude, spectrum.SNR    # SNR = [LISA, Taiji]
```

## Installation

Install the package from a clone of the repository:

```bash
git clone https://github.com/PhaseTracer/PhaseTracer.git
cd PhaseTracer
./scripts/install_dependencies.sh --python     # system packages (sudo), or:
./scripts/install_dependencies.sh --local      # build the C++ dependencies into ./.deps (no sudo)

python3 -m venv ~/.venvs/pt
source ~/.venvs/pt/bin/activate
pip install -U pip
pip install .
```

`pip install .` runs CMake for you, in a temporary build environment. It builds only the PhaseTracer
libraries and the Python module (no examples or unit tests), and installs the package into the
active environment. Due to needing to compile the entire library, the first build may take a few minutes.

The installation can be checked using:
```bash
cd ~ && python -c "import phasetracer; print(phasetracer.__file__)"
```

From then on, `import phasetracer` works from any script, directory, or notebook that uses this
environment. The repository is no longer needed at run time. The package contains the
PhaseTracer libraries:
- **With `--local`:** it also contains ALGLIB, and only needs libstdc++ and libgomp from the system.
- **With system packages:** it also needs the system NLopt, ALGLIB, GSL and Boost libraries.

Either way the install is for the machine it was built on (or an identical one).

**Options and maintenance:**
- **SoundShell GW backend (HydroGrav):** `pip install . -C cmake.define.BUILD_WITH_HG=ON`. This needs network access to clone HydroGrav, or run `./scripts/install_dependencies.sh --hydrograv` first.
- **Dependencies in a non-default prefix:** `pip install . -C cmake.define.PT_DEPS_PREFIX=/path/to/prefix`.
- **Update:** `git pull && pip install .`.
- **Uninstall:** `pip uninstall phasetracer`.
- **Jupyter:** run `pip install ipykernel && python -m ipykernel install --user --name pt`, then choose the `pt` kernel.

### Working on the code

To change PhaseTracer itself, whether its C++ or the bindings, see [Development](#development)
below. It covers a development build, rebuilding after changes, and adding C++ models to
`phasetracer.models`.

## Usage

### Config

`pt.Config()` holds every setting of the pipeline. It has one group per class, and the defaults
are those of the C++ classes:

| Group | Class | Examples |
|---|---|---|
| `config.phase_finder` | PhaseFinder | `seed`, `t_low`, `t_high`, `lower_bounds`, `check_hessian_singular` |
| `config.transition_finder` | TransitionFinder | `TC_tol_rel`, `assume_only_one_transition` |
| `config.action` | ActionCalculator | `method`, `PD_xtol`, `PD_deformation_npoints` |
| `config.thermo_finder` | ThermoFinder | `percolation_target`, `action_spline_evaluations`, `transition_filter`, `prefactor_function` |
| `config.gravwave` | GravWaveCalculator | `gw_method`, `min_frequency`, `max_frequency`, `vw` |
| `config.pipeline` | the run itself | `stop_after`, `throw_on_error`, `log_level`, `to_print`, `scale` |

- **Editing and checking:**
  - Groups are edited in place: `config.phase_finder.seed = 0`.
  - `print(config)` lists every value.
  - `config.validate()` checks for inconsistent settings before anything runs.
  - `copy.deepcopy(config)` makes an independent copy.
- **Units:** `config.pipeline.scale` (`pt.Scale.MEV`, `pt.Scale.GEV` (default) or `pt.Scale.TEV`) is the unit of the
  model's field values and temperatures. It sets the Planck mass used for the Hubble rate and the conversion
  applied to the GW spectrum, so frequencies always come out in Hz.
- **List settings** (`lower_bounds`, `upper_bounds`, `guess_points`): assign a whole new list. Calling `.append()` on them does not change the config.
- **Callbacks:** `transition_filter` and `prefactor_function` take Python callables:
  ```python
  config.thermo_finder.transition_filter = lambda transitions: [t for t in transitions if t.TC < 150]
  config.thermo_finder.prefactor_function = lambda T, S_over_T, action_result: T**4
  ```

### Running and errors

- **`runner.run()` returns a `RunStatus`**, which is truthy on success. Otherwise it gives the stage where the run stopped and the reason:
  - `status.code`, `status.stage`, `status.message`;
  - `status.warnings`, for example transitions whose thermal parameters failed.
- **The status codes:**
  - The `No*` codes (`NoTransitions`, `NoThermalParameters`, `NoSpectra`, …) are physics outcomes.
  - The `*Failed` codes mean a stage threw.
  - `InvalidConfig` means a setting was rejected.
- **Raising instead:** with `config.pipeline.throw_on_error = True`, `run()` raises `pt.RunnerError`, and the status is available as `error.status`.
- **Stopping early:** `config.pipeline.stop_after = pt.Stage.TransitionFinder`, for example, skips the expensive thermal and GW stages.
- **Logging:** `pt.set_log_level("debug")` shows PhaseTracer's log output. The default is `fatal`.

### Results and stage objects

| Call | Returns |
|---|---|
| `runner.get_phases()`, `get_transitions()`, `get_spectra()` | lists of `Phase`, `Transition`, `GravWaveSpectrum` (copies) |
| `runner.get_thermal_parameters()` | one `ThermalParameterSet` per analysed transition: `TC`, the `onset` / `nucleation` / `percolation` / `completion` milestones, `profiles`, `decay_rate()` (S(T), Γ(T), …) and `equation_of_state()` |
| `runner.phase_finder()`, `action_calculator()`, `transition_finder()`, `thermo_finder()`, `gravwave_calculator()` | the stage objects, e.g. `runner.action_calculator().get_action(phase1, phase2, T)` |

- **Plain data stays valid:** phases, transitions, milestones and spectra are copies, so they remain usable after re-runs.
- **Stage objects and thermal parameter sets belong to the Runner.** After `runner.run()` is called again, using one obtained earlier raises `RuntimeError("stale handle …")` instead of returning outdated data. They keep the Runner alive.
- **Re-running:** `runner.config` can be edited between runs, and every `run()` starts from scratch.

### Models

`phasetracer.models` contains the models shipped with PhaseTracer: compiled C++ models, and a few
written in Python. `help(phasetracer.models.<Name>)` shows each one's constructor arguments.

| Family | Model | Notes |
|---|---|---|
| Test models | `OneDimModel()` | one field, analytic TC (`get_TC_from_expression`) |
| | `TwoDimModel()` | one-loop, two fields; reaches percolation and gives a GW spectrum |
| | `Z2ScalarSingletModel()` | high-temperature xSM with analytic TC; `set_m_s`, `set_lambda_hs` |
| xSM | `xSM_MSbar.from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, ...)` | MS-bar xSM with tadpole conditions applied; also `xSM_MSbar(lambda_hs, lambda_s, ms)` |
| | `xSM_OSlike(lambda_hs, lambda_s, ms)` | on-shell-like renormalisation |
| | `xSM_HT(lambda_hs, lambda_s, m_s)` | high-temperature expansion, analytic TC |
| RSS | `RSS(muh_sq, mus_sq, lh, ls, c1, c2, b3, muh_sq_tree_EWSB, mus_sq_tree_EWSB, vs, ms, theta)` | singlet without Z2 symmetry; `RSS()` + `setParameters(...)` also works |
| | `RSS_HT(muh_sq, mus_sq, lh, ls, c1, c2, b3, vs, ms, theta)` | its high-temperature expansion |
| Other BSM | `ConformalDarkU1(g, y, v)` | conformal dark U(1), arXiv:2502.19478 |
| | `SuperCoolModel(kappa, mu_sq, lambda_, mu0_sq)` | strong supercooling; or `SuperCoolModel(input_file)` |
| | `ScalarYukawa3D()` | 3d scalar-Yukawa toy model; set `s`, `msq`, `g`, `lam`, `mf`, `y` |
| 2HDM | `THDM(tan_beta, m12_sq, lambda1, ..., lambda5, yukawa_type=2, RG_scale=-1)` | the TransitionListener2 benchmark 2HDM with counterterms |
| Python | `QuarticThermal(D, E, lam, T0)` | polynomial thermal potential, analytic `TC()` |
| | `DarkHiggs(lam, v, g, y)` | one-loop dark Higgs with a dark gauge boson and a Yukawa fermion |
| | `TL2HDM(...)` | TransitionListener's own 2HDM. Needs a TransitionListener checkout (set `TRANSITIONLISTENER_ROOT`), and is only imported when you first use it |

`OneDimModel`, `Z2ScalarSingletModel` and `QuarticThermal` are toy potentials whose transitions are
too weak or too low for the thermal stage. A full run of them ends with `NoThermalParameters`, and
the reason is in `status.warnings`.

### Writing a model in Python

Subclass `pt.Potential`, call `super().__init__()`, and implement `V(phi, T)` and
`get_n_scalars()`. `phi` is a numpy array.

```python
import numpy as np
import phasetracer as pt

class QuarticThermal(pt.Potential):
    def __init__(self, D=0.1, E=0.01, lam=0.1, T0=100.):
        super().__init__()
        self.D, self.E, self.lam, self.T0 = D, E, lam, T0
    def V(self, phi, T):
        x = phi[0]
        return self.D*(T**2 - self.T0**2)*x**2 - self.E*T*x**3 + 0.25*self.lam*x**4
    def get_n_scalars(self):
        return 1
    def dV_dx(self, phi, T):          # optional; the default is a numerical derivative
        x = phi[0]
        return np.array([2*self.D*(T**2 - self.T0**2)*x - 3*self.E*T*x**2 + self.lam*x**3])

runner = pt.Runner(QuarticThermal(), config)   # the Runner keeps the model alive
```

**Optional overrides:** `forbidden(phi)`, `apply_symmetry(phi)`, `dV_dx`, `d2V_dx2` and `get_low_t_phases()`.

**One-loop models:** subclass `pt.OneLoopPotential` and implement:
- `V0(phi)` and `get_n_scalars()`;
- the field-dependent masses and their degrees of freedom: `get_scalar_masses_sq(phi, xi)` / `get_scalar_dofs()`, `get_vector_masses_sq` / `get_vector_dofs`, `get_fermion_masses_sq` / `get_fermion_dofs` (with positive dofs), and optionally the Debye masses.

PhaseTracer adds the Coleman–Weinberg and thermal corrections. Choose the daisy resummation with
`set_daisy_method(pt.DaisyMethod.Parwani)`, and the scale with `set_renormalization_scale(Q)`.

Complete examples of both kinds are in the package: `phasetracer/models/quartic_thermal.py`
(`pt.Potential`) and `phasetracer/models/dark_higgs.py` (`pt.OneLoopPotential`). Copy one as a
starting point. `examples/custom_model.py` runs both, the dark-Higgs model through to its GW spectrum.

**Where to put your model**
- **Anywhere you like.** A model is just a Python class, so define it in your own script or module and pass an instance to `pt.Runner`. This works with any installed PhaseTracer, including one from `pip install`, and needs no rebuild. `phasetracer.models` holds the models *shipped with* PhaseTracer; you don't have to add yours there.
- **Don't edit an installed copy:** the files in your environment's `site-packages/phasetracer/` are replaced by the next `pip install`.
- **Contributing a model to PhaseTracer:** in a clone, add a file to `python/phasetracer/models/`, import the class in `python/phasetracer/models/__init__.py` and add it to `__all__`. Then rebuild (`cmake --build build --target _phasetracer`) or re-run `pip install .`. A model with an optional dependency, like `TL2HDM`, is instead returned from the module-level `__getattr__` in that file, so that importing `phasetracer` never needs the dependency.

**Performance:** a Python potential is called many thousands of times, and every call holds the
Python GIL.
- Call `pt.set_num_threads(1)`, since extra OpenMP threads would only wait.
- Provide analytic derivatives where you can.
- For production scans, a C++ model is likely faster.

## Examples and tests

The examples in `examples/` only use the package. The models they run live in `phasetracer.models`.
- `simple_example.py`: the shortest complete run.
- `run_runner.py`: the full pipeline for `TwoDimModel`. Add `--plot` to save the spectrum.
- `custom_model.py`: the Python models `QuarticThermal` and `DarkHiggs`.
- `run_conformal_u1.py`: percolation parameters of `ConformalDarkU1` on the MeV scale, and how to sort status codes into successes, physical null results and errors.
- `benchmark_THDM.py`: the 2HDM benchmark scan, with either the C++ `THDM` or TransitionListener's potential (`TL2HDM`).
- `bubbledet_interface.py`: a decay-rate prefactor computed by BubbleDet. Needs BubbleDet and scipy.

The tests in `tests/` (`test_runner.py` for the interface, `test_models.py` for the models) run
with `python -m pytest python/tests`. They need `pytest`, available as `pip install .[test]`.

## For Developers
<details>
<summary>Click me</summary>
All commands below are run from the root of the repository.

### Building a local copy

There are three ways to get an importable copy of a modified tree:

| Route | Good for | After a change |
|---|---|---|
| Development build in `build/` | working on C++ and bindings | rebuild one CMake target |
| `pip install .` | a fixed snapshot in a venv | re-run `pip install .` |
| Editable install | working in a venv | nothing; it rebuilds on `import` |

**1. Development build (recommended while editing the C++).** This builds the module next to the
C++ libraries and examples:

```bash
cmake -B build -DBUILD_PYTHON=ON                 # once; add -DCMAKE_BUILD_TYPE=Release for speed
cmake --build build -j 4 --target _phasetracer   # the libraries, the module and the .py files
```

The importable package is assembled in `build/python/phasetracer`. Make Python find it with:

```bash
export PYTHONPATH=$PWD/build/python              # for this shell
```

Or make it permanent for one virtual environment by adding a `.pth` file to that environment.
Activate the environment first, and install numpy into it, since the module needs numpy:

```bash
pip install numpy
echo $PWD/build/python > $(python -c 'import site; print(site.getsitepackages()[0])')/phasetracer-dev.pth
```

Check which copy is imported with `python -c "import phasetracer; print(phasetracer.__file__)"`.
It should print `.../build/python/phasetracer/__init__.py`.

- **pybind11:** taken from an installed copy (`pip install pybind11`) if there is one, and downloaded otherwise.
- **Which Python:** the module is compiled for the Python that CMake finds. To build for a particular one, e.g. a venv's, configure with `cmake -B build -DBUILD_PYTHON=ON -DPython_EXECUTABLE=$(which python)`. A module built for a different Python version fails to import, or isn't found.
- **`make` builds the module too:** `BUILD_PYTHON` is stored in `build/CMakeCache.txt`, so a plain `cmake --build build` builds the module as well. Turn it off with `cmake -B build -DBUILD_PYTHON=OFF`.

**2. Regular install into a venv.** `pip install .` builds a separate Release copy in a temporary
directory and installs it into the active environment. It never touches `build/` or `lib/`. The
copy does not follow later edits, so re-run `pip install .` after each change. For a debug build,
use `pip install . -C cmake.define.CMAKE_BUILD_TYPE=Debug`.

**3. Editable install.** The package is rebuilt automatically whenever it is imported after a change:

```bash
pip install scikit-build-core pybind11 ninja numpy
pip install --no-build-isolation -e . -Ceditable.rebuild=true -Cbuild-dir=build/editable
```

- **Speed:** the first install takes a few minutes. After a change, the next `import phasetracer` re-runs the build, which should now run make quicker, and prints the CMake output (add `-Ceditable.verbose=false` to the install command to hide it).
- **The build directory:** `build/editable` must stay where it is; delete it and re-run the install if it gets into a bad state.

**Troubleshooting**
- **The wrong copy is imported:** check `phasetracer.__file__`. A `pip install .` in the active environment takes precedence over `PYTHONPATH`; remove it with `pip uninstall phasetracer`.
- **`ImportError: libphasetracer.so: cannot open shared object file`:** the development module finds the C++ libraries through absolute paths in `lib/` and `EffectivePotential/lib/`. Rebuild after moving the repository, or after deleting `lib/`.
- **`ImportError: ... undefined symbol`:** the module and the C++ libraries are out of step, e.g. after switching branches. Rebuild with `cmake --build build --target _phasetracer`.

### Rebuilding after changes

| You changed | Development build | `pip install .` | Editable install |
|---|---|---|---|
| A `.py` file in `python/phasetracer/` | `cmake --build build --target _phasetracer` | `pip install .` | nothing |
| C++ in `src/`, `include/` or `EffectivePotential/` | `cmake --build build --target _phasetracer`; this rebuilds the libraries first | `pip install .` | nothing |
| A binding file in `python/src/` | `cmake --build build --target _phasetracer` | `pip install .` | nothing |
| **Added** a `.cpp` file to `python/src/`, a `.py` file to `python/phasetracer/`, or a `.cpp` to `src/` | `cmake build`, then `cmake --build build --target _phasetracer`. The file lists are globbed, so CMake has to re-run to see new files | `pip install .` | re-run the `pip install -e` command |
| A CMake file or an option | `cmake build` (or `cmake -B build -D...`), then build | `pip install .` | re-run the `pip install -e` command |

- **Always restart Python afterwards** (the interpreter, script or Jupyter kernel). A compiled extension module cannot be reloaded, so `importlib.reload` keeps using the old one.
- **Test both sides after a change:**
  ```bash
  PYTHONPATH=build/python python -m pytest python/tests   # Python interface
  cmake --build build -j 4 --target unit_tests && ./bin/unit_tests   # C++ library
  ```
- **Some C++ changes need a matching binding change:**
  - A new `PROPERTY` setting in PhaseFinder, TransitionFinder, ActionCalculator, ThermoFinder or GravWaveCalculator also needs:
    - the field in the matching `Config` struct (`include/config.hpp`), with the same default;
    - a line in `apply()` and `operator<<` (`src/config.cpp`);
    - the drift check in `unit_tests/test_config.cpp`;
    - a `.def_readwrite` in `python/src/bind_config.cpp`.
  - A changed result struct (`Phase`, `Transition`, `TransitionMilestone`, `GravWaveSpectrum`, ...) needs its `def_readonly` lines in `python/src/bind_results.cpp` updated.
  - A new method on a stage class needs exposing in `python/src/bind_stages.cpp`, as a lambda that starts with `h.get()`.
  - A new `StatusCode` or `Stage` needs adding to `to_string` (`src/run_status.cpp`) and to `python/src/bind_status.cpp`.
  - A field or method that isn't bound is simply not visible from Python. Nothing fails, so the tests won't catch a missing binding.

### Adding a C++ model to `phasetracer.models`

The models are bound one family per file, in `python/src/bind_models_<family>.cpp`:
- `test`: the test models;
- `xsm`: the xSM family;
- `rss`: RSS and RSS_HT;
- `bsm`: other BSM models;
- `thdm`: the 2HDM.

They are called from `python/src/bind_models.cpp`, and re-exported by
`python/phasetracer/models/__init__.py`. Once a model works in C++, these steps make it available as
`from phasetracer.models import MyModel`.

**1. The model header.** Put it in `EffectivePotential/include/models/`, as a header-only class
deriving from `EffectivePotential::Potential` or `EffectivePotential::OneLoopPotential`, like the
existing models. That directory is already on the include path, so no CMake change is needed.
Give the header its own include guard. If you copy an existing model, rename the guard too.
Otherwise the second header included is silently skipped, and the build fails with
"`MyModel` is not a member of `EffectivePotential`". `xSM_HT.hpp` once shared its guard with
`Z2_scalar_singlet_model.hpp`, so the two could not be bound in the same file.

Suppose it looks like this:

```cpp
// EffectivePotential/include/models/my_model.hpp
#ifndef POTENTIAL_MY_MODEL_HPP_INCLUDED
#define POTENTIAL_MY_MODEL_HPP_INCLUDED
#include "one_loop_potential.hpp"
#include "property.hpp"

namespace EffectivePotential {
class MyModel : public OneLoopPotential {
public:
  MyModel(double lambda_s, double ms) : lambda_s(lambda_s), ms(ms) {}
  void init_params(double Q, bool use_daisy = true) { /* ... */ }
  double V0(Eigen::VectorXd phi) const override { /* ... */ }
  size_t get_n_scalars() const override { return 2; }
  // ...masses and dofs...
  PROPERTY(double, mass, 125.)
private:
  double lambda_s, ms;
};
}
#endif
```

**2. Bind it** in the matching `python/src/bind_models_<family>.cpp`. Include the header at the top
of the file, and add a `py::class_` to that file's `bind_models_<family>()` function:

```cpp
#include "models/my_model.hpp"
...
void bind_models_bsm(py::module_ &m) {
  ...
  using EffectivePotential::MyModel;
  py::class_<MyModel, OneLoopPotential>(m, "MyModel", "One-line description shown by help()")
      .def(py::init<double, double>(), py::arg("lambda_s"), py::arg("ms"))
      .def("init_params", &MyModel::init_params, py::arg("Q"), py::arg("use_daisy") = true)
      .def_property("mass", &MyModel::get_mass, &MyModel::set_mass);
}
```

**How each C++ feature maps to a binding:**

| C++ | Binding |
|---|---|
| Base class | The second template argument must be a class that is already bound: `Potential`, `OneLoopPotential`, or a bound intermediate class (`xSM_MSbar` and `xSM_OSlike` use `xSM_base`; `RSS_HT` uses `RSS`). An intermediate class does not have to be bound: naming `OneLoopPotential` also works. Everything the base provides (`V`, `V0`, masses, `set_daisy_method`, ...) is then inherited automatically. |
| Constructors | One `.def(py::init<argument types...>(), py::arg(...)...)` per constructor; `py::init<>()` for a default constructor. A class with pure virtual methods cannot be constructed from Python. |
| Default arguments | C++ defaults are not visible to pybind, so repeat them: `py::arg("use_daisy") = true`. |
| Methods | `.def("name", &MyModel::name, py::arg(...))`. Overloaded methods need `py::overload_cast<argument types...>(&MyModel::name)`, adding `py::const_` for const ones. |
| `PROPERTY(type, name, default)` | `.def_property("name", &MyModel::get_name, &MyModel::set_name)`. For a name that is a Python keyword (`lambda`, `from`, ...), use a trailing underscore, as `OneDimModel.lambda_` does. |
| Static factories | `.def_static("from_tadpoles", &MyModel::from_tadpoles, py::arg(...))` |
| `operator<<` | `.def("__repr__", &repr_from_stream<MyModel>)` |

Arguments and return values of type `Eigen::VectorXd`, `std::vector`, `std::map` and the bound
result types are converted automatically, to and from numpy arrays, Python lists and dicts.

**A new family file:** for a model that fits no existing family, create
`python/src/bind_models_<family>.cpp` with a `void bind_models_<family>(py::module_ &m)`:
- declare that function in `python/src/bindings.hpp`;
- call it from `bind_models()` in `python/src/bind_models.cpp`;
- re-run CMake (`cmake build`), because the source list is globbed.

Keep each model header in a single binding file.

**3. Export it from `python/phasetracer/models/__init__.py`.** Assign it with the other compiled
models, and add its name to `__all__`:

```python
MyModel = _models.MyModel

__all__ = [..., "MyModel"]
```

**4. Rebuild and check it** (restart Python first):

```bash
cmake --build build --target _phasetracer
PYTHONPATH=build/python python -c "
from phasetracer.models import MyModel
m = MyModel(lambda_s=1.0, ms=60.)
print(m.get_n_scalars(), m.V([0., 0.], 100.))"
```

**5. Add a test** to `python/tests/test_models.py`:
- add an instance to `all_models()`, which checks that it constructs, evaluates and is exported;
- add a physics check, e.g. run it with `config.pipeline.stop_after = pt.Stage.TransitionFinder` and compare the critical temperature with a value from the C++ code or the literature (see `test_xsm_ht_analytic_tc`).

**Things to watch for**
- **Headers outside `EffectivePotential/include/models/`:** move the header into that directory, as was done for `THDM.hpp` and its helpers (which used to live in `example/Benchmarks/benchmark_models/`). Don't point the bindings at `example/`: the Python package should not depend on the examples.
- **Duplicate class names:** two headers that define a class with the same name cannot both be bound (e.g. the two `THDM` headers, or `RSS.hpp` and `RSS_old.hpp`). Bind one of them.
- **Optional dependencies:** a model that needs BSMPT or FlexibleSUSY must only be included and bound inside the matching `#ifdef BUILD_WITH_...` block.
- **Keep namespace imports contained:** a model header that puts `using namespace std;` at global scope (as `THDM.hpp` does) should have a binding file of its own (`bind_models_thdm.cpp`), so it doesn't leak into other bindings. Never add such lines to `bindings.hpp`, which every binding file includes.
- **Subclassing in Python:** binding a C++ model as above is all that is needed to use it. Only if users should *subclass* the model in Python, overriding its methods, does it also need a trampoline class, like `PyOneLoopPotential` in `python/src/bind_potential.cpp`.
<!-- <details> -->
