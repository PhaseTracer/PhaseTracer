# PhaseTracer ↔ BubbleDet validation harness (step 6.1)

A throwaway-friendly scaffold that de-risks the BubbleDet interface **before**
any in-C++ plumbing is committed. It checks that PhaseTracer's bounce action and
profile mean the same thing BubbleDet/CosmoTransitions expect, using the 1D
analytic test model (exact potential → any disagreement is in the interface, not
the physics).

## Pieces

| File | Role |
|------|------|
| [bounce_bridge.cpp](bounce_bridge.cpp) | pybind11 module `pt_bounce`: runs the real PhaseTracer pipeline (PhaseFinder → TransitionFinder → ActionCalculator) for `OneDimModel` and exposes the bounce `action`, profile `(R, Φ, dΦ)`, the vacua, and `V/dV/d2V` to Python. |
| [validate_bubbledet.py](validate_bubbledet.py) | The harness. Cross-checks the PhaseTracer action against CosmoTransitions on the *same* C++ potential, then (if BubbleDet is importable) runs the one-loop determinant on PhaseTracer's own profile. |

## Build & run

The core libs must be built **shared** first (the main project build does this):
`libeffectivepotential.so`, `libphasetracer.so`.

```bash
cd example/BubbleDetBridge
cmake -S . -B build
cmake --build build -j
python3 validate_bubbledet.py        # finds build/pt_bounce automatically
```

## What it proves

1. **Action mapping (always runs — needs only CosmoTransitions).**
   `S_PT` vs `S_CT` for `dim = 3` (O(3) thermal bounce `S₃`, CosmoTransitions
   `alpha = dim-1`). Agreement ⇒ PhaseTracer's `S` and `S/T` are the same
   quantity in the same dimension/units BubbleDet assumes. Current agreement is
   ~10⁻⁴–10⁻³ relative (degrades mildly at large `S/T`, i.e. thick barriers —
   a profile-resolution effect, not a units mismatch).

2. **Prefactor (runs if `BubbleDet` is importable).**
   Feeds PhaseTracer's profile into `BubbleConfig.fromCosmoTransitions` with a
   single Higgs scalar (`W_Phi = d²V/dφ²`, `spin=0`, `dof_internal=1`,
   `zero_modes="Higgs"`), `thermal=True`, and reports `S1 = −log A`. It also
   runs BubbleDet on CosmoTransitions' own profile as a control, so the two `S1`
   values can be compared. This is exactly what the future C++ prefactor functor
   (installed via `FalseVacuumDecayRate::set_prefactor_function`) will compute.

> BubbleDet is not bundled. Use a Python env with it installed (here:
> `~/.venvs/phasetracer_venv`); the harness skips step 2 cleanly otherwise.

## Findings from running it (validated)

- **Action mapping holds:** `S_PT` vs `S_CT` agree to ~10⁻⁴–10⁻³ across the
  temperature window. Units/dimension are correct.
- **BubbleDet requires the profile to include `r = 0`.** It evaluates the
  fluctuation operator at the origin, but PhaseTracer's grid starts at a small
  finite `rmin` (PhaseTracer's `R[0]` is exactly CosmoTransitions' `R[1]`).
  The fix is to **prepend `(R=0, Φ=Φ(rmin), dΦ=0)`** (spherical regularity).
  *This is a hard requirement the production C++ bridge must satisfy* — either
  prepend the origin when handing the profile to BubbleDet, or have
  `ActionCalculator` emit an origin-inclusive profile.
- **Prefactor agrees within profile error:** at `T ≈ 46.16`,
  `S1(PT profile) = −2.99 ± 0.036` vs `S1(CT profile) = −2.91 ± 0.0003`
  (`A = e^{−S1} ≈ 19.8`). Central values agree to ~2.5%; PhaseTracer's larger
  quoted error comes from its coarser/less-smooth grid near the origin. For
  production precision, consider a denser PhaseTracer profile (raise
  `PD_npoints` / lower `PD_rmin`) when BubbleDet is in use.

## Production interface (built)

The validated mapping is now implemented as a drop-in prefactor for
`FalseVacuumDecayRate`, with **no embedded Python** — BubbleDet runs in a
persistent subprocess and only sampled arrays cross the boundary.

| Piece | Role |
|-------|------|
| [bubbledet_worker.py](bubbledet_worker.py) | Persistent worker. Imports BubbleDet once, then serves one bounce per stdin line (JSON), replying one line per response. Builds `BubbleConfig`/`ParticleConfig` from C++-sampled `V/dV/d2V/W` splines. |
| [bubbledet_prefactor.hpp](../../include/bubbledet_prefactor.hpp) / [.cpp](../../src/bubbledet_prefactor.cpp) | `PhaseTracer::BubbleDetPrefactor` — a copyable functor matching `FalseVacuumDecayRate::PrefactorFunction`. Spawns/owns the worker (POSIX pipes), prepends the `r=0` point, samples `V` and each particle's `W(φ)` over the field range, sources particles from `Potential::get_fluctuation_spectrum(T)` (or a `spectrum_provider` override), tags the tunneling scalar as the Higgs zero mode, and returns `A = exp(−S1)`. Thread-safe (mutex over the pipe), and **falls back to the analytic prefactor** on any failure. |
| [run_bubbledet_prefactor.cpp](../run_bubbledet_prefactor.cpp) | End-to-end test: installs `BubbleDetPrefactor` on a `FalseVacuumDecayRate` for the 1D model and compares `get_gamma` against analytic. |

Usage:

```cpp
PhaseTracer::BubbleDetConfig cfg;
cfg.python_executable = "/home/.../phasetracer_venv/bin/python";
cfg.worker_script     = "/abs/path/example/BubbleDetBridge/bubbledet_worker.py";
PhaseTracer::BubbleDetPrefactor A(potential, cfg);
PhaseTracer::FalseVacuumDecayRate rate(trans, ac, t_min, t_max, n_knots, A);
// rate.get_gamma(T) now uses BubbleDet's one-loop prefactor.
```

For an `OneLoopPotential`, `cfg.spectrum_provider` is unnecessary — the
spectrum is auto-synthesised. The bare 1D test model has no spectrum, so the
test supplies a single Higgs scalar (`W = d²V/dφ²`) via `spectrum_provider`.

Verified end-to-end (`bin/run_bubbledet_prefactor`):
```
       T          gamma_BD    gamma_analytic           ratio
    43.0      8.882e-01       3.899e+04        2.28e-05
    46.0      9.400e-02       1.158e+04        8.12e-06
```
Ratio ≪ 1 (not 1) ⇒ BubbleDet's real determinant `A ≈ 20` replaces the
analytic dimensional estimate `T⁴(S/2π)^{3/2} ≈ 4×10⁶`. With a missing/broken
worker the ratio is exactly 1 (clean fallback, no crash).

> **Scope:** single-field bounces (`n_scalars == 1`) for now; multifield needs
> the field-space path projected onto BubbleDet's 1D `W(φ)`, which is deferred.

## How this maps onto the real interface

- `pt_bounce` is a stand-in for what the production bridge will expose from the
  `ActionCalculator` / `FalseVacuumDecayRate` pipeline.
- The single Higgs `ParticleConfig` here is hand-built. In production the
  particle list comes from `Potential::get_fluctuation_spectrum(T)`
  (auto-synthesised by `OneLoopPotential`), with the tunneling field tagged
  `ZeroModeType::Higgs` by the decay-rate driver.
- The `S1` produced here is `−log A`; in `FalseVacuumDecayRate::get_splines`
  that becomes `log(prefactor) = −S1`, i.e. `log_gamma = −S1 − S`.
