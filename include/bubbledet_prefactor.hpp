// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef PHASETRACER_BUBBLEDET_PREFACTOR_HPP_
#define PHASETRACER_BUBBLEDET_PREFACTOR_HPP_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "effectivepotential/particle_spec.hpp"
#include "potential.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp" // not self-contained; needs the two above first

namespace PhaseTracer {

/**
 * Configuration for BubbleDetPrefactor.
 *
 * The C++ side never embeds Python: it launches `python_executable
 * worker_script` as a persistent subprocess and talks to it over pipes. This
 * keeps BubbleDet's NumPy/SciPy stack out of the PhaseTracer process and,
 * crucially, keeps the GIL off PhaseTracer's OpenMP action loop.
 */
struct BubbleDetConfig {
  /** Python interpreter that can import BubbleDet (e.g. a venv's python). */
  std::string python_executable = "python3";

  /** Absolute path to bubbledet_worker.py. Required. */
  std::string worker_script;

  /** Spacetime dimension of the bounce (3 for the thermal O(3) bounce). */
  int dim = 3;

  /** Pass BubbleDet's thermal dynamical prefactor. */
  bool thermal = true;

  /** Number of points used to sample V and the masses over the field range. */
  int phi_grid_points = 200;

  /** Fractional extension of the [min(Phi), max(Phi)] sampling range. */
  double phi_grid_margin = 0.05;

  /**
   * Degrees of freedom of the longitudinal (Higgs) fluctuation that the functor
   * always emits with W = d2V along the path. The bounce field itself is one
   * scalar dof.
   */
  double higgs_dof = 1.0;

  /**
   * Number of transverse CP-even scalar modes to add, each with
   * W = n_hat^T (d2V/dphi^2) n_hat (the effective-potential Hessian projected
   * onto a path-normal direction). -1 means "auto" = n_scalars - 1. Set 0 to
   * disable. Only the 2-field case (one transverse mode) is implemented.
   */
  int transverse_scalar_modes = -1;

  /** Degrees of freedom assigned to each transverse CP-even mode. */
  double transverse_dof = 1.0;

  /**
   * Optional override of the fluctuation spectrum. If unset, the potential's
   * own get_fluctuation_spectrum(T) is used. Handy for single-field models that
   * derive only from Potential (no built-in spectrum).
   */
  std::function<std::vector<EffectivePotential::ParticleSpec>(double T)> spectrum_provider;
};

/**
 * A FalseVacuumDecayRate prefactor that evaluates the full one-loop functional
 * determinant via BubbleDet, returning A = exp(-S1).
 *
 * Copyable and storable in a std::function (state lives behind a shared_ptr),
 * so it can be installed with FalseVacuumDecayRate::set_prefactor_function or
 * passed to a constructor. Thread-safe: concurrent calls (from the OpenMP
 * action loop) are serialised onto the single worker.
 *
 * If BubbleDet is unavailable, the spectrum is empty, or a determinant fails,
 * it transparently falls back to the analytic prefactor
 * (FalseVacuumDecayRate::default_decay_rate_prefactor) and logs a warning, so a
 * pipeline never breaks because BubbleDet hiccupped on one temperature.
 */
class BubbleDetPrefactor {
public:
  BubbleDetPrefactor(const EffectivePotential::Potential &potential, BubbleDetConfig config);

  /** Matches FalseVacuumDecayRate::PrefactorFunction. */
  double operator()(double temperature, double action_on_T, const ActionResult &bounce) const;

private:
  class Impl;
  std::shared_ptr<Impl> impl_;
};

} // namespace PhaseTracer

#endif // PHASETRACER_BUBBLEDET_PREFACTOR_HPP_
