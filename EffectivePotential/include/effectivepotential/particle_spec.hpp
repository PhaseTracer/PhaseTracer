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

#ifndef EFFECTIVEPOTENTIAL_PARTICLE_SPEC_HPP_
#define EFFECTIVEPOTENTIAL_PARTICLE_SPEC_HPP_

#include <functional>
#include <string>

#include <Eigen/Core>

namespace EffectivePotential {

/** Quantum statistics of a fluctuating field. */
enum class FieldStatistics { Boson,
                             Fermion };

/**
 * Treatment of zero modes when evaluating the one-loop fluctuation determinant
 * around the bounce. Mirrors BubbleDet's `ParticleConfig.zero_modes`:
 *  - Higgs:     the field providing the translational (and negative) modes
 *               along the tunneling direction. Exactly one species per bounce.
 *  - Goldstone: a field whose mass vanishes in the broken phase.
 *  - None:      no special zero-mode treatment.
 */
enum class ZeroModeType { None,
                          Higgs,
                          Goldstone };

/**
 * A single fluctuating field species participating in the one-loop bounce
 * determinant. Every quantity here is also required for the Coleman-Weinberg /
 * thermal loop corrections, so OneLoopPotential can synthesise the spectrum
 * automatically from its existing mass and dof accessors.
 *
 * This is the C++ analogue of BubbleDet's `ParticleConfig`.
 */
struct ParticleSpec {
  /** Human-readable label (e.g. "Higgs", "W", "top"). */
  std::string name;

  /**
   * Field-dependent mass-squared evaluated on the bounce background, m^2(phi).
   * Supplied as a closure so any temperature or gauge dependence is already
   * baked in at the point the spectrum is requested. For scalars this is an
   * eigenvalue of the field-space Hessian d2V/dphi^2.
   */
  std::function<double(const Eigen::VectorXd &phi)> mass_sq;

  /** Spin of the field: 0, 0.5 or 1. */
  double spin = 0.0;

  /** Internal degrees of freedom (BubbleDet `dof_internal`). */
  double dof = 1.0;

  /** Quantum statistics. */
  FieldStatistics statistics = FieldStatistics::Boson;

  /** Zero-mode treatment for the fluctuation determinant. */
  ZeroModeType zero_mode = ZeroModeType::None;
};

} // namespace EffectivePotential

#endif // EFFECTIVEPOTENTIAL_PARTICLE_SPEC_HPP_
