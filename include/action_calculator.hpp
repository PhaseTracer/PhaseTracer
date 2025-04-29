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

#ifndef PHASETRACER_ACTION_CALCULATOR_HPP_
#define PHASETRACER_ACTION_CALCULATOR_HPP_

#include <algorithm>
#include <limits>
#include <memory>
#include <ostream>
#include <vector>

#ifdef BUILD_WITH_BP
#include "bubble_profiler.hpp"
#endif
#include "path_deformation.hpp"

namespace PhaseTracer {

/** Method for action calculation */
enum class ActionMethod { None,
                          BubbleProfiler,
                          PathDeformation,
                          All };

class ActionCalculator {

private:
  EffectivePotential::Potential &potential;

  /** Number of dimensions */
  PROPERTY(size_t, num_dims, 3)

  // Save the profile and path
  mutable Profile1D bubble_profile;
  mutable std::vector<Eigen::VectorXd> tunneling_path;

  /** Choose method to calculate the action */
  ActionMethod action_method = ActionMethod::PathDeformation;

  /** Set parameters in BubbleProfiler */
  /** Use perturbative method or not*/
  PROPERTY(bool, BP_use_perturbative, false)
  PROPERTY(double, BP_initial_step_size, 1.e-2)
  PROPERTY(double, BP_interpolation_points_fraction, 1.0)
  // TODD: Add more settings if needed

  /** Set parameters for PathDeformation */
  /** The precision of field values after taking the logarithm */
  PROPERTY(double, PD_xtol, 1e-4)
  /** The fractional error tolerance in integration*/
  PROPERTY(double, PD_phitol, 1e-4)
  /** The cut off for finding the initial conditions for integration */
  PROPERTY(double, PD_thin_cutoff, .01)
  /** Number of points to return in the profile */
  PROPERTY(double, PD_npoints, 500)
  /** The smallest starting radius */
  PROPERTY(double, PD_rmin, 1e-4)
  /** The maximum allowed integration distance */
  PROPERTY(double, PD_rmax, 1e4)
  /** The maximum number of points to be positioned during the integration process. */
  PROPERTY(boost::uintmax_t, PD_max_iter, 100)
  /** Number of basis splines to use */
  PROPERTY(size_t, PD_nb, 10);
  /** Order of basis splines */
  PROPERTY(size_t, PD_kb, 3);
  /** Get each step saved */
  PROPERTY(bool, PD_save_all_steps, false);
  /** The smallest the square of dphidr is allowed to be */
  PROPERTY(double, PD_v2min, 0.0);
  /** Maximum number of steps to take in a deformation */
  PROPERTY(size_t, PD_step_maxiter, 500);
  /** Maximum number of allowed deformation iterations */
  PROPERTY(size_t, PD_path_maxiter, 20);
  /** Number of samples to take along the path to create the spline
   interpolation functions */
  PROPERTY(size_t, PD_V_spline_samples, 100);
  /** Flag to extend the path to minimums*/
  PROPERTY(bool, PD_extend_to_minima, true);

public:
  explicit ActionCalculator(EffectivePotential::Potential &potential_) : potential(potential_) {
  }

  void set_action_calculator(ActionMethod am) {
    action_method = am;

#ifndef BUILD_WITH_BP
    if (action_method == ActionMethod::BubbleProfiler || action_method == ActionMethod::All) {
      LOG(fatal) << "Enable BubbleProfiler in CMake configuration before using it.";
      throw std::runtime_error("BubbleProfiler is not installed.");
    }
#endif
  }

  ActionMethod get_action_calculator() const { return action_method; }

  double get_action(Eigen::VectorXd true_vacuum, Eigen::VectorXd false_vacuum, double T) const {

    if (potential.V(true_vacuum, T) > potential.V(false_vacuum, T))
      true_vacuum.swap(false_vacuum);

#ifdef BUILD_WITH_BP
    double action_BP = std::numeric_limits<double>::quiet_NaN();
    if (action_method == ActionMethod::BubbleProfiler || action_method == ActionMethod::All) {
      V_BubbleProfiler V_BP(potential); // perturbative_profiler only accept non-const potential
      V_BP.set_T(T);

      if (potential.get_n_scalars() == 1 && !BP_use_perturbative) {

        double false_min = false_vacuum[0];
        double true_min = true_vacuum[0];
        auto barrier_ = V_BP.find_one_dimensional_barrier(true_vacuum, false_vacuum, T);
        double barrier = barrier_[0];

        LOG(debug) << "Calculate action (BP) at T=" << T << ", with false, true, barrier = " << false_min << ", " << true_min << ", " << barrier;

        try {
          BubbleProfiler::Shooting one_dim;
          one_dim.solve(V_BP, false_min,
                        true_min, barrier, num_dims, BubbleProfiler::Shooting::Solver_options::Compute_action);
          action_BP = one_dim.get_euclidean_action();
        } catch (const std::exception &e) {
          LOG(warning) << "At T=" << T << ", between[" << false_min << "] and [" << true_min << "]: " << e.what();
        }
      } else {
        LOG(debug) << "Calculate action (BP) at T = " << T << ", between " << false_vacuum << " and " << true_vacuum;
        BubbleProfiler::RK4_perturbative_profiler profiler;

        profiler.set_initial_step_size(BP_initial_step_size);
        profiler.set_interpolation_points_fraction(BP_interpolation_points_fraction);

        profiler.set_false_vacuum_loc(false_vacuum);
        profiler.set_true_vacuum_loc(true_vacuum);
        profiler.set_number_of_dimensions(num_dims);
        auto root_finder = std::make_shared<BubbleProfiler::GSL_root_finder<Eigen::Dynamic>>();
        profiler.set_root_finder(root_finder);
        std::shared_ptr<BubbleProfiler::Profile_guesser> guesser;
        guesser = std::make_shared<BubbleProfiler::Kink_profile_guesser>();
        profiler.set_initial_guesser(guesser);
        auto convergence_tester = std::make_shared<BubbleProfiler::Relative_convergence_tester>(
            1.e-3, 1.e-3);
        profiler.set_convergence_tester(convergence_tester);

        try {
          profiler.calculate_bubble_profile(V_BP);
          action_BP = profiler.get_euclidean_action();
        } catch (const std::exception &e) {
          LOG(warning) << "At T = " << T << ", between " << false_vacuum << " and " << true_vacuum << ": " << e.what();
        }
      }
      LOG(debug) << "S(BP) = " << action_BP;
    }
#endif
    double action_PD = std::numeric_limits<double>::quiet_NaN();
    if (action_method == ActionMethod::PathDeformation || action_method == ActionMethod::All) {
      LOG(debug) << "Calculate action (PD) at T = " << T << ", between " << false_vacuum << " and " << true_vacuum;
      if (potential.get_n_scalars() == 1) {
        OneDimPotentialForShooting ps(potential);
        ps.set_T(T);
        Shooting st(ps, num_dims - 1);

        st.set_xtol(PD_xtol);
        st.set_phitol(PD_phitol);
        st.set_thin_cutoff(PD_thin_cutoff);
        st.set_rmin(PD_rmin);
        st.set_rmax(PD_rmax);
        st.set_max_iter(PD_max_iter);

        try {
          bubble_profile = st.findProfile(false_vacuum[0], true_vacuum[0]);
          action_PD = st.calAction(bubble_profile);
        } catch (const std::exception &e) {
          LOG(warning) << "At T = " << T << ", between " << false_vacuum << " and " << true_vacuum << ": " << e.what();
        }
      } else {
        PathDeformation pd(potential);

        pd.set_nb(PD_nb);
        pd.set_kb(PD_kb);
        pd.set_save_all_steps(PD_save_all_steps);
        pd.set_v2min(PD_v2min);
        pd.set_step_maxiter(PD_step_maxiter);
        pd.set_path_maxiter(PD_path_maxiter);
        pd.set_V_spline_samples(PD_V_spline_samples);
        pd.set_extend_to_minima(PD_extend_to_minima);

        /** Pass through the shooting settings */
        pd.set_xtol(PD_xtol);
        pd.set_phitol(PD_phitol);
        pd.set_thin_cutoff(PD_thin_cutoff);
        pd.set_rmin(PD_rmin);
        pd.set_rmax(PD_rmax);
        pd.set_max_iter(PD_max_iter);

        pd.set_T(T);
        std::vector<Eigen::VectorXd> path_pts;
        path_pts.push_back(true_vacuum);
        path_pts.push_back(false_vacuum);
        try {
          FullTunneling full_tunneling = pd.full_tunneling(path_pts);
          bubble_profile = full_tunneling.profile1D;
          tunneling_path = full_tunneling.phi;
        } catch (const std::exception &e) {
          LOG(warning) << "At T = " << T << ", between " << false_vacuum << " and " << true_vacuum << ": " << e.what();
        }
        action_PD = pd.get_action();
      }
      LOG(debug) << "S(PD) = " << action_PD;
    }

#ifdef BUILD_WITH_BP
    if (std::isnan(action_BP)) {
      return action_PD;
    }
    if (std::isnan(action_PD)) {
      return action_BP;
    }
    return std::min(action_BP, action_PD);
#else
    return action_PD;
#endif
  }
  Profile1D get_bubble_profile() {
    return bubble_profile;
  }
  std::vector<Eigen::VectorXd> get_tunneling_path() {
    return tunneling_path;
  }
};

} // namespace PhaseTracer

#endif // PHASETRACER_ACTION_CALCULATOR_HPP_
