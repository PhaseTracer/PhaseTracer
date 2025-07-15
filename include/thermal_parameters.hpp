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

#ifndef PHASETRACER_THERMAL_PARAMETERS_HPP_
#define PHASETRACER_THERMAL_PARAMETERS_HPP_

#include <cmath>
#include <vector>
#include <optional>
#include <stdexcept>
#include <interpolation.h>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {

/**
 * @class Thermodynamics
 * @brief Computes the thermodynamic quantities p, e, w, s for a given phase
 */
class Thermodynamics {

private:
  Phase phase;
  std::vector<double> potential_values;
  std::vector<double> temperature_values;
  alglib::spline1dinterpolant potential;

  double t_min, t_max;

  alglib::spline1dinterpolant p;
  alglib::spline1dinterpolant e;
  alglib::spline1dinterpolant w;
  alglib::spline1dinterpolant s;

  std::vector<double> temperature;
  std::vector<double> pressure;
  std::vector<double> entropy;
  std::vector<double> energy;
  std::vector<double> enthalpy;

  int n_temp;
  double dof;

public:
  Thermodynamics(Phase phase_in, int n_temp_in, double dof_in);

  void get_thermodynamic_splines();

	/**
	 * @brief Computes the sound speed.
	 * @param T Temperature at which to evaluate
	 * @return Sound speed at T.
	 */
  double get_cs(double T) const;

	/**
	 * @brief Computes the trace anomoly.
	 * @param T Temperature at which to evaluate.
	 * @return Trace anomoly at T.
	 */
  double get_theta(double T) const;

	/**
	 * @brief Computes the pressure density in the phase.
	 * @param T Temperature at which to evaluate.
	 * @return Pressure density at T.
	 */
  double get_pressure(double T) const;

	/**
	 * @brief Computes the energy density in the phase.
	 * @param T Temperature at which to evaluate.
	 * @return Energy density at T.
	 */
  double get_energy(double T) const;

	/**
	 * @brief Computes the enthalpy density in the phase.
	 * @param T Temperature at which to evaluate.
	 * @return Enthalpy density at T.
	 */
  double get_enthalpy(double T) const;

	/**
	 * @brief Computes the entropy density in the phase.
	 * @param T Temperature at which to evaluate.
	 * @return Entropy density at T.
	 */
  double get_entropy(double T) const;

private:
  /**
   * @brief Checks if temperature is within valid range
   * @param T Temperature to check
   * @param caller Name of calling function for error message
   * @throws std::out_of_range if temperature is out of bounds
   */
  void check_temperature_range(double T, const char* caller) const {
    if (T < t_min || T > t_max) {
      throw std::out_of_range(std::string("Temperature out of interpolation bounds in ") + caller);
    }
  }
};

/**
 * @class Bounce
 * @brief Evaluates bounce action
 */
class Bounce {

private:
  TransitionFinder tf;
	Transition t;
  int spline_evaluations;

public:
	double minimum_temp, maximum_temp;
	alglib::spline1dinterpolant action_spline, gamma_spline;

	Bounce();
	Bounce(Transition t_in, TransitionFinder tf_in, double minimum_temp_in, double maximum_temp_in, int spline_evaluations_in);

	/**
	 * @brief Computes the bounce and gamma splines for given transition.
	 * @return void.
	 */
	void get_splines();

private:
  /** Calculate action at given temperature */
  std::optional<double> calculate_action(double temperature) const {
    double action = tf.get_action(t.true_phase, t.false_phase, temperature) / temperature;
    if (std::isnan(action) || std::isinf(action) || action > 1e150) {
      return std::nullopt;
    }
    return action;
  }

  /** Calculate gamma at given temperature and action */
  std::optional<double> calculate_gamma(double temperature, double action) const {
    double log_gamma = 4*log(temperature) + 1.5*log(action/(2.*M_PI)) - action;
    return log_gamma < -700 ? std::nullopt : std::optional<double>(exp(log_gamma));
  }
};

/**
 * @struct Contains all thermal_parameters
 */
struct ThermalParams {
  size_t key;
  std::vector<double> pressure_density_true, pressure_density_false;
  std::vector<double> energy_density_true, energy_density_false;
	alglib::spline1dinterpolant hubble_spline;
	double TC;
  double TN;
	double TP;
	double TF;
	double alpha_tp, alpha_tn;
	double betaH_tp, betaH_tn;
	double beta_tp, beta_tn;
	double H_tp, H_tn;
  bool percolates;
  bool nucleates;
  bool completes;

  /** Pretty-printer for single phase */
  friend std::ostream &operator<<(std::ostream &o, const ThermalParams &tp) {
    o << "=== transition @ TC = " << tp.TC << " ===" << std::endl;
    if(tp.percolates) {
      o << "percolation temperature = " << tp.TP << std::endl
      << "alpha (TP) = " << tp.alpha_tp << std::endl
      << "beta (TP) = " << tp.beta_tp << std::endl
      << "beta/H (TP) = " << tp.betaH_tp << std::endl
      << "H (TP) = " << tp.H_tp << std::endl;
    } else {
      o << "transition does not percolate." << std::endl;
    }
    if(tp.nucleates) {
      o << "nucleation temperature = " << tp.TN << std::endl
      << "alpha (TN) = " << tp.alpha_tn << std::endl
      << "beta (TN) = " << tp.beta_tn << std::endl
      << "beta/H (TN) = " << tp.betaH_tn << std::endl
      << "H (TN) = " << tp.H_tn << std::endl;
    } else {
      o << "transition does not nucleate." << std::endl;
    }
    if(tp.completes) {
      o << "completion temperature = " << tp.TF << std::endl;
    } else {
      o << "transition does not complete." << std::endl;
    }
    return o;
  }

};

/**
 * @class ThermalParameters
 * @brief Performs detailed calculations of thermal parameters and EOS
 */
class ThermalParameters {

private:
  TransitionFinder tf;
  std::vector<ThermalParams> thermalparameterContainer;

  bool calculated_thermal_params = false;

  /** number of temperature values for thermo splines  **/
  PROPERTY(int, n_temp, 250);

  /** relativistic degrees of freedom  **/
  PROPERTY(double, dof, 106.75);

  /** number of temp values for the action spline */
  PROPERTY(int, spline_evaluations, 50);

  /** uses percolation temp if true, otherwise nucleation temp */
  PROPERTY(bool, use_percolation_temp, false);

  /** minimum temperature interval a phase must exist. */
  PROPERTY(double, dt_tol, 1.0);

  /** Newtons constant in GeV  */
  PROPERTY(double, G, 1/((1.22 * 1e19)*(1.22 * 1e19)));

public:

  /** Retrieve all transitions between all phases */
  std::vector<ThermalParams> get_thermal_params() const { return thermalparameterContainer; }

  friend std::ostream &operator<<(std::ostream &o, const ThermalParameters &tp) {
    auto thermal_params = tp.thermalparameterContainer;

    o << "found " << thermal_params.size() << " thermal param set";
    if (thermal_params.size() != 1) {
      o << "s";
    }
    o << std::endl
      << std::endl;

    for (auto &p : thermal_params) {
      o << p << std::endl;
    }

    return o;
  }

  ThermalParameters(TransitionFinder tf_in) : tf(tf_in) {}

  void find_thermal_parameters();

  std::vector<ThermalParams> get_thermal_parameters();

private:
  /**
   * @brief Calculates alpha parameter
   * @param T Temperature at which to evaluate
   * @param true_thermo True vacuum phase thermodynamics
   * @param false_thermo False vacuum phase thermodynamics
   * @return Alpha parameter at T
   */
  double get_alpha(double T,
                  Thermodynamics true_thermo,
                  Thermodynamics false_thermo);

  double get_betaH(double T,
                  Bounce bounce);
  /**
   * @brief Calculates the Hubble rate at a given temperature
   * @param T Temperature at which to evaluate
   * @param true_phase True vacuum phase thermodynamics
   * @param false_phase False vacuum phase thermodynamics
   * @return Hubble rate at T
   */
  double get_hubble_rate(double T,
                        Thermodynamics true_thermo, 
                        Thermodynamics false_thermo);

  /**
   * @brief Creates a spline for the Hubble parameter
   * @param hubble_spline Output parameter for the spline
   * @param true_phase True vacuum phase thermodynamics
   * @param false_phase False vacuum phase thermodynamics
   * @param t_min Minimum temperature
   * @param t_max Maximum temperature
   * @param n_temp Number of temperature points for spline
   */
  void make_hubble_spline(alglib::spline1dinterpolant& hubble_spline, 
                          Thermodynamics true_phase, 
                          Thermodynamics false_phase, 
                          double t_min, double t_max, 
                          double n_temp = 50);

  /**
   * @brief Calculates the integral of the Hubble parameter
   * @param hubble_spline The Hubble parameter spline
   * @param T Lower temperature bound
   * @param Tdash Upper temperature bound
   * @return The integral value
   */
  double hubble_integral(alglib::spline1dinterpolant& hubble_spline, double T, double Tdash);

  /**
   * @brief Calculates the integrand for the false vacuum fraction
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @param T Temperature
   * @param Tdash Integration variable temperature
   * @return The integrand value
   */
  double false_vacuum_fraction_integrand(alglib::spline1dinterpolant& hubble_spline, 
                                        const Bounce& bounce, double T, double Tdash);

  /**
   * @brief Calculates the false vacuum fraction
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @param T Temperature
   * @param vw Wall velocity
   * @return The false vacuum fraction
   */
  double false_vacuum_fraction(alglib::spline1dinterpolant& hubble_spline, 
                              const Bounce& bounce, double T, double vw);

  /**
   * @brief Helper function for binary search to find critical temperatures
   * @param init_T Initial (upper) temperature
   * @param end_T End (lower) temperature
   * @param tol_rel Relative tolerance for convergence
   * @param calc_value Function that calculates the value to match against target
   * @param target Target value to find
   * @param max_iterations Maximum number of iterations
   * @return The found temperature
   * @throws std::runtime_error if the search fails to converge
   */
  static double find_temperature_binary_search(double init_T, double end_T, 
                                              double tol_rel,
                                              const std::function<double(double)>& calc_value,
                                              double target, 
                                              int max_iterations = 1000);

  /**
   * @brief Gets the percolation temperature
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @param vw Wall velocity (must be between 0 and 1)
   * @return The percolation temperature
   * @throws std::runtime_error if the transition does not percolate or parameters are invalid
   */
  double get_percolation_temperature(alglib::spline1dinterpolant& hubble_spline, 
                                    const Bounce& bounce, double vw);

  /**
   * @brief Gets the completion temperature
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @param vw Wall velocity (must be between 0 and 1)
   * @return The completion temperature
   * @throws std::runtime_error if the transition does not complete or parameters are invalid
   */
  double get_completion_temperature(alglib::spline1dinterpolant& hubble_spline, 
                                  const Bounce& bounce, double vw);

  /**
   * @brief Calculates the nucleation rate integrand
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @param T Temperature
   * @return The integrand value
   */
  double nucleation_rate_integrand(alglib::spline1dinterpolant& hubble_spline, 
                                  const Bounce& bounce, double T);

  /**
   * @brief Calculates the nucleation rate
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @param T Temperature
   * @return The nucleation rate
   */
  double nucleation_rate(alglib::spline1dinterpolant& hubble_spline, 
                        const Bounce& bounce, double T);

  /**
   * @brief Gets the nucleation temperature
   * @param hubble_spline The Hubble parameter spline
   * @param bounce The bounce solution
   * @return The nucleation temperature
   * @throws std::runtime_error if the temperature range is invalid or calculation fails
   */
  double get_nucleation_temperature(alglib::spline1dinterpolant& hubble_spline, 
                                  const Bounce& bounce);
};

} // namespace PhaseTracer

#endif // PHASETRACER_THERMAL_PARAMETERS_HPP_