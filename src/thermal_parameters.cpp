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

#include <cmath>
#include "logger.hpp"
#include "thermal_parameters.hpp"

namespace PhaseTracer {

/*======================================
	Below are all definitions pertaining to
	the Thermodynamics class
	======================================*/

Thermodynamics::Thermodynamics(Phase phase_in, int n_temp_in, double dof_in)
  : phase(phase_in),
	potential_values(phase.V),
	temperature_values(phase.T),
	t_min(phase.T.front()),
	t_max(phase.T.back()),
	n_temp(n_temp_in),
	dof(dof_in) {

  alglib::real_1d_array t, v;
  v.setcontent(potential_values.size(), potential_values.data());
  t.setcontent(temperature_values.size(), temperature_values.data());
  alglib::spline1dbuildcubic(t, v, this->potential);

  temperature.resize(n_temp);
  double dT = std::abs(t_max - t_min) / (n_temp - 1);
  for (int i = 0; i < n_temp; ++i) {
	temperature[i] = t_min + i * dT;
  }

	if (temperature.size() != n_temp) {
		throw std::runtime_error("Temperature vector size does not match n_temp.");
	}

	get_thermodynamic_splines();
}

void Thermodynamics::get_thermodynamic_splines() {

  pressure.resize(n_temp);
  energy.resize(n_temp);
  enthalpy.resize(n_temp);
  entropy.resize(n_temp);

  for (int i = 0; i < n_temp; ++i) {
	double temp = temperature[i];
	double v, dvdT, ddvdT;
	alglib::spline1ddiff(this->potential, temp, v, dvdT, ddvdT);

	double temp2 = temp * temp;
	double temp4 = temp2 * temp2;
	double pres_bag = dof * M_PI * M_PI / 90.0 * temp4;

	pressure[i] = pres_bag - v;
	energy[i] = 3.0 * pres_bag + v - temp * dvdT;
	enthalpy[i] = 4.0 * pres_bag - temp * dvdT;
	entropy[i] = 4.0 * pres_bag / temp - dvdT;
  }

  alglib::real_1d_array t_array, p_array, e_array, w_array, s_array;
  t_array.setcontent(n_temp, temperature.data());
  p_array.setcontent(n_temp, pressure.data());
  e_array.setcontent(n_temp, energy.data());
  w_array.setcontent(n_temp, enthalpy.data());
  s_array.setcontent(n_temp, entropy.data());

  alglib::spline1dbuildcubic(t_array, p_array, this->p);
  alglib::spline1dbuildcubic(t_array, e_array, this->e);
  alglib::spline1dbuildcubic(t_array, w_array, this->w);
  alglib::spline1dbuildcubic(t_array, s_array, this->s);
}

double Thermodynamics::get_cs(double T) const {
  check_temperature_range(T, "get_cs");
  
  double pp, dp, ddp;
  double ee, de, dde;

  alglib::spline1ddiff(this->p, T, pp, dp, ddp);
  alglib::spline1ddiff(this->e, T, ee, de, dde);

  return sqrt(dp / de);
}

double Thermodynamics::get_theta(double T) const {
  check_temperature_range(T, "get_theta");

  double ee = alglib::spline1dcalc(this->e, T);
  double pp = alglib::spline1dcalc(this->p, T);

  return (ee - 3.0 * pp) / 4.0;
}

double Thermodynamics::get_pressure(double T) const {
  check_temperature_range(T, "get_pressure");
  return alglib::spline1dcalc(this->p, T);
}

double Thermodynamics::get_energy(double T) const {
  check_temperature_range(T, "get_energy");
  return alglib::spline1dcalc(this->e, T);
}

double Thermodynamics::get_enthalpy(double T) const {
  check_temperature_range(T, "get_enthalpy");
  return alglib::spline1dcalc(this->w, T);
}

double Thermodynamics::get_entropy(double T) const {
  check_temperature_range(T, "get_entropy");
  return alglib::spline1dcalc(this->s, T);
}

/*======================================
	Below are all definitions pertaining to
	the Bounce class
	======================================*/

Bounce::Bounce(Transition t_in, TransitionFinder tf_in, double minimum_temp_in, double maximum_temp_in, int spline_evaluations_in)
  : t(t_in),
	tf(tf_in),
	minimum_temp(minimum_temp_in),
	maximum_temp(maximum_temp_in),
	spline_evaluations(spline_evaluations_in) {}

void Bounce::get_splines() {
	//alglib::spline1dbuildcubic(temp_array, log_gamma_array, this->gamma_spline);

	// Use vectors to collect only valid points
	std::vector<double> valid_temps, valid_actions, valid_log_gammas;
	double dt = (maximum_temp - minimum_temp) / (spline_evaluations - 1);
	const double gamma_min = 1e-100;

	for (int i = 0; i < spline_evaluations; i++) {
		double tt = minimum_temp + i * dt;
		auto action_opt = calculate_action(tt);
		if (action_opt.has_value()) {
			double action = action_opt.value();
			auto gamma_opt = calculate_gamma(tt, action);
			double gamma = gamma_min;
			if (gamma_opt.has_value()) {
				gamma = std::max(gamma_opt.value(), gamma_min);
			}
			// std::cout << "t = " << tt << ", action = " << action << ", gamma = " << gamma << std::endl;
			valid_temps.push_back(tt);
			valid_actions.push_back(action);
			valid_log_gammas.push_back(std::log(gamma));
		} else {
			std::cout << "t = " << tt << ", action = INVALID" << std::endl;
		}
	}

	if (valid_temps.size() < 2) {
		throw std::runtime_error("Not enough valid action points to build spline.");
	}

	alglib::real_1d_array temp_array, action_array, log_gamma_array;
	temp_array.setcontent(valid_temps.size(), valid_temps.data());
	action_array.setcontent(valid_actions.size(), valid_actions.data());
	log_gamma_array.setcontent(valid_log_gammas.size(), valid_log_gammas.data());

	alglib::spline1dbuildcubic(temp_array, action_array, this->action_spline);
	alglib::spline1dbuildcubic(temp_array, action_array, this->action_spline);
	alglib::spline1dbuildcubic(temp_array, log_gamma_array, this->gamma_spline);
}

/*======================================
	Below are all definitions pertaining to
	the ThermalParameters class
	======================================*/

void ThermalParameters::find_thermal_parameters() {

	if (calculated_thermal_params) {
    return;
  }
	
	calculated_thermal_params = true;
	std::vector<Transition> trans = tf.get_transitions();

	for (auto t : trans) {
		ThermalParams tp_local;

		double minimum_temp = t.false_phase.T.front();
		double maximum_temp = t.TC;
		tp_local.TC = maximum_temp;

		if(minimum_temp > maximum_temp) { 
			throw std::runtime_error("Minimum temp exceeds maximum temp.");
		}

		if(std::abs(maximum_temp - minimum_temp) < dt_tol) {
			throw std::runtime_error("Temperature range is below set tolerance. Phase may not nucleate.");
		}

		Thermodynamics thermo_true(t.true_phase, n_temp, dof);
		Thermodynamics thermo_false(t.false_phase, n_temp, dof);
		Bounce bounce_class(t, tf, minimum_temp, maximum_temp, spline_evaluations);
		bounce_class.get_splines();

		alglib::spline1dinterpolant hubble_spline;
		this->make_hubble_spline(hubble_spline, thermo_true, thermo_false, minimum_temp, maximum_temp);
		tp_local.hubble_spline = hubble_spline;

		double tp = 0.0;
		try { 
			tp = get_percolation_temperature(hubble_spline, bounce_class, 0.35); 
			tp_local.TP = tp;
			tp_local.alpha_tp = get_alpha(tp, thermo_true, thermo_false);
			tp_local.betaH_tp = get_betaH(tp, bounce_class);
			tp_local.H_tp = get_hubble_rate(tp, thermo_true, thermo_false);
			tp_local.beta_tp = tp_local.betaH_tp * tp_local.H_tp * time_conversion;
			tp_local.percolates = true;
		} catch (const std::runtime_error &e) {
			tp_local.percolates = false;
			std::cout << "Failed to find percolation temperature: " << e.what() << std::endl;
		}
		

		double tf = 0.0;
		try { 
			tf = get_completion_temperature(hubble_spline, bounce_class, 0.35); 
			tp_local.TF = tf;
			tp_local.completes = true;
		} catch (const std::runtime_error &e) {
			tp_local.completes = false;
			std::cout << "Failed to find percolation temperature: " << e.what() << std::endl;
		}

		double tn = 0.0;
		try { 
			tn = get_nucleation_temperature(hubble_spline, bounce_class); 
			tp_local.TN = tn;
			tp_local.alpha_tn = get_alpha(tn, thermo_true, thermo_false);
			tp_local.betaH_tn = get_betaH(tn, bounce_class);
			tp_local.H_tn = get_hubble_rate(tn, thermo_true, thermo_false);
			tp_local.beta_tn = tp_local.betaH_tn * tp_local.H_tn * time_conversion;
			tp_local.nucleates = true;
		} catch (const std::runtime_error &e) {
			tp_local.nucleates = false;
			std::cout << "Failed to find percolation temperature: " << e.what() << std::endl;
		}

		thermalparameterContainer.push_back(tp_local);
	}
}

double ThermalParameters::get_alpha(double T, Thermodynamics true_thermo, Thermodynamics false_thermo) {
	return abs(true_thermo.get_theta(T) - false_thermo.get_theta(T))/false_thermo.get_enthalpy(T) * 4./3.;
}

double ThermalParameters::get_hubble_rate(double T, Thermodynamics true_thermo, Thermodynamics false_thermo) {
	double rho =  dof/30. * M_PI*M_PI * T*T*T*T + (true_thermo.get_energy(T) - false_thermo.get_energy(T));
	double h = sqrt( 8. * M_PI * G/3. * rho);
	return h;
}

double ThermalParameters::get_betaH(double T, Bounce bounce) {
	double y, dy, ddy;
	alglib::spline1ddiff(bounce.action_spline, T, y, dy, ddy);
	std::cout << "T = " << T << ", y = " << y << ", dy = " << dy << std::endl;
	return T * dy;
}

void ThermalParameters::make_hubble_spline(alglib::spline1dinterpolant& hubble_spline, Thermodynamics true_thermo, Thermodynamics false_thermo, double t_min, double t_max, double n_temp) {

	alglib::real_1d_array temp_array, integrand_array;
	temp_array.setlength(n_temp);
	integrand_array.setlength(n_temp);

	double dt = (t_max - t_min) / (n_temp - 1);
	for (int i = 0; i < n_temp; i++) {
		double tt = t_min + i * dt;
		double h = get_hubble_rate(tt, true_thermo, false_thermo);
		temp_array[i] = tt;
		integrand_array[i] = 1./h;
	}

	alglib::spline1dbuildcubic(temp_array, integrand_array, hubble_spline);
}

double ThermalParameters::hubble_integral(alglib::spline1dinterpolant& hubble_spline, double T, double Tdash) {
	double result = spline1dintegrate(hubble_spline, Tdash) - spline1dintegrate(hubble_spline, T);
	return result;
}

double ThermalParameters::false_vacuum_fraction_integrand(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double T, double Tdash) {
	double h = alglib::spline1dcalc(hubble_spline, Tdash); // 1/h !!!
	double h_int = hubble_integral(hubble_spline, T, Tdash);
	double gamma = std::exp(alglib::spline1dcalc(bounce.gamma_spline, Tdash)); // exp since we store log(gamma)
	double Tfac = 1/(Tdash*Tdash*Tdash*Tdash);

	return Tfac * gamma * h * h_int*h_int*h_int;
}

double ThermalParameters::false_vacuum_fraction(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double T, double vw) {
	double tc = bounce.maximum_temp;
	double tmin = T; // lower bound on integral
	int n_temp = 250;

	alglib::real_1d_array temp_array, integrand_array;
	alglib::spline1dinterpolant integrand_spline;
	temp_array.setlength(n_temp);
	integrand_array.setlength(n_temp);

	double dt = (tc - tmin) / (n_temp - 1);
	for (int i = 0; i < n_temp; i++) {
		double Tdash = tmin + i * dt;
		double result = false_vacuum_fraction_integrand(hubble_spline, bounce, T, Tdash);
		temp_array[i] = Tdash;
		integrand_array[i] = result;
	}

	alglib::spline1dbuildcubic(temp_array, integrand_array, integrand_spline);

	double result = - 4 * M_PI * vw*vw*vw / 3. * spline1dintegrate(integrand_spline, tc);

	return exp(result);
}

double ThermalParameters::nucleation_rate_integrand(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double T) {
	double gamma = std::exp(alglib::spline1dcalc(bounce.gamma_spline, T)); // exp since we store log(gamma)
	double h = alglib::spline1dcalc(hubble_spline, T); // 1/h !!!
	return gamma * h * h * h;
}

double ThermalParameters::nucleation_rate(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double T) {

	double tc = bounce.maximum_temp;
	double tmin = T;
	int n_temp = 250;

	alglib::real_1d_array temp_array, integrand_array;
	alglib::spline1dinterpolant integrand_spline;
	temp_array.setlength(n_temp);
	integrand_array.setlength(n_temp);

	double dt = (tc - tmin) / (n_temp - 1);
	for (int i = 0; i < n_temp; i++) {
		double tt = tmin + i * dt;
		double h = nucleation_rate_integrand(hubble_spline, bounce, tt);
		temp_array[i] = tt;
		integrand_array[i] = h;
	}

	alglib::spline1dbuildcubic(temp_array, integrand_array, integrand_spline);

	double result = 4. * M_PI / 3. * spline1dintegrate(integrand_spline, tc);
	return result;
}


double ThermalParameters::find_temperature_binary_search(double init_T, double end_T, double tol_rel, const std::function<double(double)>& calc_value, double target, int max_iterations) {
	int iterations = 0;
	double value_init = calc_value(init_T) - target;
	double value_end = calc_value(end_T) - target;
	
	while ((init_T - end_T) > tol_rel && iterations < max_iterations) {
		double mid_T = (init_T + end_T) / 2.0;
		double fmid = calc_value(mid_T) - target;

		if (fabs(fmid) < tol_rel) {
			return mid_T;
		}

		if (value_end * fmid < 0) {
			init_T = mid_T;
			value_init = fmid;
		} else {
			end_T = mid_T;
			value_end = fmid;
		}
		iterations++;
	}

	if (iterations >= max_iterations) {
		throw std::runtime_error("Failed to converge when finding temperature");
	}

	return (init_T + end_T) / 2.0;
}

double ThermalParameters::get_percolation_temperature(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double vw) {
	LOG(debug) << "Beginning percolation temperature search with vw = " << vw;

	if (vw <= 0.0 || vw >= 1.0) {
		throw std::runtime_error("Wall velocity vw must be between 0 and 1");
	}

	const double target = 0.71;
	const double tol_abs = 1e-4;
	const double tol_rel = 1e-4;
	double init_T = bounce.maximum_temp - tol_abs;
	double end_T = bounce.minimum_temp;

	if (init_T <= end_T) {
		throw std::runtime_error("Invalid temperature range for percolation search");
	}

	try {
		// Check if percolation is possible
		if (false_vacuum_fraction(hubble_spline, bounce, end_T, vw) > target + tol_abs) {
			throw std::runtime_error("Transition does not percolate");
		}

		// Create lambda for the binary search
		auto calc_value = [&](double T) {
			return false_vacuum_fraction(hubble_spline, bounce, T, vw);
		};

		return find_temperature_binary_search(init_T, end_T, tol_rel, calc_value, target);
	} catch (const alglib::ap_error& e) {
		throw std::runtime_error(std::string("ALGLIB error in percolation temperature calculation: ") + e.msg);
	}
}

double ThermalParameters::get_completion_temperature(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double vw) {
	LOG(debug) << "Beginning completion temperature search with vw = " << vw;

	if (vw <= 0.0 || vw >= 1.0) {
		throw std::runtime_error("Wall velocity vw must be between 0 and 1");
	}

	const double target = 0.01;
	const double tol_abs = 1e-4;
	const double tol_rel = 1e-4;
	double init_T = bounce.maximum_temp - tol_abs;
	double end_T = bounce.minimum_temp;

	if (init_T <= end_T) {
		throw std::runtime_error("Invalid temperature range for completion search");
	}

	try {
		// Check if completion is possible
		if (false_vacuum_fraction(hubble_spline, bounce, end_T, vw) > target + tol_abs) {
			throw std::runtime_error("Transition does not complete");
		}

		// Create lambda for the binary search
		auto calc_value = [&](double T) {
			return false_vacuum_fraction(hubble_spline, bounce, T, vw);
		};

		return find_temperature_binary_search(init_T, end_T, tol_rel, calc_value, target);
	} catch (const alglib::ap_error& e) {
		throw std::runtime_error(std::string("ALGLIB error in completion temperature calculation: ") + e.msg);
	}
}

double ThermalParameters::get_nucleation_temperature(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce) {
	LOG(debug) << "Beginning nucleation temperature search";

	const double target = 1.0;
	const double tol_rel = 1e-4;
	double init_T = bounce.maximum_temp - tol_rel;
	double end_T = bounce.minimum_temp;

	if (init_T <= end_T) {
		throw std::runtime_error("Invalid temperature range for nucleation search");
	}

	try {
		// Create lambda for the binary search
		auto calc_value = [&](double T) {
			return nucleation_rate(hubble_spline, bounce, T);
		};

		return find_temperature_binary_search(init_T, end_T, tol_rel, calc_value, target);
	} catch (const alglib::ap_error& e) {
		throw std::runtime_error(std::string("ALGLIB error in nucleation temperature calculation: ") + e.msg);
	}
}

} // namespace PhaseTracer