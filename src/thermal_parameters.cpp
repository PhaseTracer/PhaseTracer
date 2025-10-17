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

	if (temperature.front() != t_min || temperature.back() != t_max) {
		throw std::runtime_error("Temperature vector does not match phase temperature bounds.");
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
		try {
			alglib::spline1ddiff(this->potential, temp, v, dvdT, ddvdT);
		} catch (const std::exception& e) {
			LOG(error) << "Error in spline1ddiff: " << e.what() << "for temperature " << temp << std::endl;
			continue;
		}

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

	try {
		alglib::spline1dbuildcubic(t_array, p_array, this->p);
		alglib::spline1dbuildcubic(t_array, e_array, this->e);
		alglib::spline1dbuildcubic(t_array, w_array, this->w);
		alglib::spline1dbuildcubic(t_array, s_array, this->s);
	} catch (const std::exception& e) {
		LOG(error) << "Error in spline1dbuildcubic: " << e.what() << std::endl;
		throw std::runtime_error("Failed to build thermodynamic splines");
	}
}

std::vector<double> Thermodynamics::get_potential(double T) {
  check_temperature_range(T, "get_potential");
  double v, dvdT, ddvdT;
  alglib::spline1ddiff(this->potential, T, v, dvdT, ddvdT);
  return {v, dvdT, ddvdT};
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

std::vector<double> Thermodynamics::get_pressure_derivs(double T) const {
  check_temperature_range(T, "get_pressure_derivs");

  double pp, dp, ddp;
  alglib::spline1ddiff(this->p, T, pp, dp, ddp);

  return {pp, dp, ddp};
}

std::vector<double> Thermodynamics::get_energy_derivs(double T) const {
  check_temperature_range(T, "get_energy_derivs");

  double ee, de, dde;
  alglib::spline1ddiff(this->e, T, ee, de, dde);

  return {ee, de, dde};
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

	if (spline_evaluations < 2) {
		throw std::runtime_error("Spline evaluations must be at least 2.");
	}

	std::vector<double> valid_temps, valid_actions, valid_log_gammas;
	double dt = (maximum_temp - minimum_temp) / (spline_evaluations - 1);
	const double log_gamma_min = -700;

	for (int i = 0; i < spline_evaluations; i++) {
		double tt = minimum_temp + i * dt;
		auto action_opt = calculate_action(tt);
		if (action_opt.has_value()) {
			double action = action_opt.value();
	auto log_gamma_opt = calculate_log_gamma(tt, action);
	double log_gamma = log_gamma_min;
	if (log_gamma_opt.has_value()) {
		log_gamma = std::max(log_gamma_opt.value(), log_gamma_min);
	}
	valid_temps.push_back(tt);
	valid_actions.push_back(action);
	valid_log_gammas.push_back(log_gamma);
		} else {
			LOG(debug) << "t = " << tt << ", action = INVALID" << std::endl;
		}
	}

	if (valid_temps.size() < 2) {
		throw std::runtime_error("Not enough valid action points to build spline.");
	}

	alglib::real_1d_array temp_array, action_array, log_gamma_array;
	temp_array.setcontent(valid_temps.size(), valid_temps.data());
	action_array.setcontent(valid_actions.size(), valid_actions.data());
	log_gamma_array.setcontent(valid_log_gammas.size(), valid_log_gammas.data());

	LOG(debug) << "Building action and gamma arrays";
	alglib::spline1dbuildcubic(temp_array, action_array, this->action_spline);
	alglib::spline1dbuildcubic(temp_array, log_gamma_array, this->log_gamma_spline);
	LOG(debug) << "Action and log(gamma) arrays built.";
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
		try 
		{
			ThermalParams tp_local;
			try 
			{
				tp_local = process_transition(t);
			} catch(...) 
			{
				tp_local.success = false;
			}

			if (tp_local.success);
			{
				thermalparameterContainer.push_back(tp_local);
			}
		} catch (...) {
			continue;
		}
	}
}

ThermalParams ThermalParameters::process_transition(Transition t) {

	ThermalParams output;

	auto [minimum_temp, maximum_temp] = setup_temperature_range(t, output);
	if(output.success = false) {return output;}

	auto [thermo_true, thermo_false] = setup_thermodynamics_classes(t);

	auto bounce_class = setup_bounce_class(t, tf, minimum_temp, maximum_temp, spline_evaluations);

	alglib::spline1dinterpolant hubble_spline;
	this->make_hubble_spline(hubble_spline, thermo_true, thermo_false, minimum_temp, maximum_temp);
	output.hubble_spline = hubble_spline;
	LOG(debug) << "Created Hubble spline.";

	alglib::spline1dinterpolant dtdT_spline;
	this->make_dtdT_spline(dtdT_spline, thermo_true, thermo_false, minimum_temp, maximum_temp);
	output.dtdT_spline = dtdT_spline;
	LOG(debug) << "Created dtdT spline.";

	if(compute_debug) {
		add_debug_information(output, minimum_temp, maximum_temp, hubble_spline, bounce_class, thermo_true, thermo_false, dtdT_spline);
	}

	output.eos = get_eos(minimum_temp, maximum_temp, thermo_true, thermo_false);
	
	calculate_completion_milestone(output, t, hubble_spline, dtdT_spline, bounce_class, thermo_true, thermo_false);
	calculate_percolation_milestone(output, t, hubble_spline, dtdT_spline, bounce_class, thermo_true, thermo_false);
	calculate_nucleation_milestone(output, t, hubble_spline, dtdT_spline, bounce_class, thermo_true, thermo_false);

	calculate_durations(output);

	return output;
}

std::pair<double, double>
ThermalParameters::setup_temperature_range(const Transition& t, ThermalParams& output)
{
	LOG(debug) << "Setting up temperature ranges for transition " << t.key;
	double minimum_temp = t.false_phase.T.front();
	double maximum_temp = t.TC;
	output.TC = maximum_temp;

	if(minimum_temp > maximum_temp) { 
		std::swap(minimum_temp, maximum_temp);
	}

	if(std::abs(maximum_temp - minimum_temp) < dt_tol) {
		output.success = false;
	}

	return {minimum_temp, maximum_temp};
}

std::pair<Thermodynamics, Thermodynamics>
ThermalParameters::setup_thermodynamics_classes(const Transition& t)
{
	LOG(debug) << "Creating true phase Thermodynamics class for transition " << t.key;
	Thermodynamics thermo_true(t.true_phase, n_temp, dof);
	LOG(debug) << "Creating false phase Thermodynamics class for transition " << t.key;
	Thermodynamics thermo_false(t.false_phase, n_temp, dof);

	return {thermo_true, thermo_false};
}

Bounce
ThermalParameters::setup_bounce_class(
    const Transition&t, 
    TransitionFinder tf, 
    const double& minimum_temp, 
    const double& maximum_temp, 
    const int& spline_evaluations
)
{
	LOG(debug) << "Creating bounce class for transition " << t.key;
	Bounce bounce_class(t, tf, minimum_temp, maximum_temp, spline_evaluations);
	LOG(debug) << "Creating bounce spline for transition " << t.key;
	bounce_class.get_splines();
	LOG(debug) << "Bounce.get_splines() successful for transition " << t.key;
	return bounce_class;
}

void 
ThermalParameters::add_debug_information(
	ThermalParams& output,
	const double& minimum_temp,
	const double& maximum_temp,
	alglib::spline1dinterpolant hubble_spline,
	Bounce bounce_class,
	Thermodynamics& thermo_true,
	Thermodynamics& thermo_false,
	alglib::spline1dinterpolant dtdT_spline
)
{
	LOG(debug) << "Creating debug info for thermal parameters";
	ThermalDEBUG debug;
	int n = 75;
	double dt = (maximum_temp - minimum_temp)/(n-1.);
	for (double tt = minimum_temp; tt < maximum_temp; tt += dt) {
		debug.temp.push_back(tt);
		double pf;
		try {pf = false_vacuum_fraction(hubble_spline, bounce_class, tt, vw);} catch (...) {pf = 1.0;}
		debug.pf.push_back(pf);
		double nt;
		try {nt = nucleation_rate(hubble_spline, bounce_class, tt);} catch (...) {nt = 0.0;}
		debug.nt.push_back(nt);
		double gamma;
		try {gamma = std::exp(alglib::spline1dcalc(bounce_class.log_gamma_spline, tt));} catch (...) {gamma = 0.0;}
		debug.gamma.push_back(gamma);
		double hubble_rate = get_hubble_rate(tt, thermo_true, thermo_false);
		debug.hubble.push_back(hubble_rate);
		double t;
		try {t = get_duration(dtdT_spline, maximum_temp, tt);} catch (...) {t = 0.0;}
		debug.t.push_back(t);
		LOG(debug) << "Debug info for T = " << tt << ": pf = " << pf << ", nt = " << nt << ", gamma = " << gamma << ", hubble_rate = " << hubble_rate << ", t = " << t;
	}
	debug.tmin = minimum_temp;
	debug.tmax = maximum_temp;
	output.debug_info = debug;
	LOG(debug) << "Debug info created.";
}

EoS ThermalParameters::get_eos(double Tmin, double Tmax, Thermodynamics true_thermo, Thermodynamics false_thermo) {

	EoS out;
	int N = 300;
	if(N % 2 == 1) { N += 1; }
	double dt = (Tmax-Tmin)/(N-1);
	for (double tt = Tmin; tt < Tmax + 1e-6; tt += dt) 
	{
		out.temp.push_back(tt);
		out.pressure_true.push_back(true_thermo.get_pressure(tt));
		out.energy_true.push_back(true_thermo.get_energy(tt));
		out.enthalpy_true.push_back(true_thermo.get_enthalpy(tt));
		out.entropy_true.push_back(true_thermo.get_entropy(tt));
		out.pressure_false.push_back(false_thermo.get_pressure(tt));
		out.energy_false.push_back(false_thermo.get_energy(tt));
		out.enthalpy_false.push_back(false_thermo.get_enthalpy(tt));
		out.entropy_false.push_back(false_thermo.get_entropy(tt));
	}
	return out;
}

ThermalParameters::ThermalParamSet
ThermalParameters::calculate_thermal_parameters(
	const double& tref,
	Thermodynamics& thermo_true, 
	Thermodynamics& thermo_false,
	Bounce& bounce_class
)
{	
	const auto alpha = get_alpha(tref, thermo_true, thermo_false);
	const auto betaH = get_betaH(tref, bounce_class);
	const auto H = get_hubble_rate(tref, thermo_true, thermo_false);
	const auto we = get_we(tref, thermo_true);
	const auto cs_true = get_sound_speed(tref, thermo_true);
	const auto cs_false =get_sound_speed(tref, thermo_false);
	ThermalParamSet output(alpha, betaH, H, we, cs_true, cs_false);
	return output;
}

void 
ThermalParameters::calculate_completion_milestone(
	ThermalParams& output,
	Transition& t,
    alglib::spline1dinterpolant& hubble_spline, 
    alglib::spline1dinterpolant& dtdT_spline,
    Bounce& bounce_class,
    Thermodynamics& thermo_true,
    Thermodynamics& thermo_false
)
{
	try { 
		double tfin = get_completion_temperature(hubble_spline, bounce_class, vw);
		LOG(debug) << "Found completion temp " << tfin;
		output.TF = tfin;

		if(abs(tfin-t.TC) <= 1e-3) {
			LOG(debug) << "Completion temperature is close to critical temperature.";
			output.completes = MilestoneStatus::FAST;
			return;
		}
		
		output.completes = MilestoneStatus::YES;
	} catch (const std::runtime_error &e) {
		output.completes = MilestoneStatus::NO;
		LOG(debug) << "Failed to find completion temperature: " << e.what() << std::endl;
	}
}

void 
ThermalParameters::calculate_percolation_milestone(
	ThermalParams& output,
	Transition& t,
    alglib::spline1dinterpolant& hubble_spline, 
    alglib::spline1dinterpolant& dtdT_spline, 
    Bounce& bounce_class,
    Thermodynamics& thermo_true,
    Thermodynamics& thermo_false
)
{
	try { 
		double tp = get_percolation_temperature(hubble_spline, bounce_class, vw);
		LOG(debug) << "Found percolation temp " << tp;
		output.TP = tp;

		if(abs(tp-t.TC) <= 1e-3) {
			LOG(debug) << "Percolation temperature is close to critical temperature.";
			output.percolates = MilestoneStatus::FAST;
			return;
		}
		
		const auto parameter_set = calculate_thermal_parameters(tp, thermo_true, thermo_false, bounce_class);
		output.alpha_tp = parameter_set.alpha;
		output.betaH_tp = parameter_set.betaH;
		output.H_tp = parameter_set.H;
		output.we_tp = parameter_set.we;
		output.cs_true_tp = parameter_set.cs_true;
		output.cs_false_tp = parameter_set.cs_false;

		output.percolates = MilestoneStatus::YES;
	} catch (const std::runtime_error &e) {
		output.percolates = MilestoneStatus::NO;
		LOG(debug) << "Failed to find percolation temperature: " << e.what() << std::endl;
	}
}

void 
ThermalParameters::calculate_nucleation_milestone(
	ThermalParams& output,
	Transition& t,
    alglib::spline1dinterpolant& hubble_spline, 
    alglib::spline1dinterpolant& dtdT_spline, 
    Bounce& bounce_class,
    Thermodynamics& thermo_true,
    Thermodynamics& thermo_false
)
{
	try { 
		double tn = get_nucleation_temperature(hubble_spline, bounce_class); 
		LOG(debug) << "Found nucleation temp " << tn;
		output.TN = tn;

		if(abs(tn-t.TC) <= 1e-3) {
			output.nucleates = MilestoneStatus::FAST;
			LOG(debug) << "Nucleation temperature is close to critical temperature.";
			return;
		}
	
		const auto parameter_set = calculate_thermal_parameters(tn, thermo_true, thermo_false, bounce_class);
		output.alpha_tn = parameter_set.alpha;
		output.betaH_tn = parameter_set.betaH;
		output.H_tn = parameter_set.H;
		output.we_tn = parameter_set.we;
		output.cs_true_tn = parameter_set.cs_true;
		output.cs_false_tn = parameter_set.cs_false;

		if ( output.completes == MilestoneStatus::YES && output.TF > tn ) 
		{
			output.nucleates = MilestoneStatus::INVALID;
		} else 
		{
			output.nucleates = MilestoneStatus::YES;
		}
	} catch (const std::runtime_error &e) {
		output.nucleates = MilestoneStatus::NO;
		LOG(debug) << "Failed to find nucleation temperature: " << e.what() << std::endl;
	}
}

void
ThermalParameters::calculate_durations(ThermalParams& output)
{
	if (output.completes == MilestoneStatus::YES)
	{
		output.dtf_tc = get_duration(output.dtdT_spline, output.TC, output.TF) * output.H_tp;
		if (output.percolates == MilestoneStatus::YES)
		{
			output.dtf_tp = get_duration(output.dtdT_spline, output.TP, output.TF) * output.H_tp;
		}
		if (output.nucleates == MilestoneStatus::YES)
		{
			output.dtf_tn = get_duration(output.dtdT_spline, output.TN, output.TF) * output.H_tn;
		}
	}
	
	
}

double ThermalParameters::get_sound_speed(double T, Thermodynamics phase_thermo) {
	double dedT = phase_thermo.get_energy_derivs(T)[1];
	double dPdT = phase_thermo.get_pressure_derivs(T)[1];

	return sqrt(dPdT/dedT);
}

double ThermalParameters::get_we(double Tref, Thermodynamics true_thermo) {
	return true_thermo.get_enthalpy(Tref)/true_thermo.get_energy(Tref);
}

double ThermalParameters::get_duration(alglib::spline1dinterpolant& dtdT_spline, double tStart, double tEnd) {

	double duration = alglib::spline1dintegrate(dtdT_spline, tStart) - alglib::spline1dintegrate(dtdT_spline, tEnd);
	return duration;
}

double ThermalParameters::get_alpha(double T, Thermodynamics true_thermo, Thermodynamics false_thermo) {
	return abs(true_thermo.get_theta(T) - false_thermo.get_theta(T))/false_thermo.get_enthalpy(T) * 4./3.;
}

double ThermalParameters::get_hubble_rate(double T, Thermodynamics true_thermo, Thermodynamics false_thermo) {
	double rho = dof/30. * M_PI*M_PI * T*T*T*T + abs(true_thermo.get_energy(T) - false_thermo.get_energy(T));
	if(isnan(rho) || rho < 0.0) { rho = dof/30. * M_PI*M_PI * T*T*T*T; } // fall back to RD result
	double h = sqrt( 8. * M_PI * G/3. * rho);
	return h;
}

double ThermalParameters::get_betaH(double T, Bounce bounce) {
	double y, dy, ddy;
	alglib::spline1ddiff(bounce.action_spline, T, y, dy, ddy);
	LOG(debug) << "For betaH calc: T = " << T << ", y = " << y << ", dy = " << dy << std::endl;
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

void ThermalParameters::make_dtdT_spline(alglib::spline1dinterpolant& dtdT_spline, Thermodynamics true_thermo, Thermodynamics false_thermo, double t_min, double t_max, double n_temp) {

	alglib::real_1d_array temp_array, dtdT_array;
	temp_array.setlength(n_temp);
	dtdT_array.setlength(n_temp);

	double dt = (t_max - t_min) / (n_temp - 1);
	for (int i = 0; i < n_temp; i++) {
		double tt = t_min + i * dt;
		std::vector<double> pressure_vals = true_thermo.get_pressure_derivs(tt);
		double dp = pressure_vals[1];
		double ddp = pressure_vals[2];
		double dTdt = 3. * get_hubble_rate(tt, true_thermo, false_thermo) * dp/ddp;
		temp_array[i] = tt;
		dtdT_array[i] = 1/dTdt;
	}

	alglib::spline1dbuildcubic(temp_array, dtdT_array, dtdT_spline);
}

double ThermalParameters::hubble_integral(alglib::spline1dinterpolant& hubble_spline, double T, double Tdash) {
	double result = spline1dintegrate(hubble_spline, Tdash) - spline1dintegrate(hubble_spline, T);
	return result;
}

double ThermalParameters::false_vacuum_fraction_integrand(alglib::spline1dinterpolant& hubble_spline, const Bounce& bounce, double T, double Tdash) {
  double h = alglib::spline1dcalc(hubble_spline, Tdash); // 1/h !!!
  double h_int = hubble_integral(hubble_spline, T, Tdash);
  double log_gamma = alglib::spline1dcalc(bounce.log_gamma_spline, Tdash); // log(gamma)
  double gamma = std::exp(log_gamma);
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
  double log_gamma = alglib::spline1dcalc(bounce.log_gamma_spline, T); // log(gamma)
  double gamma = std::exp(log_gamma);
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