#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <sys/stat.h>
#include <filesystem>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <array>
#include <memory>
#include <unistd.h>

#include "dataanalysis.h"
#include <interpolation.h>

#include "DRalgo_IDM.hpp"
#include "phasetracer.hpp"
#include "helperIncludes/gwTools.hpp"
#include "helperIncludes/wallGoWrapper.hpp"
#include "deep_phase.hpp"

using json = nlohmann::json;

/*
  Reads .json input parameters
*/
json readFile(std::string fileName){
  std::ifstream file;
  file.open(fileName);

  if(file.fail()){
    throw std::runtime_error("error loading model parameters file!");
  }

  json data = json::parse(file);

  file.close();

  return data;
}

/*
  Checks the output file exists and is non-empty
*/
bool file_exists_and_not_empty(const std::string& filename) {
  struct stat buffer;
  if (stat(filename.c_str(), &buffer) != 0) return false;
  return buffer.st_size > 0;
}

std::string
format_thermal_params(std::vector<double> in, int flag, double TC, PhaseTracer::TransitionMilestone milestone)
{
  std::string output = "";

  for ( double val : in ) { output += std::to_string(val) + ","; }

  output += std::to_string(flag) + ",";
  double Tc = TC;
  output += std::to_string(Tc) + ",";

  if(milestone.status == PhaseTracer::MilestoneStatus::YES) {
    output += std::to_string(milestone.temperature) + ","
    + std::to_string(milestone.alpha) + ","
    + std::to_string(milestone.betaH) + ","
    + std::to_string(milestone.Rs) + ","
    + std::to_string(milestone.dt) + ","
    + std::to_string(milestone.we) + ","
    + std::to_string(milestone.cs_plus) + ","
    + std::to_string(milestone.cs_minus) + ",";
  } else {
    output += "0,0,0,0,0,0,0,";
  }

  if ( !output.empty() ) { output.pop_back(); }
  output += "\n";

  return output;
}

std::string
format_fail_case(std::vector<double> in, int flag)
{
  std::string output = "";

  for ( double val : in ) { output += std::to_string(val) + ","; }

  output += std::to_string(flag) + ",";

  output += "0,0,0,0,0,0,0,0,0\n";

  return output;
}

/*
  Main
*/
int main(int argc, char* argv[]) {

	// default output
	std::string output_filename = "comparisonData/nucleationType/scan/idm_lambdaL.csv";
	bool file_has_data = file_exists_and_not_empty(output_filename);

	std::ofstream output_file(output_filename, std::ios::app);
	if (!file_has_data) {
		output_file << "mh,mH,mA,mHpm,lambda2,lambdaL,";
		output_file << "success,Tc,";
		output_file << "T,alpha,betaH,RsH,dtH,we,cs_plus,cs_minus\n";
	}

	double m12, mh, mH, mA, mHpm, lambda2, lambdaL;
  	double tHigh, tLow;
	double vw;
	double dtauRs;

	LOGGER(debug);

	std::string json_filename;
	if ( argc == 1 ) 
	{
		json_filename = "example/TestDeepPhase/modelParams_IDM.json";
	} else 
	{
		json_filename = argv[1];
		std::cout << json_filename << std::endl;
	}
	json modelParams = readFile(json_filename);

	std::cout << "Read model params\n";

	try {
		mh = modelParams["mh"].get<double>();
		mH = modelParams["mH"].get<double>();
		mA = modelParams["mA"].get<double>();
		mHpm = modelParams["mHpm"].get<double>();
		lambda2 = modelParams["lambda2"].get<double>();
		lambdaL = modelParams["lambdaL"].get<double>();
		tHigh = modelParams["High T"].get<double>();
		tLow = modelParams["Low T"].get<double>();
		vw = 1/sqrt(3.);
		dtauRs = 10.0;
	} catch (...) {
		std::cerr << "Input file might not be formatted correctly! Quitting now." << std::endl;
		output_file << format_fail_case({-1.,-1.,-1.,-1.,-1.,-1.}, -10);
    	output_file.close();
		return 1;
	}

	std::vector<double> in = {mh, mH, mA, mHpm, lambda2, lambdaL};

	std::unique_ptr<EffectivePotential::DR_idm> model;
	try {
		model = std::make_unique<EffectivePotential::DR_idm>(in);
	} catch (const std::exception& e) {
		std::cerr << "Error constructing DR_idm model: " << e.what() << std::endl;
		output_file << format_fail_case(in, -1);
    	output_file.close();
		return 1;
	} catch (...) {
		std::cerr << "Unknown error during DR_idm construction." << std::endl;
		output_file << format_fail_case(in, -1);
    	output_file.close();
		return 1;
	}

	PhaseTracer::PhaseFinder pf(*model);
	pf.set_seed(0);
	pf.set_check_hessian_singular(false);
	pf.set_check_vacuum_at_high(false);
	pf.set_check_vacuum_at_low(false);
	pf.set_t_high(tHigh);
	pf.set_t_low(tLow);

	try {
	  	pf.find_phases();
	} catch (const std::exception& e) {
		std::cerr << "Exception in pf.find_phases(): " << e.what() << std::endl;
		output_file << format_fail_case(in, -2);
		output_file.close();
	  return 1;
	} catch (...) {
		std::cerr << "Unknown exception in pf.find_phases()" << std::endl;
		output_file << format_fail_case(in, -2);
		output_file.close();
		return 1;
	}
	auto p1 = pf.get_phases();

	std::cout << pf;

	// check for phases
	if (p1.size() == 0) {
	  	return 1;
	}

	PhaseTracer::ActionCalculator ac(pf);
	PhaseTracer::TransitionFinder tf(pf);
	try {
	  	tf.find_transitions();
	} catch (const std::exception& e) {
		std::cerr << "Exception in tf.find_transitions(): " << e.what() << std::endl;
		output_file << format_fail_case(in, -3);
		output_file.close();
		return 1;
	} catch (...) {
		std::cerr << "Unknown exception in tf.find_transitions()" << std::endl;
		output_file << format_fail_case(in, -4);
		output_file.close();
		return 1;
	}
	std::cout << tf;

	auto t = tf.get_transitions();
	// check for transitions
	if (t.size() == 0) {
		output_file << format_fail_case(in, -5);
		output_file.close();
		return 1;
	}

	PhaseTracer::ThermoFinder tm(ac);
	tm.set_vw(vw);
	tm.set_dof(110.75); 
	tm.set_background_dof(110.75);
	tm.set_percolation_print_setting(PhaseTracer::PrintSettings::VERBOSE);
	tm.set_nucleation_print_setting(PhaseTracer::PrintSettings::VERBOSE);
  	tm.set_compute_profiles(true);

	PhaseTracer::TransitionMilestone milestone;
  	std::optional<PhaseTracer::EquationOfState> eos_opt;

	int percolates, nuc_type;
	try {
		auto tps = tm.get_thermal_parameter_set(t[0]);

		milestone = tps.percolation;
		eos_opt = std::move(tps.eos);
		
		if (milestone.status == PhaseTracer::MilestoneStatus::YES)
		{
			percolates = 1;
			nuc_type = (milestone.nucleation_type == PhaseTracer::NucleationType::EXPONENTIAL) ? 1 : 0;
		} else {
			percolates = 0;
			nuc_type = 2;
		}
		std::string output_file_path = "comparisonData/nucleationType/idm/thermal_profiles/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + ".csv";
		tps.profiles.write(output_file_path);
	} catch(...) {
		std::cout << "Error during calculation of thermal parameters." << std::endl;
		output_file << format_fail_case(in, -6);
		output_file.close();
		return 1;
	}
  
	if(milestone.status != PhaseTracer::MilestoneStatus::YES) {
		std::cout << "Milestone not achieved." << std::endl;
		output_file << format_fail_case(in, -7);
		output_file.close();
		return 1;
	}

	std::cout << "============================" << "\n";
  	std::cout << milestone;
	std::cout << "============================" << "\n";

	// creata data for effective potential at T = Tref
	if (percolates == 1) {
		const double tt = milestone.temperature;
		const double tc = t[1].TC;
		const double phi_min = -300.0;
		const double phi_max = 300.0;
		const int Npoints = 1000;
		const std::string path = "comparisonData/nucleationType/idm/potential/" + std::to_string(lambdaL) + ".csv";
		std::ofstream pot_file(path);
		for ( int i = 0; i <= Npoints; ++i ) {
			double phi = phi_min + i * (phi_max - phi_min) / Npoints;
			Eigen::VectorXd field_vec(1);
			field_vec(0) = phi;
			double V_tref = model->V(field_vec, tt);
			double V_tc = model->V(field_vec, tc);
			pot_file << phi << "," << V_tref << "," << V_tc << "\n";
		}
		pot_file.close();
	}

	output_file << format_thermal_params(in, 1, t[0].TC, milestone);
  	output_file.close();

	// return 0;

	std::string eos_path = "comparisonData/nucleationType/idm/eos/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + ".csv";
  	eos_opt->write(eos_path);

  	namespace DPInt = PhaseTracer::DeepPhaseInterface;

	if(vw > 0) {
		try {
			auto deflagaration_bag = DPInt::DeepPhaseResults(milestone, *eos_opt, 110.75, vw, dtauRs, DPInt::EoSModel::BAG);
			auto def_bag_profile = deflagaration_bag.get_profile_bag();
			std::string def_bag_path = "comparisonData/nucleationType/idm/fluid_profiles/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + "_bag.csv";
			def_bag_profile.write(def_bag_path);
			auto def_bag_spectrum = deflagaration_bag.get_spectrum_bag();
			std::string def_bag_spectrum_path = "comparisonData/nucleationType/idm/spectra/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + "_bag.csv";
			def_bag_spectrum.write(def_bag_spectrum_path);
		} catch (const std::exception& e) {
			std::cerr << "BAG model failed for lambdaL=" << lambdaL << ": " << e.what() << std::endl;
		} catch (...) {
			std::cerr << "BAG model failed for lambdaL=" << lambdaL << ": unknown error" << std::endl;
		}

		try {
			auto deflagaration_munu = DPInt::DeepPhaseResults(milestone, *eos_opt, 110.75, vw, dtauRs, DPInt::EoSModel::MUNU);
			auto def_munu_profile = deflagaration_munu.get_profile_munu();
			std::string def_munu_path = "comparisonData/nucleationType/idm/fluid_profiles/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + "_munu.csv";
			def_munu_profile.write(def_munu_path);
			auto def_munu_spectrum = deflagaration_munu.get_spectrum_munu();
			std::string def_munu_spectrum_path = "comparisonData/nucleationType/idm/spectra/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + "_munu.csv";
			def_munu_spectrum.write(def_munu_spectrum_path);
		} catch (const std::exception& e) {
			std::cerr << "MUNU model failed for lambdaL=" << lambdaL << ": " << e.what() << std::endl;
		} catch (...) {
			std::cerr << "MUNU model failed for lambdaL=" << lambdaL << ": unknown error" << std::endl;
		}

		try {
			auto deflagaration_veff = DPInt::DeepPhaseResults(milestone, *eos_opt, 106.75 + 8.00, vw, dtauRs, DPInt::EoSModel::VEFF);
			auto def_veff_profile = deflagaration_veff.get_profile_veff();
			std::string def_veff_path = "comparisonData/nucleationType/idm/fluid_profiles/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + "_veff.csv";
			def_veff_profile.write(def_veff_path);
			auto def_veff_spectrum = deflagaration_veff.get_spectrum_veff();
			std::string def_veff_spectrum_path = "comparisonData/nucleationType/idm/spectra/" + std::to_string(lambdaL) + "_" + std::to_string(percolates) + "_" + std::to_string(nuc_type) + "_veff.csv";
			def_veff_spectrum.write(def_veff_spectrum_path);
		} catch (const std::exception& e) {
			std::cerr << "VEFF model failed for lambdaL=" << lambdaL << ": " << e.what() << std::endl;
		} catch (...) {
			std::cerr << "VEFF model failed for lambdaL=" << lambdaL << ": unknown error" << std::endl;
		}
	}


	return 0;
}