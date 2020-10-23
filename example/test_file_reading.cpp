#include "models/xSM_MSbar_noZ2.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "transition_graph_util.hpp"

#include <iostream>
#include <iomanip>
#include <Eigen/Eigenvalues>
#include <sstream>
#include <vector>
#include <iterator>
#include <fstream>
#include <limits>
#include <assert.h>


// Based on https://stackoverflow.com/questions/236129/how-do-i-iterate-over-the-words-of-a-string
template <typename DataType>
void split(const std::string &inputString, char delimiter, DataType result, int keepCount)
{
	std::istringstream inputStringStream(inputString);
	std::string substring;
	int i = 0;

	while(std::getline(inputStringStream, substring, delimiter) && i++ < keepCount)
	{
		*result++ = std::stod(substring);
	}
}

std::vector<double> split(const std::string &inputString, char delimiter, int keepCount)
{
	std::vector<double> elements;

	split(inputString, delimiter, std::back_inserter(elements), keepCount);

	return elements;
}

// Copied from transition_finder.cpp (called changed()), because it's a private method in TransitionFinder.
std::vector<bool> hasVEVChanged(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd&
	false_vacuum, const PhaseTracer::TransitionFinder& tf)
{
	auto delta = true_vacuum - false_vacuum;
	std::vector<bool> changed;

	for (unsigned int i = 0; i < delta.size(); ++i)
	{
		changed.push_back(std::abs(delta[i]) > tf.get_change_abs_tol()
			+ tf.get_change_rel_tol()*std::max(true_vacuum.norm(), false_vacuum.norm()));
	}

	return changed;
}

// Transitions isn't const because we modify it.
void appendSubcriticalTransitions(const std::vector<PhaseTracer::Phase>& phases, std::vector<PhaseTracer::Transition>&
	transitions, const PhaseTracer::PhaseFinder& pf, const PhaseTracer::TransitionFinder& tf)
{
	// Loop through phases and check if any phase has lower energy than any other phase when it first appears (i.e. at
	// its Tmax).
	for(int i = 0; i < phases.size()-1; ++i)
	{
		double Tmax = phases[i].T.back();
		double energyAtTmax = phases[i].V.back();

		for(int j = j+1; j < phases.size(); ++j)
		{
			PhaseTracer::Point phasejAtTmax = pf.phase_at_T(phases[j], Tmax);

			// If the phases do not overlap in temperature, there cannot be a transition between them.
			if(abs(phasejAtTmax.t - Tmax) > 1)
			{
				continue;
			}

			// If this second phase has a higher energy, then we have a subcritical transition.
			// We flag subcritical transitions with a negative critical temperature and gamma.
			if(phasejAtTmax.potential > energyAtTmax)
			{
				Eigen::VectorXd trueVacuum = phases[i].X.back();
				Eigen::VectorXd falseVacuum = phasejAtTmax.x;
				double gamma = (trueVacuum - falseVacuum).norm() / Tmax;
				std::vector<bool> changed = hasVEVChanged(trueVacuum, falseVacuum, tf);
				double deltaPotential = energyAtTmax - phasejAtTmax.potential;

				transitions.push_back({PhaseTracer::SUCCESS, -Tmax, phases[i], phases[j], trueVacuum, falseVacuum,
					-gamma, changed, deltaPotential, transitions.size()});
			}
		}
	}
}

void generateCriticalTemperatureData(std::string inputFileName, std::string outputFileName)
{
	std::ifstream inputFile(inputFileName);

	std::cout << "is_open(): " << inputFile.is_open() << std::endl;

	if(!inputFile)
	{
		std::cerr << "Cannot open the file: " << inputFileName << std::endl;
	}

	std::string line;
	std::vector<std::vector<double>> data;

	std::cout.precision(std::numeric_limits<double>::max_digits10);

	while(std::getline(inputFile, line))
	{
		if(line.size() > 0)
		{
			// Keep only the first 12 values.
			data.push_back(split(line, ' ', 12));
		}
	}

	inputFile.close();

	EffectivePotential::xSM_MSbar_noZ2 model;
	model.set_renormalization_scale(model.get_vh());
	model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
	model.set_xi(0);

	std::vector<PhaseTracer::Transition> transitions;

	std::ofstream outputFile(outputFileName);

	if(!outputFile.is_open())
	{
		std::cerr << "Could not create or open '" << outputFileName << "'" << std::endl;
		return;
	}

	// Data is correct. Now use the data to construct a potential and perform phase tracing.
	for(int i = 0; i < data.size(); ++i)
	{
		//std::cout << "===================================================================" << std::endl;
		std::cout << "Parameter point " << i << " / " << data.size()-1 << std::endl;
		//std::cout << "===================================================================" << std::endl;

		model.setParameters(data[i]);
		//model.printParameters();

		PhaseTracer::PhaseFinder pf(model);
		pf.set_t_high(300);
		pf.set_check_vacuum_at_high(false);

		PhaseTracer::TransitionFinder tf(pf);

		pf.find_phases();
		tf.find_transitions();
		//std::cout << tf;

		transitions = tf.get_transitions();

		for(int j = 0; j < transitions.size(); ++j)
		{
			outputFile << transitions[j].TC << " " << transitions[j].gamma << " " << transitions[j].false_vacuum[0]
				<< " " << transitions[j].false_vacuum[1] << " " << transitions[j].true_vacuum[0] << " "
				<< transitions[j].true_vacuum[1] << " ";
		}

		outputFile << std::endl;
	}

	outputFile.close();
}

void getCriticalTemperatureDataForPoint(std::string inputFileName, int scanPointIndex, bool bPlot)
{
	std::ifstream inputFile(inputFileName);

	std::cout << "is_open(): " << inputFile.is_open() << std::endl;

	if(!inputFile)
	{
		std::cerr << "Cannot open the file: " << inputFileName << std::endl;
	}

	std::string line;
	std::vector<double> data;

	std::cout.precision(std::numeric_limits<double>::max_digits10);

	int lineIndex = 0;

	while(std::getline(inputFile, line))
	{
		if(lineIndex++ == scanPointIndex)
		{
			if(line.size() == 0)
			{
				std::cerr << "Line " << scanPointIndex << " is empty!" << std::endl;
				return;
			}

			// Keep only the first 12 values.
			data = split(line, ' ', 12);
		}
	}

	inputFile.close();

	EffectivePotential::xSM_MSbar_noZ2 model;
	model.set_renormalization_scale(model.get_vh());
	model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
	model.set_xi(0);

	std::vector<PhaseTracer::Transition> transitions;

	model.setParameters(data);
	//model.printParameters();

	PhaseTracer::PhaseFinder pf(model);
	pf.set_t_high(300);
	pf.set_check_vacuum_at_high(false);
	// Fixes cases where there are slight discontinuities in otherwise presumably second-order transitions that are
	// treated as two distinct phases with no transition between them. This causes issues where no path is found from
	// the high temperature phase to the EW VEV due to this discontinuity.
	pf.set_check_hessian_singular(false);

	PhaseTracer::TransitionFinder tf(pf);

	pf.find_phases();
	tf.find_transitions();

	std::cout << pf;
	std::cout << tf;

	std::cout << std::endl << "============================================" << std::endl;
	std::cout << "<getTransitionPaths> Generating symmetric partners..." << std::endl;

	std::vector<PhaseTracer::Phase> symmetrisedPhases;
	std::vector<PhaseTracer::Transition> symmetrisedTransitions;

	TransitionGraph::extractExplicitSymmetricPhasesAndTransitions(pf.get_phases(), tf.get_transitions(), {0},
		symmetrisedPhases, symmetrisedTransitions);

	std::cout << "Checking if the phases at T=0 are consistent with our Universe..." << std::endl;

	Eigen::VectorXd EWVEV = model.get_EW_VEV();

	//TransitionGraph::PhaseStructureData phaseStructureData = TransitionGraph::extractPhaseStructureData(
	//	pf.get_phases(), tf.get_transitions(), EWVEV);
	TransitionGraph::PhaseStructureData phaseStructureData = TransitionGraph::extractPhaseStructureData(
		symmetrisedPhases, symmetrisedTransitions, EWVEV);

	if(phaseStructureData.validAtZeroT)
	{
		std::cout << "Phase structure is valid at T=0." << std::endl;
		std::cout << "Phase " << phaseStructureData.EWVEVIndex << " is the EW VEV phase at T=0." << std::endl;
		std::cout << "Phase " << phaseStructureData.highTPhaseIndex << " is the high temperature phase at T="
			<< pf.get_phases()[phaseStructureData.highTPhaseIndex].T.back() << "." << std::endl;

		std::cout << std::endl;
		std::cout << "Finding all transition paths from the high temperature phase to the EW VEV..." << std::endl;

		std::vector<TransitionGraph::Path> paths = TransitionGraph::getTransitionPaths(symmetrisedPhases,
			symmetrisedTransitions, phaseStructureData);
		std::cout << "Found " << paths.size() << " paths:" << std::endl;

		for(int i = 0; i < paths.size(); ++i)
		{
			bool valid = paths[i].phases.back() == phaseStructureData.EWVEVIndex;
			std::cout << (valid ? "[Valid]  " : "[Invalid]") << " Path " << i+1  << ": " << paths[i] << std::endl;
			//std::cout << paths[i] << std::endl;
		}

		std::cout << std::endl;
	}
	else
	{
		std::cout << "Phase structure is invalid at T=0. It does not describe our Universe." << std::endl;
	}

	if(bPlot)
	{
		std::string fileName = "xSM_MSbar_noZ2_point" + std::to_string(scanPointIndex);
		PhaseTracer::potential_plotter(model, tf.get_transitions().front().TC, fileName, 0., 2., 0.01, -2, 0., 0.01);
		PhaseTracer::potential_line_plotter(model, tf.get_transitions(), fileName);
		PhaseTracer::phase_plotter(tf, fileName);
	}
}

int main(int argc, char* argv[])
{
	LOGGER(fatal);

	std::string inputFileName = std::string("../../../Monash/PhD/Software/Pot2GW/Output/scans/xSM_msbar_mixing/17/")
		+ "criticalTemperatureScanPoints.txt";
	std::string outputFileName = "criticalTemperatureData2.txt";

	if(argc < 2)
	{
		std::cerr << "Missing first argument. Valid options are '-all' or '-one'." << std::endl;
		return 1;
	}

	if(strcmp(argv[1], "-all") == 0)
	{
		std::cout << "Generating critical temperature data from input file..." << std::endl;
		generateCriticalTemperatureData(inputFileName, outputFileName);
	}
	else if(strcmp(argv[1], "-one") == 0)
	{
		if(argc > 2)
		{
			int scanPointIndex = std::atoi(argv[2]);

			std::cout << "Getting critical temperature data for scan point " << scanPointIndex << " from input file..."
				<< std::endl;

			bool bPlot = argc > 3 && strcmp(argv[3], "-plot") == 0;

			getCriticalTemperatureDataForPoint(inputFileName, scanPointIndex,  bPlot);
			
			if(!bPlot)
			{
				std::cout << "(Include a third argument '-plot' if you would like the results to be plotted.)"
					<< std::endl;
			}
		}
		else
		{
			std::cerr << "No scan point index specified!" << std::endl;
		}
	}
	else
	{
		std::cerr << "Unrecognised argument '" << argv[1] << "'. Valid options are '-all' or '-one'." << std::endl;
	}
	
	/*std::string argv_str(argv[0]);
	std::string base = argv_str.substr(0, argv_str.find_last_of("/"));

	std::cout << "Base directory: " << base << std::endl;*/
	
	//std::ifstream inputFile("E:/Monash/PhD/Software/Pot2GW/Output/scans/xSM_msbar_mixing/17/data.txt");

	
}