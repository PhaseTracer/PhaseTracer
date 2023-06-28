/**
	Generates and outputs the phase structure for a toy model.
*/

#include <iostream>
#include <string.h>

#include "models/ToyModel.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"
#include "transition_graph_util.hpp"

#include <Eigen/Eigenvalues>

void printPaths(std::vector<TransitionGraph::Path> paths)
{
	if (paths.size() == 0)
	{
		std::cout << "Found no paths!" << std::endl << std::endl;
		return;
	}

	std::cout << "Found " << paths.size() << " paths:" << std::endl;

	for (int i = 0; i < paths.size(); ++i)
	{
		std::cout << "Path " << i + 1 << ": " << paths[i] << std::endl;
	}

	std::cout << std::endl;
}

int main(int argc, char* argv[])
{
	bool check_subcrit = false;
	bool allow_phase_oscillation = false;
	bool bDebug = false;
	bool bTrace = false;
	bool bPlot = false;
	bool bNoTransitionPathFinding = false;
	double dx = -1;
	double dt = -1;

	std::vector<std::string> args;

	if(argc > 1)
	{
		// See https://stackoverflow.com/questions/15344714/convert-command-line-argument-to-string
		// +1 for the starting pointer so we skip over the executable name argv[0].
		args.assign(argv, argv + argc);
	}

	std::string inputFileName = args[1];
	std::string outputFolderName = args[2];

	// Check for additional input configuration settings.
	for(int i = 3; i < argc; ++i)
	{
		//std::cout << i << ": " << argv[i] << std::endl;

		if(!bDebug && !bTrace && strcmp(argv[i], "-debug") == 0)
		{
			LOGGER(debug);
			bDebug = true;
			continue;
		}

		if(!bTrace && strcmp(argv[i], "-trace") == 0)
		{
			LOGGER(trace);
			bTrace = true;
			continue;
		}

		if(!bPlot && strcmp(argv[i], "-plot") == 0)
		{
			bPlot = true;
			continue;
		}

		if(!check_subcrit && strcmp(argv[i], "-subcrit") == 0)
		{
			check_subcrit = true;
			continue;
		}

		if(!allow_phase_oscillation && strcmp(argv[i], "-osc") == 0)
		{
			allow_phase_oscillation = true;
			continue;
		}

		if(!bNoTransitionPathFinding && strcmp(argv[i], "-no_tpf") == 0)
		{
			bNoTransitionPathFinding = true;
			continue;
		}

		if(dx < 0 && args[i].compare(0, 3, "dx=") == 0)
		{
			dx = std::stod(args[i].substr(3, args[i].size()-3));
			continue;
		}

		if(dt < 0 && args[i].compare(0, 3, "dt=") == 0)
		{
			dt = std::stod(args[i].substr(3, args[i].size()-3));
			continue;
		}
	}
	
	// Set level of screen output
	if(!bDebug && !bTrace)
	{
		LOGGER(fatal);
	}
	
	// Construct model
	EffectivePotential::ToyModel model(inputFileName);
	
	// Make PhaseFinder object and find the phases
	PhaseTracer::PhaseFinder pf(model);
	pf.set_t_high(140);
	pf.set_check_vacuum_at_high(false);
	pf.set_check_hessian_singular(false);
	pf.find_phases();
	
	// Make TransitionFinder object and find the transitions
	PhaseTracer::TransitionFinder tf(pf);
	tf.set_check_subcritical_transitions(check_subcrit);
	tf.set_assume_only_one_transition(!allow_phase_oscillation);
	tf.find_transitions();
	tf.find_transition_paths(model, false);

	//std::cout << pf;
	//std::cout << tf;

	if(!bNoTransitionPathFinding)
	{
		LOG(debug) << pf;
		LOG(debug) << "Finding transition paths...";
		tf.find_transition_paths(model, true);
	}

	if(bDebug)
	{
		std::cout << pf;
		std::cout << tf;

		if (!bNoTransitionPathFinding)
		{
			printPaths(tf.get_transition_paths());
		}
	}

	PhaseTracer::phase_plotter(tf, outputFolderName, "phase_structure", bPlot);
	
	return 0;
}
