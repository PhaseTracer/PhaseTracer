/**
  1D example program for PhaseTracer.
*/

#include <iostream>

#include "models/SubcriticalTransitionExample.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"

int main(int argc, char* argv[])
{
	const bool debug_mode = argc > 1 && strcmp(argv[1], "-d") == 0;
	
	// Set level of screen output
	if (debug_mode)
	{
		LOGGER(debug);
	}
	else
	{
		LOGGER(fatal);
	}
	
	// Construct model
	EffectivePotential::SubcriticalTransitionExample model;
	
	// Make PhaseFinder object and find the phases
	PhaseTracer::PhaseFinder pf(model);
	pf.set_t_high(2);
	pf.set_check_vacuum_at_high(false);
	pf.set_check_hessian_singular(false);
	pf.set_x_abs_identical(0.1);
	pf.set_x_abs_jump(0.05);
	pf.set_dt_max_abs(0.2);
	pf.set_phase_min_length(0.01);
	pf.find_phases();
	std::cout << pf;
	
	// Make TransitionFinder object and find the transitions
	PhaseTracer::TransitionFinder tf(pf);
	tf.find_transitions();
	std::cout << tf;
	
	if(debug_mode)
	{
		PhaseTracer::phase_plotter(tf, "1D_test_model");
	}
	
	return 0;
}
