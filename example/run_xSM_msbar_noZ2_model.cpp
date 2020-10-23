#include "models/xSM_MSbar_noZ2.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "logger.hpp"

#include <iostream>
#include <iomanip>
#include <Eigen/Eigenvalues>


int main(int argc, char* argv[])
{
	const bool debug_mode = argc > 1 && strcmp(argv[1], "-d") == 0;

	// Set level of screen output.
	if(debug_mode)
	{
		LOGGER(debug);
	}
	else
	{
		LOGGER(fatal);
	}

	// Construct our model.
	/*const double muh_sq = 6.165772569552732166e+05;
	const double mus_sq = -3.919985889020391041e+05;
	const double lh = 1.014673237931146416e-01;
	const double ls = -5.712709273123119796e-02;
	const double c1 = -1.685546619011946404e+03;
	const double c2 = 1.075016951613527949e+00;
	const double b3 = 3.434443203487260234e+02;
	const double muh_sq_tree_EWSB = 6.260801362621320877e+05;
	const double mus_sq_tree_EWSB = -3.898159562881914899e+05;
	const double vs = 9.280544295376481614e+02;
	const double ms = 8.546474096471117718e+02;
	const double theta = -2.565024655979150081e-01;*/
	
	const double muh_sq = 4.395235615305439569e+05;
	const double mus_sq = 1.270896332544635225e+05;
	const double lh = 6.243549301412576469e-02;
	const double ls = 2.560973572166680978e-01;
	const double c1 = -1.564255343078616534e+03;
	const double c2 = 1.349490036530488934e+00;
	const double b3 = -3.580910027023583098e+02;
	const double muh_sq_tree_EWSB = 4.420772632609364809e+05;
	const double mus_sq_tree_EWSB = 1.281930288220542279e+05;
	const double vs = 6.315965328619600996e+02;
	const double ms = 5.871141051491841836e+02;
	const double theta = -2.814157977928305465e-01;

	EffectivePotential::xSM_MSbar_noZ2 model(muh_sq, mus_sq, lh, ls, c1, c2, b3, muh_sq_tree_EWSB, mus_sq_tree_EWSB,
		vs, ms, theta);

	model.set_renormalization_scale(model.get_vh());
	model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
	model.set_xi(0);
	
	Eigen::VectorXd origin(2);
	origin << 0., 0.;
	Eigen::VectorXd vev = model.get_EW_VEV();

	std::cout.precision(std::numeric_limits<double>::max_digits10);
	
	std::cout << "V0(0 , 0)     : " << model.V0(origin) << std::endl;
	std::cout << "V0(vh, vs)    : " << model.V0(vev) << std::endl;
	std::cout << "V(0 , 0 , 0)  : " << model.V(origin, 0.) << std::endl;
	std::cout << "V(vh, vs, 0)  : " << model.V(vev, 0.) << std::endl;
	std::cout << "V(0 , 0 , 100): " << model.V(origin, 100.) << std::endl;
	std::cout << "V(vh, vs, 100): " << model.V(vev, 100.) << std::endl;

	PhaseTracer::PhaseFinder pf(model);
	pf.set_t_high(300);
	pf.set_check_vacuum_at_high(false);
	std::cout << "before find_phases()" << std::endl;
	pf.find_phases();
	std::cout << "after find_phases()" << std::endl;
	std::cout << pf;

	PhaseTracer::TransitionFinder tf(pf);
	tf.find_transitions();
	std::cout << tf;

	if(debug_mode)
	{
		std::string fileName = "xSM_MSbar_noZ2";
		PhaseTracer::potential_plotter(model, tf.get_transitions().front().TC, fileName, 0., 2., 0.01, -2, 0., 0.01);
		PhaseTracer::potential_line_plotter(model, tf.get_transitions(), fileName);
		PhaseTracer::phase_plotter(tf, fileName);
	}
}