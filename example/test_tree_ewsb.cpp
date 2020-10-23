#include "models/xSM_MSbar.hpp"
#include "models/xSM_MSbar_solver.hpp"

#include <iostream>
#include <iomanip>
#include <Eigen/Eigenvalues>

int main(int argc, char* argv[])
{
	double vh = 246.221;
	auto model = EffectivePotential::xSM_MSbar(0.5, vh, false);

	double V0, V1, V2;
	double step = 0.1;

	Eigen::VectorXd vev(2);
	vev << vh+step, 0.;
	V0 = model.V(vev, 0.);
	vev << vh, 0.;
	V1 = model.V(vev, 0.);
	vev << vh-step, 0.;
	V2 = model.V(vev, 0.);

	double mhSq = (V0 - 2*V1 + V2) / (step*step);

	std::cout.precision(std::numeric_limits<double>::max_digits10);
	std::cout << "tree_ewsb=false: " << mhSq << std::endl;

	auto model2 = EffectivePotential::xSM_MSbar(0.5, vh, true);

	vev << vh+step, 0.;
	V0 = model2.V(vev, 0.);
	vev << vh, 0.;
	V1 = model2.V(vev, 0.);
	vev << vh-step, 0.;
	V2 = model2.V(vev, 0.);

	mhSq = (V0 - 2*V1 + V2) / (step*step);

	std::cout << "tree_ewsb=true: " << mhSq << std::endl;
}