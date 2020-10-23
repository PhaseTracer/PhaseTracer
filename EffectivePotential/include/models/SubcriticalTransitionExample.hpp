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

#ifndef POTENTIAL_SubcriticalTransitionExample_HPP_INCLUDED
#define POTENTIAL_SubcriticalTransitionExample_HPP_INCLUDED

#include "potential.hpp"
#include "property.hpp"
#include "pow.hpp"


namespace EffectivePotential
{
class SubcriticalTransitionExample : public Potential
{
	public:
	/*double V(Eigen::VectorXd x, double T) const override
	{
		return (c * square(T) - m2) * square(x[0]) + kappa * cube(x[0]) + lambda * pow_4(x[0]);
	}*/

	double V(Eigen::VectorXd phi, double T) const override
	{
		const double x = phi[0];
		const double x2 = x*x;
		const double x3 = x2*x;
		const double x4 = x3*x;
		const double x5 = x4*x;
		const double x6 = x5*x;
		const double T2 = T*T;

		return -2. + T
			+ (-12. + 9.025*T - 2.45*T2)*x
			+ (24.275 - 16.1417*T + 3.28333*T2)*x2
			+ (-22. + 11.5111*T - 0.755556*T2)*x3
			+ (9.6625 - 3.16458*T - 0.770833*T2)*x4
			+ (-2.04 + 0.193333*T + 0.413333*T2)*x5
			+ (0.166667 + 0.0277778*T - 0.0555556*T2)*x6;
	}
	
	size_t get_n_scalars() const override { return 1; }
	
	bool forbidden(Eigen::VectorXd x) const override { return x[0] < -0.1; }
	
	/*Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const override
	{
		return 4. * c * T * phi;
	}*/
	
	/*Eigen::MatrixXd d2V_dx2(Eigen::VectorXd phi, double T) const  override
	{
		Eigen::MatrixXd d(1, 1);
		d(0, 0) = 2. * (c * square(T) - m2) + 6. * kappa * phi[0] + 12. * lambda * square(phi[0]);
		return d;
	}*/
}; // class SubcriticalTransitionExample
}  // namespace EffectivePotential

#endif
