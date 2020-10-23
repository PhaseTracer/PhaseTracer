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

#ifndef POTENTIAL_XSM_MSbar_noZ2_HPP_INCLUDED
#define POTENTIAL_XSM_MSbar_noZ2_HPP_INCLUDED

//TODO:
/**
  * Set check_vacuum_at_high to false, otherwise we'll get errors for minima away from the origin at T>>Tc (which
  * happens).
  */

/**
 * Real scalar singlet extension of the Standard Model, not assuming Z2 symmetry.
 */

#include "one_loop_potential.hpp"
#include "pow.hpp"
#include "SM_parameters.hpp"

#include <vector>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>


namespace EffectivePotential
{
class xSM_MSbar_noZ2 : public OneLoopPotential
{
	public:
	/**
	* @brief Make an xSM model (with no Z2 symmetry) from Lagrangian parameters.
	*/
	xSM_MSbar_noZ2(double muh_sq_, double mus_sq_, double lh_, double ls_, double c1_, double c2_, double b3_, double
		muh_sq_tree_EWSB_, double mus_sq_tree_EWSB_, double vs_, double ms_, double theta_) :
		muh_sq(muh_sq_), mus_sq(mus_sq_), lh(lh_), ls(ls_), c1(c1_), c2(c2_), b3(b3_),
		muh_sq_tree_EWSB(muh_sq_tree_EWSB_), mus_sq_tree_EWSB(mus_sq_tree_EWSB_), vs(vs_), theta_vev(theta_)
	{
		ms2 = ms_*ms_;

		cos2Theta_vev = std::cos(2*theta_vev);
		sin2Theta_vev = std::sin(2*theta_vev);
		cosSqTheta_vev = std::cos(theta_vev);
		cosSqTheta_vev *= cosSqTheta_vev;
		sinSqTheta_vev = 1 - cosSqTheta_vev;

		init();
	}

	/**
	  * Default constructor. Use at your own risk. There are no checks when using the potential that it has been
	  * properly initialised with parameter values. If using this default constructor, it is assumed you will Use
	  * setParameters afterwards.
	  */
	xSM_MSbar_noZ2()
	{
		init();
	}

	void init()
	{
		mh2 = 125.2*125.2;
		vh = 246.221;	

		massShiftHack = 1.;
		// Setting this to 1e-10 causes isses where the Goldstone mass can be just above 1e-10 and then we include the
		// massShiftHack term when avoiding the divergence from the Goldstone at the VEV.
		minMassThreshold = 1e-5;

		// Store the squared gauge and Yukawa couplings for use in the mass functions.
		/*g2 = SM::g;
		g2 = g2*g2;
		gp2 = SM::gp;
		gp2 = gp2*gp2;
		yt2 = SM::yt_sq;*/

		g2 = 0.650809745*0.650809745;
		gp2 = 0.357606440*0.357606440;
		yt2 = 0.997918833;
	}

	/** Used to consider a new potential without creating a new model object. */
	void setParameters(double muh_sq_, double mus_sq_, double lh_, double ls_, double c1_, double c2_, double b3_,
		double muh_sq_tree_EWSB_, double mus_sq_tree_EWSB_, double vs_, double ms_, double theta_)
	{
		muh_sq = muh_sq_;
		mus_sq = mus_sq_;
		lh = lh_;
		ls = ls_;
		c1 = c1_;
		c2 = c2_;
		b3 = b3_;
		muh_sq_tree_EWSB = muh_sq_tree_EWSB_;
		mus_sq_tree_EWSB = mus_sq_tree_EWSB_;
		vs = vs_;
		ms2 = ms_*ms_;
		theta_vev = theta_;

		cos2Theta_vev = std::cos(2*theta_vev);
		sin2Theta_vev = std::sin(2*theta_vev);
		cosSqTheta_vev = std::cos(theta_vev);
		cosSqTheta_vev *= cosSqTheta_vev;
		sinSqTheta_vev = 1 - cosSqTheta_vev;
	}

	void setParameters(std::vector<double> params)
	{
		if(params.size() != 12)
		{
			std::cerr << "Called xSM_MSbar_noZ2::setParameters with the wrong number of parameters: got " <<
				params.size() << ", expected 12.";
			return;
		}

		muh_sq = params[5];
		mus_sq = params[6];
		lh = params[7];
		ls = params[8];
		c1 = params[0];
		c2 = params[9];
		b3 = params[1];
		muh_sq_tree_EWSB = params[10];
		mus_sq_tree_EWSB = params[11];
		vs = params[3];
		ms2 = params[4]*params[4];
		theta_vev = params[2];

		cos2Theta_vev = std::cos(2*theta_vev);
		sin2Theta_vev = std::sin(2*theta_vev);
		cosSqTheta_vev = std::cos(theta_vev);
		cosSqTheta_vev *= cosSqTheta_vev;
		sinSqTheta_vev = 1 - cosSqTheta_vev;
	}

	void printParameters()
	{
		auto precisionPrev = std::cout.precision(std::numeric_limits<double>::max_digits10);

		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << "muh_sq: " << muh_sq << std::endl;
		std::cout << "mus_sq: " << mus_sq << std::endl;
		std::cout << "lh: " << lh << std::endl;
		std::cout << "ls: " << ls << std::endl;
		std::cout << "c1: " << c1 << std::endl;
		std::cout << "c2: " << c2 << std::endl;
		std::cout << "b3: " << b3 << std::endl;
		std::cout << "muh_sq_tree_EWSB: " << muh_sq_tree_EWSB << std::endl;
		std::cout << "mus_sq_tree_EWSB: " << mus_sq_tree_EWSB << std::endl;
		std::cout << "vs: " << vs << std::endl;
		std::cout << "ms: " << std::sqrt(ms2) << std::endl;
		std::cout << "theta: " << theta_vev << std::endl;
		std::cout << "------------------------------------------------------" << std::endl;

		std::cout.precision(precisionPrev);
	}
	
	/*xSM_MSbar_noZ2(double c1_, double b3_, double theta_, double vs_, double ms_, double Q) :
		c1(c1_), b3(b3_), theta(theta_), vs(vs_)
	{
		mh2 = 125.2*125.2;
		vh = 246.221;
		ms2 = ms_*ms_;

		cos2Theta = std::cos(2*theta);
		sin2Theta = std::sin(2*theta);
		cosSqTheta = std::cos(theta);
		cosSqTheta *= cosSqTheta;
		sinSqTheta = std::sin(theta);
		sinSqTheta *= sinSqTheta;

		valid = true;
	}*/

	Eigen::VectorXd get_EW_VEV() const
	{
		 Eigen::VectorXd vev(2);
		 vev << vh, vs;

		 return vev;
	}

	double get_vh()
	{
		return vh;
	}

	double V0(Eigen::VectorXd phi) const override
	{
		double h = phi[0];
		double s = phi[1];
		double h2 = h*h;
		double s2 = s*s;

		return muh_sq*h2 + mus_sq*s2 + lh*h2*h2 + ls*s2*s2 + c1*h2*s + c2*h2*s2 + b3*s2*s;
	}

	std::vector<double> get_scalar_thermal_sq(double T) const override
	{
		const double T2 = T*T;
		const double debyeHiggs = 3./16.*g2 + gp2/16. + yt2/4. + 2*lh + c2/6.;
		const double debyeSinglet = 2./3.*c2 + ls;
		const double debyeGoldstone = debyeHiggs;

		return {debyeHiggs*T2, debyeSinglet*T2, debyeGoldstone*T2};
	}

	std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const override
	{
		// TODO: Need to add xi dependence!
		double h = phi[0];
		double s = phi[1];
		double h2 = h*h;
		double s2 = s*s;

		std::vector<double> scalar_eigs = get_scalar_eigenvalues(phi, xi);

		double goldstone = 2*muh_sq_tree_EWSB + 4*lh*h2 + 2*c1*s + 2*c2*s2;

		return {scalar_eigs[0], scalar_eigs[1], goldstone};
	}

	std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const override
	{
		// TODO: Need to add xi dependence!
		double h = phi[0];
		double s = phi[1];
		double h2 = h*h;
		double s2 = s*s;

		std::vector<double> scalar_eigs = get_scalar_eigenvalues_thermal(phi, xi, T);

		double goldstone = 2*muh_sq_tree_EWSB + 4*lh*h2 + 2*c1*s + 2*c2*s2;
		goldstone += (3./16.*g2 + gp2/16. + yt2/4. + 2*lh + c2/6.)*T*T;

		return {scalar_eigs[0], scalar_eigs[1], goldstone};
	}

	std::vector<double> get_scalar_eigenvalues(Eigen::VectorXd phi, double xi) const
	{
		// TODO: Need to add xi dependence! (is there xi dependence in the eigenvalues???)
		double mhh = d2V0mdh2(phi);
		double mhs = d2V0mdhds(phi);
		double mss = d2V0mds2(phi);
		double theta = get_mixing_angle_supplied(mhh, mhs, mss);

		double cosSqTheta = std::cos(theta);
		cosSqTheta *= cosSqTheta;
		double sinSqTheta = 1 - cosSqTheta;
		double sin2Theta = std::sin(2*theta);

		double mh = mhh*cosSqTheta + mhs*sin2Theta + mss*sinSqTheta;
		double ms = mhh*sinSqTheta - mhs*sin2Theta + mss*cosSqTheta;

		return {mh, ms};
	}

	std::vector<double> get_scalar_eigenvalues_thermal(Eigen::VectorXd phi, double xi, double T) const
	{
		auto thermal_corrections = get_scalar_thermal_sq(T);
		
		// TODO: Need to add xi dependence! (is there xi dependence in the eigenvalues???)
		double mhh = d2V0mdh2(phi) + thermal_corrections[0];
		double mhs = d2V0mdhds(phi);
		double mss = d2V0mds2(phi) + thermal_corrections[1];
		double theta = get_mixing_angle_supplied(mhh, mhs, mss);

		double cosSqTheta = std::cos(theta);
		cosSqTheta *= cosSqTheta;
		double sinSqTheta = 1 - cosSqTheta;
		double sin2Theta = std::sin(2*theta);

		double mh = mhh*cosSqTheta + mhs*sin2Theta + mss*sinSqTheta;
		double ms = mhh*sinSqTheta - mhs*sin2Theta + mss*cosSqTheta;

		return {mh, ms};
	}

	double get_mixing_angle_supplied(double mhh, double mhs, double mss) const
	{
		return 0.5*std::atan2(2*mhs, mhh - mss);
	}

	double get_mixing_angle_tree(Eigen::VectorXd phi) const
	{
		double mhh = d2V0mdh2(phi);
		double mhs = d2V0mdhds(phi);
		double mss = d2V0mds2(phi);

		return get_mixing_angle_supplied(mhh, mhs, mss);
	}

	std::vector<double> get_scalar_dofs() const override
	{
		return {1., 1., 3.};
	}

	std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const override
	{
		const double h2 = phi[0]*phi[0];
		const double mWSq = 0.25*g2*h2;
		const double mZSq = 0.25*(g2 + gp2)*h2;
		const double mPhSq = 0;
		
		return {mWSq, mZSq, mPhSq};
	}

	std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const override
	{
		const double h2 = phi[0]*phi[0];
		const double T2 = T*T;

		const double mWSq = 0.25*g2*h2 + 11./6.*g2*T2;
		//const double mZSq = 0.25*(g2 + gp2)*h2;
		//const double mPhSq = 0;
		
		const double a = (g2 + gp2)*(3*h2 + 22*T2);
		const double b = std::sqrt(9*square(g2 + gp2)*h2*h2 + 132*square(g2 - gp2)*h2*T2 + 484*square(g2 - gp2)
			*T2*T2);

		const double mZSq = (a + b)/24.;
		const double mPhSq = (a - b)/24.;

		return {mWSq, mZSq, mPhSq};
	}

	std::vector<double> get_vector_dofs() const override
	{
		return {6., 3., 2.};
	}

	std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const override
	{
		return {0.5*yt2*phi[0]*phi[0]};
	}

	std::vector<double> get_fermion_dofs() const override
	{
		return {12.};
	}

	size_t get_n_scalars() const override
	{
		return 2;
	}

	std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override
	{
		/*auto phi1 = phi;
		phi1[0] = -phi[0];
		auto phi2 = phi;
		phi2[1] = -phi[1];
		
		return {phi1, phi2};*/

		auto negH = phi;
		negH[0] = -phi[0];

		return {negH};
	}

	// TODO: this is ready, but the xi-dependent derivatives of the full potential need to be implemented
	// TODO: first. We never need to calculate the one-loop mixing angle, as we know it should equal the
	// TODO: theta that we inputted.
	/*double get_mixing_angle_1l(Eigen::VectorXd phi, double xi)
	{
		double mhh = d2Vdh2(phi, xi);
		double mhs = d2Vdhds(phi, xi);
		double mss = d2Vds2(phi, xi);

		return get_mixing_angle_supplied(mhh, mhs, mss);
	}*/

	/* ----------------------------------------------------------------------------------------------------
	 * Derivatives of the tree-level potential up to fourth order usiung tree-level EWSB to calculate the
	 * quadratic parameters from the quartic parameters (i.e. mu{h,s}_sq from l{h,s}). The 'subscript' 'm'
	 * in V0m signifies these forms are for use in masses that enter the Coleman-Weinberg potential, where
	 * we want the parameters to be related through EWSB at tree-level only.
	 * --------------------------------------------------------------------------------------------------*/

	// TODO: actually for the moment we only need them up to second order, since we're not writing the
	// TODO: solver yet. The only time we need these for general use of the potential is in calculating the
	// TODO: scalar mass eigenvalues.

	double dV0mdh(Eigen::VectorXd phi) const
	{
		const double h = phi[0];
		const double s = phi[1];

		return 2*muh_sq_tree_EWSB*h + 4*lh*h*h*h + 2*c1*h*s + 2*c2*h*s*s;
	}

	double dV0mds(Eigen::VectorXd phi) const
	{
		const double h = phi[0];
		const double s = phi[1];

		return 2*mus_sq_tree_EWSB*s + 4*ls*s*s*s + c1*h*h + 2*c2*h*h*s + 3*b3*s*s;
	}

	double d2V0mdh2(Eigen::VectorXd phi) const
	{
		const double h = phi[0];
		const double s = phi[1];

		return 2*muh_sq_tree_EWSB + 12*lh*h*h + 2*c1*s + 2*c2*s*s;
	}

	double d2V0mds2(Eigen::VectorXd phi) const
	{
		const double h = phi[0];
		const double s = phi[1];

		return 2*mus_sq_tree_EWSB + 12*ls*s*s + 2*c2*h*h + 6*b3*s;
	}

	double d2V0mdhds(Eigen::VectorXd phi) const
	{
		const double h = phi[0];
		const double s = phi[1];

		return 2*c1*h + 4*c2*h*s;
	}

	private:
	double muh_sq;
	double mus_sq;
	double lh;
	double ls;
	double c1;
	double c2;
	double b3;
	double muh_sq_tree_EWSB;
	double mus_sq_tree_EWSB;
	double vh;
	double vs;
	double mh2;
	double ms2;
	double theta_vev;

	double cosSqTheta_vev;
	double sinSqTheta_vev;
	double cos2Theta_vev;
	double sin2Theta_vev;

	// TODO: these are only needed in the second derivatives of VCW. Perhaps we don't need these variables in this
	// TODO: class, and can move them to the solver (when it has been implemented).
	double massShiftHack;
	double minMassThreshold;

	double g2;
	double gp2;
	double yt2;
}; // class xSM_MSbar_noZ2
}  // namespace EffectivePotential

#endif
