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

#ifndef PHASETRACER_TRANSITION_METRICS_HPP_
#define PHASETRACER_TRANSITION_METRICS_HPP_

#include <cmath>
#include <vector>
#include <optional>
#include <stdexcept>
#include <interpolation.h>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "false_vacuum_decay_rate.hpp"
#include "equation_of_state.hpp"

namespace PhaseTracer {

enum SomethingStatus
{
    YES,
    FAST,
    NO,
    INVALID,
    ERR
};

enum MilestoneType
{
    PERCOLATION,
    NUCLEATION,
    COMPLETION
};

struct TransitionMilestone
{
    MilestoneType type;
    SomethingStatus status;
    double temperature;

    TransitionMilestone(const MilestoneType& type_in)
    : type(type_in), status(SomethingStatus::ERR), temperature(0.0) {}
};

class TransitionMetrics {

    FalseVacuumDecayRate decay_rate;

    EquationOfState eos;

    double t_min, t_max;

    alglib::spline1dinterpolant a2a1_integrand_spline;

    alglib::spline1dinterpolant Pf_spline;

    PROPERTY(double, total_number_temp_steps, 200);

    PROPERTY(double, vw, 0.577);

    PROPERTY(double, dof, 106.75);

    PROPERTY(double, newtonG,  1/((1.22 * 1e19)*(1.22 * 1e19)));

public :

    TransitionMetrics(const FalseVacuumDecayRate& decay_rate_in, const EquationOfState& eos_in) :
    decay_rate(decay_rate_in), eos(eos_in), t_min(decay_rate_in.get_t_min()), t_max(decay_rate_in.get_t_max()) 
    {
        make_a2a1_integrand_spline();

        compute_Pf_spline();
    }

    const double get_hubble_rate(const double& T);

    const double get_dtdT(const double& T);

    const double get_atop_abottom(const double& Ttop, const double& Tbottom);

    const double get_extended_volume(const double& T);

    const double get_nucleation_rate(const double& T);

    const double find_pf_temperature(const double& target);

    const double find_nucleation_temperature();

private:

    void make_a2a1_integrand_spline();

    const double get_volume_term(const double& T1, const double& T2);

    const double extended_volume_integrand(const double& T1, const double& T2);

    void compute_Pf_spline();
};


} // namespace PhaseTracer

#endif // PHASETRACER_TRANSITION_METRICS_HPP_