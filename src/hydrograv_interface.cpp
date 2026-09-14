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

#ifdef BUILD_WITH_HG

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "logger.hpp"
#include "hydrograv_interface.hpp"

namespace PhaseTracer {
namespace HydroGravBridge {

HydroGrav::PhaseTransition::EquationOfState
to_hydrograv_eos(const PhaseTracer::EquationOfState &eos)
{
    const double t_min = eos.get_t_min();
    const double t_max = eos.get_t_max();
    const int n_temp = eos.get_n_temp();
    const double dt = (t_max - t_min) / (n_temp - 1);

    std::vector<double> t_vals, p_plus_vals, p_minus_vals, e_plus_vals, e_minus_vals;

    for (int i = 0; i < n_temp; ++i)
    {
        double T = t_min + i * dt;
        auto [p_plus, p_minus] = eos.get_pressure(T);
        auto [e_plus, e_minus] = eos.get_energy(T);

        if(p_plus < 0 || p_minus < 0 || e_plus < 0 || e_minus < 0) { continue; }

        t_vals.push_back(T);
        p_plus_vals.push_back(p_plus);
        p_minus_vals.push_back(p_minus);
        e_plus_vals.push_back(e_plus);
        e_minus_vals.push_back(e_minus);
    }

    return HydroGrav::PhaseTransition::EquationOfState(t_vals, p_plus_vals, p_minus_vals, e_plus_vals, e_minus_vals);
}

HydroGrav::PhaseTransition::Universe
to_hydrograv_universe(const TransitionMilestone &milestone, double dof)
{
    if (milestone.status != MilestoneStatus::YES)
    {
        return HydroGrav::PhaseTransition::Universe();
    }

    if (milestone.temperature <= 0 || milestone.H <= 0)
    {
        throw std::invalid_argument("Invalid thermal parameters for universe creation");
    }

    return HydroGrav::PhaseTransition::Universe(milestone.temperature, dof, milestone.H);
}

const char *
to_hydrograv_nuc_type(NucleationType nucleation_type)
{
    switch (nucleation_type)
    {
        case NucleationType::SIMULTANEOUS: return "sim";
        case NucleationType::EXPONENTIAL:  return "exp";
    }
    return "exp";
}

HydroGrav::PhaseTransition::PTParams_Veff
to_hydrograv_pt_params(const ThermalParameterSet &tps,
                       const TransitionMilestone &milestone,
                       double vw, double dof)
{
    return HydroGrav::PhaseTransition::PTParams_Veff(
        vw,
        milestone.alpha_munu,
        milestone.temperature,
        milestone.betaH_eff * milestone.H,
        milestone.Rs / milestone.H,
        to_hydrograv_nuc_type(milestone.nucleation_type),
        to_hydrograv_universe(milestone, dof),
        to_hydrograv_eos(tps.get_equation_of_state())
    );
}

FluidProfile
to_phasetracer_profile(const HydroGrav::Hydrodynamics::FluidProfile &fp)
{
    FluidProfile out;

    out.xi = fp.xi_vals();
    out.v = fp.v_vals();
    out.w = fp.w_vals();
    out.lambda = fp.la_vals();
    out.T = fp.T_vals();

    out.xi_min = fp.xi_min();
    out.xi_max = fp.xi_max();
    out.cs_plus_sq = fp.cpsq();
    out.cs_minus_sq = fp.cmsq();
    out.mode = fp.mode();
    out.shock_converged = fp.shock_flag();

    return out;
}

void
fill_spectrum(GravWaveSpectrum &sp, const HydroGrav::Spectrum::PowerSpec &spec_in)
{
    sp.method = GravWaveMethod::SoundShell;

    sp.frequency = spec_in.freq();
    sp.kRs = spec_in.K();
    sp.sound_wave = spec_in.P();
    sp.total_amplitude = spec_in.P();
    sp.dtau = spec_in.dtau();

    sp.turbulence = std::vector<double>(sp.frequency.size(), 0.0);
    sp.bubble_collision = std::vector<double>(sp.frequency.size(), 0.0);

    const auto &peaks = spec_in.peak_vals();
    sp.peak_frequency = peaks.first;
    sp.peak_amplitude = peaks.second;

    sp.profile = to_phasetracer_profile(spec_in.profile());

    sp.SNR.push_back(LISA_snr(sp.frequency, sp.total_amplitude));
    sp.SNR.push_back(std::numeric_limits<double>::quiet_NaN());
}

} // namespace HydroGravBridge
} // namespace PhaseTracer

#endif // BUILD_WITH_HG
