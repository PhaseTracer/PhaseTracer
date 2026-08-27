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

#ifndef PHASETRACER_EQUATION_OF_STATE_HPP_
#define PHASETRACER_EQUATION_OF_STATE_HPP_

#include <cmath>
#include <array>
#include <vector>
#include <optional>
#include <stdexcept>
#include <interpolation.h>

#include "property.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"

namespace PhaseTracer {

class EquationOfState 
{
private:
    class EquationOfStateInPhase
    {
    private:
        Phase phase;
        std::vector<double> potential_values;
        std::vector<double> temperature_values;
        alglib::spline1dinterpolant potential_spline;
        std::vector<double> temperature;
        std::vector<double> pressure;
        std::vector<double> entropy;
        std::vector<double> energy;
        std::vector<double> enthalpy;
        double t_min, t_max;
        int n_temp;
        double background_dof;
        double energy_norm;

    public:
        alglib::spline1dinterpolant pressure_spline;
        alglib::spline1dinterpolant energy_spline;
        alglib::spline1dinterpolant enthalpy_spline;
        alglib::spline1dinterpolant entropy_spline;

        const alglib::spline1dinterpolant& get_potential_spline() const { return potential_spline; }

        EquationOfStateInPhase(Phase phase_in, int n_temp_in, double background_dof_in, double energy_norm_in = 0.0) :
        phase(phase_in),
        potential_values(phase_in.V),
        temperature_values(phase_in.T),
        t_min(phase_in.T.front()),
        t_max(phase_in.T.back()),
        n_temp(n_temp_in),
        background_dof(background_dof_in),
        energy_norm(energy_norm_in)
        {
            alglib::real_1d_array t, v;
            v.setcontent(potential_values.size(), potential_values.data());
            t.setcontent(temperature_values.size(), temperature_values.data());
            alglib::spline1dbuildcubic(t, v, this->potential_spline);

            temperature.resize(n_temp);
            double dT = std::abs(t_max - t_min) / (n_temp - 1);
            for (int i = 0; i < n_temp; ++i) 
            {
                temperature[i] = t_min + i * dT;
            }

            if (abs(temperature.front() - t_min) > 1e-6 || abs(temperature.back() - t_max) > 1e-6) {
                throw std::runtime_error("Temperature vector does not match phase temperature bounds.");
            }

            if (temperature.size() != n_temp) {
                throw std::runtime_error("Temperature vector size does not match n_temp.");
            }

            get_thermodynamic_splines();
        }

        void get_thermodynamic_splines();

    }; // class EquationOfStateInPhase

    Transition transition;
    double t_min = 0.0, t_max = 0.0;
    alglib::spline1dinterpolant false_potential_spline;
    alglib::spline1dinterpolant true_potential_spline;
    alglib::spline1dinterpolant p_plus_spline;
    alglib::spline1dinterpolant p_minus_spline;
    alglib::spline1dinterpolant e_plus_spline;
    alglib::spline1dinterpolant e_minus_spline;
    alglib::spline1dinterpolant w_plus_spline;
    alglib::spline1dinterpolant w_minus_spline;
    alglib::spline1dinterpolant s_plus_spline;
    alglib::spline1dinterpolant s_minus_spline;
    int n_temp = 200;
    double background_dof = 0.0;
    std::optional<double> energy_norm;
    bool calculated = false;

public:

    EquationOfState() = default;
    
    EquationOfState(const EquationOfState&) = delete;
    EquationOfState& operator=(const EquationOfState&) = delete;
    
    EquationOfState(EquationOfState&&) = default;
    EquationOfState& operator=(EquationOfState&&) = default;

    explicit EquationOfState(Transition transition_in) :
    transition(transition_in),
    t_min(transition_in.false_phase.T.front()),
    t_max(transition_in.TC)
    {}

    /**
     * @brief Builds the pressure, energy, enthalpy and entropy splines in both
     *        phases, along with the potential splines.
     *
     * Must be called before any accessor methods. Failing to do so will throw
     * a logic error.
     *
     * If no energy normalisation has been set via set_energy_norm(), it is
     * computed internally from the true vacuum via find_normalisation(). This 
     * keeps the Hubble rate well defined.
     */
    void calculate();

    /** @brief Whether calculate() has completed successfully. */
    bool is_calculated() const { return calculated; }

    /**
     * @brief Set the energy normalisation explicitly (negative by definition).
     *
     * Use when the true phase is not the zero-temperature ground state -- e.g.
     * an intermediate transition -- so the internally derived normalisation is
     * not the one you want. Leave unset to derive it from the true vacuum.
     */
    void set_energy_norm(double energy_norm_in) { energy_norm = energy_norm_in; }

    /** @brief The explicitly set energy normalisation, if any. */
    std::optional<double> get_energy_norm() const { return energy_norm; }

    void set_n_temp(int n_temp_in) { n_temp = n_temp_in; }

    void set_background_dof(double background_dof_in) { background_dof = background_dof_in; }

    double get_background_dof() const { return background_dof; }


    std::array<double, 3> eval_false_potential(double T) const;
    std::array<double, 3> eval_true_potential(double T) const;

    std::pair<double, double> get_energy(double T) const;
    std::pair<double, double> get_pressure(double T) const;
    std::pair<double, double> get_enthalpy(double T) const;
    std::pair<double, double> get_entropy(double T) const;

    std::pair<std::vector<double>, std::vector<double>> get_energy_derivs(double T) const;
    std::pair<std::vector<double>, std::vector<double>> get_pressure_derivs(double T) const;

    std::pair<double, double> get_sound_speed(double T) const;
    std::pair<double, double> get_theta(double T, bool use_munu = false) const;

    double get_energy_plus(double T) const;
    double get_energy_minus(double T) const;
    double get_pressure_plus(double T) const;
    double get_pressure_minus(double T) const;
    double get_enthalpy_plus(double T) const;
    double get_enthalpy_minus(double T) const;
    double get_entropy_plus(double T) const;
    double get_entropy_minus(double T) const;

    std::vector<double> get_energy_derivs_plus(double T) const;
    std::vector<double> get_energy_derivs_minus(double T) const;
    std::vector<double> get_pressure_derivs_plus(double T) const;
    std::vector<double> get_pressure_derivs_minus(double T) const;

    double get_sound_speed_plus(double T) const;
    double get_sound_speed_minus(double T) const;
    double get_theta_plus(double T) const;
    double get_theta_minus(double T) const;

    void write(const std::string path) const;

    double get_t_min() const { return t_min; }
    double get_t_max() const { return t_max; }
    int get_n_temp() const { return n_temp; }

private : 

    double find_normalisation(Phase true_vacuum);

    void
    check_temperature_range(double T, const char* caller) const
    {
        require_calculated(caller);
        if (T < t_min || T > t_max)
        {
            throw std::out_of_range(std::string("Temperature out of interpolation bounds in ") + caller);
        }
    }

    void
    require_calculated(const char* caller) const
    {
        if (!calculated)
        {
            throw std::logic_error(std::string("EquationOfState::") + caller + " called before calculate().");
        }
    }

}; // class EquationOfState

} // namespace PhaseTracer


#endif // PHASETRACER_EQUATION_OF_STATE_HPP_