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
        /** Reference to phase for EoS calculations */
        const Phase& phase;

        /** Thermodynamic quantities and splines */
        std::vector<double> potential_values;
        std::vector<double> temperature_values;
        alglib::spline1dinterpolant potential_spline;
        std::vector<double> temperature;
        std::vector<double> pressure;
        std::vector<double> entropy;
        std::vector<double> energy;
        std::vector<double> enthalpy;

        /** Temperature range and other parameters */
        double t_min, t_max;
        int n_temp;
        double background_dof;
        double energy_norm;

    public:
        
        /** Thermodynamic splines */
        alglib::spline1dinterpolant pressure_spline;
        alglib::spline1dinterpolant energy_spline;
        alglib::spline1dinterpolant enthalpy_spline;
        alglib::spline1dinterpolant entropy_spline;

        /** @brief Returns the potential spline for this phase. */
        const alglib::spline1dinterpolant& get_potential_spline() const { return potential_spline; }

        EquationOfStateInPhase(const Phase& phase_in, int n_temp_in, double background_dof_in, double energy_norm_in = 0.0) :
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

        /** @brief Creates splines for the EoS */
        void get_thermodynamic_splines();

    }; // class EquationOfStateInPhase

public:

    // TODO is this needed?
    EquationOfState() = default;
    
    // Delete copy constructor and copy assignment to prevent shallow copies of ALGLIB splines
    EquationOfState(const EquationOfState&) = delete;
    EquationOfState& operator=(const EquationOfState&) = delete;
    
    // Allow move semantics
    EquationOfState(EquationOfState&&) = default;
    EquationOfState& operator=(EquationOfState&&) = default;

    explicit EquationOfState(const Transition& transition_in) :
    transition(transition_in), t_min(transition_in.false_phase.T.front()), t_max(transition_in.TC) {}

    /**
     * @brief Builds the pressure, energy, enthalpy and entropy splines in both phases, along with the potential splines.
     */
    void calculate();

    /** @brief Whether calculate() has completed successfully. */
    bool is_calculated() const { return calculated; }

    /** 
     * @brief Evaluates the potential in the false vacuum.
     * @param T Temperature at which to evaluate the potential
     * @return Array containing V, dVdT, and d2VdT
     */
    std::array<double, 3> eval_false_potential(double T) const;

    /** 
     * @brief Evaluates the potential in the true vacuum.
     * @param T Temperature at which to evaluate the potential
     * @return Array containing V, dVdT, and d2VdT
     */
    std::array<double, 3> eval_true_potential(double T) const;

    /** 
     * @brief Evaluates the energy density.
     * @param T Temperature at which to evaluate the energy density.
     * @return Pair containing energy density in the false and true vacua.
     */
    std::pair<double, double> get_energy(double T) const;

    /** 
     * @brief Evaluates the pressure density.
     * @param T Temperature at which to evaluate the pressure.
     * @return Pair containing pressure in the false and true vacua.
     */
    std::pair<double, double> get_pressure(double T) const;

    /** 
     * @brief Evaluates the enthalpy density.
     * @param T Temperature at which to evaluate the enthalpy.
     * @return Pair containing enthalpy in the false and true vacua.
     */
    std::pair<double, double> get_enthalpy(double T) const;

    /** 
     * @brief Evaluates the entropy density.
     * @param T Temperature at which to evaluate the entropy.
     * @return Pair containing entropy in the false and true vacua.
     */
    std::pair<double, double> get_entropy(double T) const;

    /**
     * @brief Evaluates the first and second derivatives of the energy density.
     * @param T Temperature at which to evaluate the derivatives.
     * @return Pair of vectors containing the first and second derivatives of the energy density in the false and true vacua.
     */
    std::pair<std::vector<double>, std::vector<double>> get_energy_derivs(double T) const;

    /**
     * @brief Evaluates the first and second derivatives of the pressure density.
     * @param T Temperature at which to evaluate the derivatives.
     * @return Pair of vectors containing the first and second derivatives of the pressure density in the false and true vacua.
     */
    std::pair<std::vector<double>, std::vector<double>> get_pressure_derivs(double T) const;

    /**
     * @brief Evaluates the sound speed in both phases.
     * @param T Temperature at which to evaluate the sound speed.
     * @return Pair containing the sound speed in the false and true vacua.
     */
    std::pair<double, double> get_sound_speed(double T) const;

    /**
     * @brief Evaluates the trace anomaly in both phases.
     * @param T Temperature at which to evaluate the trace anomaly.
     * @param use_munu Whether to use the mu-nu definition of the trace anomaly.
     * @return Pair containing the trace anomaly in the false and true vacua.
     */
    std::pair<double, double> get_theta(double T, bool use_munu = false) const;

    /**
     * @brief Evaluates the energy density in the false vacuum.
     * @param T Temperature at which to evaluate the energy density.
     * @return Energy density in the false vacuum.
     */
    double get_energy_plus(double T) const;

    /**
     * @brief Evaluates the energy density in the true vacuum.
     * @param T Temperature at which to evaluate the energy density.
     * @return Energy density in the true vacuum.
     */
    double get_energy_minus(double T) const;

    /**
     * @brief Evaluates the pressure in the false vacuum.
     * @param T Temperature at which to evaluate the pressure.
     * @return Pressure in the false vacuum.
     */
    double get_pressure_plus(double T) const;

    /**
     * @brief Evaluates the pressure in the true vacuum.
     * @param T Temperature at which to evaluate the pressure.
     * @return Pressure in the true vacuum.
     */
    double get_pressure_minus(double T) const;

    /**
     * @brief Evaluates the enthalpy in the false vacuum.
     * @param T Temperature at which to evaluate the enthalpy.
     * @return Enthalpy in the false vacuum.
     */
    double get_enthalpy_plus(double T) const;

    /**
     * @brief Evaluates the enthalpy in the true vacuum.
     * @param T Temperature at which to evaluate the enthalpy.
     * @return Enthalpy in the true vacuum.
     */
    double get_enthalpy_minus(double T) const;

    /**
     * @brief Evaluates the entropy in the false vacuum.
     * @param T Temperature at which to evaluate the entropy.
     * @return Entropy in the false vacuum.
     */
    double get_entropy_plus(double T) const;

    /**
     * @brief Evaluates the entropy in the true vacuum.
     * @param T Temperature at which to evaluate the entropy.
     * @return Entropy in the true vacuum.
     */
    double get_entropy_minus(double T) const;

    /**
     * @brief Evaluates the first and second derivatives of the energy density in the false vacuum.
     * @param T Temperature at which to evaluate the derivatives.
     * @return Vector containing the first and second derivatives of the energy density in the false vacuum.
     */
    std::vector<double> get_energy_derivs_plus(double T) const;

    /**
     * @brief Evaluates the first and second derivatives of the energy density in the true vacuum.
     * @param T Temperature at which to evaluate the derivatives.
     * @return Vector containing the first and second derivatives of the energy density in the true vacuum.
     */
    std::vector<double> get_energy_derivs_minus(double T) const;

    /**
     * @brief Evaluates the first and second derivatives of the pressure density in the false vacuum.
     * @param T Temperature at which to evaluate the derivatives.
     * @return Vector containing the first and second derivatives of the pressure density in the false vacuum.
     */
    std::vector<double> get_pressure_derivs_plus(double T) const;

    /**
     * @brief Evaluates the first and second derivatives of the pressure density in the true vacuum.
     * @param T Temperature at which to evaluate the derivatives.
     * @return Vector containing the first and second derivatives of the pressure density in the true vacuum.
     */
    std::vector<double> get_pressure_derivs_minus(double T) const;

    /**
     * @brief Evaluates the sound speed in the false vacuum.
     * @param T Temperature at which to evaluate the sound speed.
     * @return Sound speed in the false vacuum.
     */
    double get_sound_speed_plus(double T) const;

    /**
     * @brief Evaluates the sound speed in the true vacuum.
     * @param T Temperature at which to evaluate the sound speed.
     * @return Sound speed in the true vacuum.
     */
    double get_sound_speed_minus(double T) const;

    /**
     * @brief Evaluates the trace anomaly in the false vacuum.
     * @param T Temperature at which to evaluate the trace anomaly.
     * @param use_munu Whether to use the mu-nu definition of the trace anomaly.
     * @return Trace anomaly in the false vacuum.
     */
    double get_theta_plus(double T, bool use_munu = true) const;

    /**
     * @brief Evaluates the trace anomaly in the true vacuum.
     * @param T Temperature at which to evaluate the trace anomaly.
     * @param use_munu Whether to use the mu-nu definition of the trace anomaly.
     * @return Trace anomaly in the true vacuum.
     */
    double get_theta_minus(double T, bool use_munu = true) const;

    /**
     * @brief Writes the equation of state to a CSV file.
     * @param path Path to the output CSV file.
     */
    void write(const std::string path) const;

private : 

    /** Reference to transition for which the equation of state is computed */
    const Transition& transition;

    /** Minimum and maximum temperatures for which the equation of state is computed. */
    PROPERTY(double, t_min, 0.0);
    PROPERTY(double, t_max, 0.0);

    /** Splines for the potential and equation of state */
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

    /** Number of spline evaluations for building the splines */
    PROPERTY(int, n_temp, 200);

    /** Background degrees of freedom */
    PROPERTY(double, background_dof, 0.0);

    /** Energy normalisation for the EoS */
    PROPERTY(std::optional<double>, energy_norm, {});

    /** Set by calculate() once the splines are built */
    bool calculated = false;

    /** Finds the normalisation from the zero-temp true vacuum energy */
    double find_normalisation(Phase true_vacuum);

    /** Checks if the temperature is within the valid range */
    void 
    check_temperature_range(double T, const char* caller) const
    {
        require_calculated(caller);
        if (T < t_min || T > t_max)
        {
            throw std::out_of_range(std::string("Temperature out of interpolation bounds in ") + caller);
        }
    }

    /** Checks if the equation of state has been calculated */
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