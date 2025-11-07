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

#include <cmath>
#include "logger.hpp"
#include "equation_of_state.hpp"

namespace PhaseTracer {

    void 
    EquationOfState::EquationOfStateInPhase::get_thermodynamic_splines() 
    {
        pressure.resize(n_temp);
        energy.resize(n_temp);
        enthalpy.resize(n_temp);
        entropy.resize(n_temp);

        for (int i = 0; i < n_temp; ++i) 
        {
            double temp = temperature[i];
            double v, dvdT, ddvdT;
            try {
                alglib::spline1ddiff(this->potential_spline, temp, v, dvdT, ddvdT);
            } catch (const std::exception& e) {
                LOG(error) << "Error in spline1ddiff: " << e.what() << "for temperature " << temp << std::endl;
                continue;
            }

            double temp2 = temp * temp;
            double temp4 = temp2 * temp2;
            double background_pressure = background_dof * M_PI * M_PI / 90.0 * temp4;

            pressure[i] = background_pressure - v;
            energy[i] = 3.0 * background_pressure + v - temp * dvdT;
            enthalpy[i] = 4.0 * background_pressure - temp * dvdT;
            entropy[i] = 4.0 * background_pressure / temp - dvdT;
        }

        alglib::real_1d_array t_array, p_array, e_array, w_array, s_array;
        t_array.setcontent(n_temp, temperature.data());
        p_array.setcontent(n_temp, pressure.data());
        e_array.setcontent(n_temp, energy.data());
        w_array.setcontent(n_temp, enthalpy.data());
        s_array.setcontent(n_temp, entropy.data());

        try {
            alglib::spline1dbuildcubic(t_array, p_array, this->pressure_spline);
            alglib::spline1dbuildcubic(t_array, e_array, this->energy_spline);
            alglib::spline1dbuildcubic(t_array, w_array, this->enthalpy_spline);
            alglib::spline1dbuildcubic(t_array, s_array, this->entropy_spline);
            LOG(debug) << "EquationOfStateInPhase splines built.";
        } catch (const std::exception& e) {
            LOG(error) << "Error in spline1dbuildcubic: " << e.what() << std::endl;
            throw std::runtime_error("Failed to build thermodynamic splines");
        }
    }

    std::pair<double, double>
    EquationOfState::get_energy(double T)  const
    {
        check_temperature_range(T, "get_energy");
        double e_plus = alglib::spline1dcalc(this->e_plus_spline, T);
        double e_minus = alglib::spline1dcalc(this->e_minus_spline, T);
        return {e_plus, e_minus};
    }

    std::pair<double, double>
    EquationOfState::get_pressure(double T)  const
    {
        check_temperature_range(T, "get_pressure");
        double p_plus = alglib::spline1dcalc(this->p_plus_spline, T);
        double p_minus = alglib::spline1dcalc(this->p_minus_spline, T);
        return {p_plus, p_minus};
    }

    std::pair<double, double>
    EquationOfState::get_enthalpy(double T)  const
    {
        check_temperature_range(T, "get_enthalpy");
        double w_plus = alglib::spline1dcalc(this->w_plus_spline, T);
        double w_minus = alglib::spline1dcalc(this->w_minus_spline, T);
        return {w_plus, w_minus};
    }

    std::pair<double, double>
    EquationOfState::get_entropy(double T)  const
    {
        check_temperature_range(T, "get_entropy");
        double s_plus = alglib::spline1dcalc(this->s_plus_spline, T);
        double s_minus = alglib::spline1dcalc(this->s_minus_spline, T);
        return {s_plus, s_minus};
    }

    std::pair<std::vector<double>, std::vector<double>>
    EquationOfState::get_energy_derivs(double T)  const
    {
        check_temperature_range(T, "get_energy_derivs");

        double e_plus, de_plus, dde_plus;
        double e_minus, de_minus, dde_minus;

        alglib::spline1ddiff(this->e_plus_spline, T, e_plus, de_plus, dde_plus);
        alglib::spline1ddiff(this->e_minus_spline, T, e_minus, de_minus, dde_minus);

        return {{e_plus, de_plus, dde_plus}, {e_minus, de_minus, dde_minus}};
    }

    std::pair<std::vector<double>, std::vector<double>>
    EquationOfState::get_pressure_derivs(double T) const
    {
        check_temperature_range(T, "get_pressure_derivs");

        double p_plus, dp_plus, ddp_plus;
        double p_minus, dp_minus, ddp_minus;

        alglib::spline1ddiff(this->p_plus_spline, T, p_plus, dp_plus, ddp_plus);
        alglib::spline1ddiff(this->p_minus_spline, T, p_minus, dp_minus, ddp_minus);

        return {{p_plus, dp_plus, ddp_plus}, {p_minus, dp_minus, ddp_minus}};
    }

    std::pair<double, double>
    EquationOfState::get_sound_speed(double T)  const
    {
        check_temperature_range(T, "get_sound_speed");
        double dp_plus = get_pressure_derivs(T).first[1];
        double de_plus = get_energy_derivs(T).first[1];
        double dp_minus = get_pressure_derivs(T).second[1];
        double de_minus = get_energy_derivs(T).second[1];
        
        double sound_speed_plus = std::sqrt(dp_plus / de_plus);
        double sound_speed_minus = std::sqrt(dp_minus / de_minus);

        return {sound_speed_plus, sound_speed_minus};
    }

    std::pair<double, double>
    EquationOfState::get_theta(double T)  const
    {
        check_temperature_range(T, "get_theta");
        double e_plus = alglib::spline1dcalc(this->e_plus_spline, T);
        double p_plus = alglib::spline1dcalc(this->p_plus_spline, T);
        double e_minus = alglib::spline1dcalc(this->e_minus_spline, T);
        double p_minus = alglib::spline1dcalc(this->p_minus_spline, T);

        double theta_plus = (e_plus - 3.0 * p_plus) / 4.0;
        double theta_minus = (e_minus - 3.0 * p_minus) / 4.0;
        return {theta_plus, theta_minus};
    }

    double
    EquationOfState::get_energy_plus(double T) const
    {
        check_temperature_range(T, "get_energy_plus");
        return alglib::spline1dcalc(this->e_plus_spline, T);
    }

    double
    EquationOfState::get_energy_minus(double T) const
    {
        check_temperature_range(T, "get_energy_minus");
        return alglib::spline1dcalc(this->e_minus_spline, T);
    }

    double
    EquationOfState::get_pressure_plus(double T) const
    {
        check_temperature_range(T, "get_pressure_plus");
        return alglib::spline1dcalc(this->p_plus_spline, T);
    }

    double
    EquationOfState::get_pressure_minus(double T) const
    {
        check_temperature_range(T, "get_pressure_minus");
        return alglib::spline1dcalc(this->p_minus_spline, T);
    }

    double
    EquationOfState::get_enthalpy_plus(double T) const
    {
        check_temperature_range(T, "get_enthalpy_plus");
        return alglib::spline1dcalc(this->w_plus_spline, T);
    }

    double
    EquationOfState::get_enthalpy_minus(double T) const
    {
        check_temperature_range(T, "get_enthalpy_minus");
        return alglib::spline1dcalc(this->w_minus_spline, T);
    }

    double
    EquationOfState::get_entropy_plus(double T) const
    {
        check_temperature_range(T, "get_entropy_plus");
        return alglib::spline1dcalc(this->s_plus_spline, T);
    }

    double
    EquationOfState::get_entropy_minus(double T) const
    {
        check_temperature_range(T, "get_entropy_minus");
        return alglib::spline1dcalc(this->s_minus_spline, T);
    }

    std::vector<double>
    EquationOfState::get_energy_derivs_plus(double T) const
    {
        check_temperature_range(T, "get_energy_derivs_plus");
        double e, de, dde;
        alglib::spline1ddiff(this->e_plus_spline, T, e, de, dde);
        return {e, de, dde};
    }

    std::vector<double>
    EquationOfState::get_energy_derivs_minus(double T) const
    {
        check_temperature_range(T, "get_energy_derivs_minus");
        double e, de, dde;
        alglib::spline1ddiff(this->e_minus_spline, T, e, de, dde);
        return {e, de, dde};
    }

    std::vector<double>
    EquationOfState::get_pressure_derivs_plus(double T) const
    {
        check_temperature_range(T, "get_pressure_derivs_plus");
        double p, dp, ddp;
        alglib::spline1ddiff(this->p_plus_spline, T, p, dp, ddp);
        return {p, dp, ddp};
    }

    std::vector<double>
    EquationOfState::get_pressure_derivs_minus(double T) const
    {
        check_temperature_range(T, "get_pressure_derivs_minus");
        double p, dp, ddp;
        alglib::spline1ddiff(this->p_minus_spline, T, p, dp, ddp);
        return {p, dp, ddp};
    }

    double
    EquationOfState::get_sound_speed_plus(double T) const
    {
        check_temperature_range(T, "get_sound_speed_plus");
        double dp = get_pressure_derivs_plus(T)[1];
        double de = get_energy_derivs_plus(T)[1];
        return std::sqrt(dp / de);
    }

    double
    EquationOfState::get_sound_speed_minus(double T) const
    {
        check_temperature_range(T, "get_sound_speed_minus");
        double dp = get_pressure_derivs_minus(T)[1];
        double de = get_energy_derivs_minus(T)[1];
        return std::sqrt(dp / de);
    }

    double
    EquationOfState::get_theta_plus(double T) const
    {
        check_temperature_range(T, "get_theta_plus");
        double e = alglib::spline1dcalc(this->e_plus_spline, T);
        double p = alglib::spline1dcalc(this->p_plus_spline, T);
        return (e - 3.0 * p) / 4.0;
    }

    double
    EquationOfState::get_theta_minus(double T) const
    {
        check_temperature_range(T, "get_theta_minus");
        double e = alglib::spline1dcalc(this->e_minus_spline, T);
        double p = alglib::spline1dcalc(this->p_minus_spline, T);
        return (e - 3.0 * p) / 4.0;
    }

    void
    EquationOfState::write(const std::string path) const 
    {
        std::ofstream file(path);

        if (!file.is_open()) { throw std::runtime_error("Could not open file for writing EoS data"); }

        file << "T,P_plus,P_minus,e_plus,e_minus,w_plus,w_minus,s_plus,s_minus\n";

        double dtt = (t_max - t_min)/(n_temp-1);
        for (double tt = t_min; tt < t_max; tt += dtt) 
        {
        file << tt << ","
            << get_pressure_plus(tt) << "," << get_pressure_minus(tt) << ","
            << get_energy_plus(tt) << "," << get_energy_minus(tt) << ","
            << get_enthalpy_plus(tt) << "," << get_enthalpy_minus(tt) << ","
            << get_entropy_plus(tt) << "," << get_entropy_minus(tt) << "\n";
        }
    }

} // namespace PhaseTracer