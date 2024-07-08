#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>
#include <random>
#include "gravwave_calculator.hpp"

namespace PhaseTracer {

double GravWaveCalculator::GW_bubble_collision(double f, double alpha, double beta_H, double T_ref) {
    double omega_env;
    double s_env;
    double kappa = 1 / (1 + 0.715 * alpha) * (0.715 * alpha + 4. / 27 * sqrt(3 * alpha / 2));
    double f_peak_beta = 0.35 / (1 + 0.069 * vw + 0.69 * pow(vw, 4.));
    double f_env = 1.65e-5 * f_peak_beta * beta_H * (T_ref / 100) * pow(dof / 100, 1./ 6);
    double delta = 0.48 * pow(vw, 3.) / (1 + 5.3 * pow(vw, 2.) + 5 * pow(vw, 4.));
    s_env = pow(0.064 * pow(f/f_env, -3.) + (1 - 0.064 - 0.48) * pow(f / f_env, -1.) + 0.48 * (f / f_env),-1);
    omega_env = 1.67e-5 * delta * pow(beta_H, -2.) * pow(kappa * alpha / (1 + alpha), 2.) * pow(100 / dof, 1 / 3.) * s_env;
    return omega_env;
}

double GravWaveCalculator::GW_sound_wave(double f, double alpha, double beta_H, double T_ref) {
    double omega_sw;
    double omega_sw_peak;
    double Hstar_R = pow(beta_H, -1.) * pow(8 * 3.1415926, 1./ 3) * vw;
    double f_peak_sw = 2.6e-5 / Hstar_R * (T_ref / 100) * pow(dof / 100, 1./ 6);
    double K_sw = Kappa_sound_wave(alpha) * alpha / (1 + alpha);
    double H_tau = Hstar_R / sqrt(3 * K_sw /4);
    omega_sw_peak = 2.061 * 0.678 * 0.678 * 3.57e-5 * pow(100 / dof, 1./ 3) * 0.012 * Hstar_R * pow(K_sw, 2.);
    omega_sw = omega_sw_peak * pow(f / f_peak_sw, 3.) *
            pow(7. / (4. + 3. * pow(f / f_peak_sw, 2.)),7. / 2);
    omega_sw = omega_sw * std::min(1.0, H_tau);
    return omega_sw;
}

double GravWaveCalculator::GW_turbulence(double f, double alpha, double beta_H, double T_ref) {
    double omega_turb;
    double hn = 1.65e-5 * (T_ref / 100) * pow(dof/100, 1./ 6);
    double f_peak_turb = 2.7e-5 / vw * beta_H * (T_ref / 100) * pow(dof / 100, 1./ 6);
    double kappa_turb = Kappa_sound_wave(alpha) * epsilon;
    omega_turb = 3.35e-4 * pow(beta_H, -1.) * pow(kappa_turb * alpha / (1 + alpha), 3./2) * pow(100 / dof, 1./3)
            * vw * pow(f / f_peak_turb, 3) / (pow(1 + f / f_peak_turb, 11./3) * (1 + 8 * 3.1415926 * f / hn));
    return omega_turb;
}

void GravWaveCalculator::Set_frequency_list(double begin_log_frequency, double end_log_frequency, double num_frequency) {
    if (num_frequency <= 1) {
        std::cerr << "Error: Number of frequencies must be greater than 1." << std::endl;
        std::exit(1);
    }

    if (end_log_frequency < begin_log_frequency){
        std::cerr << "Error: The end log frequency must be larger than begin log frequency" << std::endl;
        std::exit(1);
    }


    double step = (end_log_frequency - begin_log_frequency) / (num_frequency - 1);
    for (int i = 0; i < num_frequency; ++i) {
        double frequency = begin_log_frequency + i * step;
        double value = std::pow(10, frequency);
        frequency_list.push_back(value);
    }

}

void GravWaveCalculator::Print_parameter() {
    std::string dashes(10, '-');
    for (std::size_t i=0; i < transitions_params.size(); ++i){
    	std::cout << "The GW Parameters of the "<< i + 1 << " th transition is: " << std::endl;
    	std::cout << "alpha is: " << transitions_params[i].alpha <<std::endl;
    	std::cout << "beta over H is: " << transitions_params[i].beta_H <<std::endl;
    	std::cout << "vw is: " << vw <<std::endl;
    	std::cout << "T_ref is: " << transitions_params[i].T_ref <<std::endl;
    	std::cout << "epsilon is: " << epsilon <<std::endl;
    	std::cout << "Degree of freedom is: " << dof <<std::endl;
    	std::cout << dashes <<std::endl;
    	std::cout << std::endl;
    }
}

double GravWaveCalculator::Kappa_sound_wave(double alpha) {
    double kappa_sw;
    double cs = sqrt(1 / 3.);
    double v_cj = 1 / (1 + alpha) * (cs + sqrt(pow(alpha,2.) + 2./3 * alpha));
    double kappa_a = pow(vw, 6./5) * 6.9 * alpha / (1.36 - 0.037 * sqrt(alpha) + alpha);
    double kappa_b = pow(alpha, 2./5) / (0.017 + pow(0.997 + alpha, 2./5));
    double kappa_c = sqrt(alpha) / (0.135 + sqrt(0.98 + alpha ));
    double kappa_d = alpha / (0.73 + 0.083 * sqrt(alpha) + alpha);
    double delta_kappa = -0.9 * log(sqrt(alpha) / (1 + sqrt(alpha)));

    if (0. < vw && vw <= cs)
    {
        kappa_sw = pow(cs, 11./ 5) * kappa_a * kappa_b / ((pow(cs, 11./ 5) - pow(vw, 11./ 5)) * kappa_b + vw * pow(cs, 6./ 5) * kappa_a);
    }
    else if (cs < vw && vw < v_cj)
    {
        kappa_sw = kappa_b + (vw - cs) * delta_kappa + pow(vw-cs, 3.) / pow(v_cj-cs, 3.) * (kappa_c - kappa_b - (v_cj -cs) * delta_kappa);
    }
    else if (v_cj <= vw && vw <= 1.)
    {
        kappa_sw = pow(v_cj-1, 3.) * pow(v_cj, 5./ 2) * pow(vw, -5./ 2) * kappa_c * kappa_d / ((pow(v_cj-1, 3.) - pow(vw-1,3.)) *
                pow(v_cj, 5./ 2) * kappa_c + pow(vw-1,3.) * kappa_d);
    }
    else
    {
        std::cerr<<"Error: Wrong bubble velocity, please check its values !";
        std::exit(1);
    }
    return kappa_sw;

}

void GravWaveCalculator::GW_total_spectrum() {
    std::vector<double> omega_tot;
    for (std::size_t i=0; i < transitions_params.size(); ++i){
    	double omega_env = 0;
    	double omega_sw;
    	double omega_turb;
    	double tot_value;
    	double peak_frequency = 0;
    	double peak_spectrum = 0;	
    	double alpha = transitions_params[i].alpha;
    	double beta_H = transitions_params[i].beta_H;
    	double T_ref = transitions_params[i].T_ref;
    	
    	for (double value : frequency_list){
        	if (T_ref < 10){
            	omega_env = GW_bubble_collision(value, alpha, beta_H, T_ref);
        	}        	
        	omega_sw = GW_sound_wave(value, alpha, beta_H, T_ref);
        	omega_turb = GW_turbulence(value, alpha, beta_H, T_ref);
        	tot_value = omega_env + omega_sw + omega_turb;
        	omega_tot.push_back(tot_value);
        	
        	if (tot_value > peak_spectrum){
            	peak_spectrum = tot_value;
            	peak_frequency = value;}
      }
      std::stringstream filename_string;
      filename_string<<"GW_results_of_transition_" << i << ".csv";
      std::string filename = filename_string.str();
      GravWaveCalculator::Write_to_csv(std::make_tuple(peak_frequency, peak_spectrum, frequency_list, omega_tot), filename);
      omega_tot.clear(); //reset vector
}	
}

void GravWaveCalculator::Write_to_csv(const std::tuple<double, double, std::vector<double>, std::vector<double>> &data, const std::string &filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file: " << filename << std::endl;
        return;
    }
    
    file << "Peak frequency, Peak spectrum, Frequency list, Spectrum list" << std::endl;
    
    double peak_frequency = std::get<0>(data);
    double peak_spectrum = std::get<1>(data);
    std::vector<double> f_list = std::get<2>(data);
    std::vector<double> omega_list = std::get<3>(data);

    size_t max_size = f_list.size();
    
    for (size_t i = 0; i < max_size; ++i) {
        if (i == 0) {
            file << peak_frequency << ",";
            file << peak_spectrum << ",";
        } else {
            file << ",";
            file << ",";
        }

        if (i < f_list.size()) {
            file << f_list[i];
        }
        file << ",";

        if (i < omega_list.size()) {
            file << omega_list[i];
        }
        file << std::endl;
    }

    std::cout << "GW results has been written to " << filename << std::endl;
    file.close();
}


}
