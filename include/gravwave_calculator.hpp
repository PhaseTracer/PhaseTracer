#ifndef PHASETRACER_GRAVWAVECALCULATOR_H
#define PHASETRACER_GRAVWAVECALCULATOR_H
#include <vector>
#include <tuple>
#include <fstream>
#include <cmath>
#include "transition_finder.hpp"
#include <random>


namespace PhaseTracer {


struct GravWavetrans_params {
  double alpha;
  double beta_H;
  double T_ref;
};

class GravWaveCalculator {
private:
    std::vector<double> frequency_list;
    std::vector<GravWavetrans_params> transitions_params;
    PROPERTY(int, num_random_iterations, 5)
    PROPERTY(double, dof, 106.75)
    PROPERTY(double, epsilon, 0.1)
    /** The step-size in numerical derivatives */
    PROPERTY(double, h, 1e-2)
    /** The number of points used in numerical derivative of action */
    PROPERTY(double, np_dSdT, 5)
  
    PROPERTY(double, vw, 0.3)
    const double G = 6.7088e-39;
    TransitionFinder tf;
    
    const double rho_R(double T){
      return std::pow(T, 4) * M_PI * M_PI * dof / 30.;
    }

    const double H(double T){
      return std::sqrt(8. * M_PI * G / 3. * rho_R(T));
    }

    /* Calcualte alpha with fixed phi */
    const double V(const Eigen::VectorXd phi, const double T) {
      return tf.pf.P.V(phi,T);
    }
    const double dVdT(const Eigen::VectorXd phi, const double T) {
      return (-V(phi,T+2*h) + 8*V(phi,T+h) - 8*V(phi,T-h) + V(phi,T-2*h)) / (12.0 * h);  
    }  
    const double rho(const Eigen::VectorXd phi, const double T){
      return V(phi,T) - 0.25*T*dVdT(phi,T);
    }
    const double alpha(const Eigen::VectorXd vacuum_1, const Eigen::VectorXd vacuum_2, double T){
      return ( rho(vacuum_1,T) - rho(vacuum_2,T) ) / rho_R(T);
    }
    
    /* Calcualte alpha along the phase */
    Eigen::VectorXd vacuum_of_phase_at_T(Phase phase1, double T){
      return tf.pf.phase_at_T(phase1, T).x;
    }
    const double V(Phase phase1, const double T) {
      return tf.pf.P.V(vacuum_of_phase_at_T(phase1,T),T);
    }
    const double dVdT(Phase phase1, const double T) {
      return (-V(phase1,T+2*h) + 8*V(phase1,T+h) - 8*V(phase1,T-h) + V(phase1,T-2*h)) / (12.0 * h);  
    }  
    const double rho(Phase phase1, const double T){
      return V(phase1,T) - 0.25*T*dVdT(phase1,T);
    }
    const double alpha(Phase phase1, Phase phase2, double T){
      return ( rho(phase1,T) - rho(phase1,T) ) / rho_R(T);
    }

    /* Using least squares method to fit a quadratic polynomial for solving the derivative of the action functional. */
    Eigen::Vector3d quadraticFit(Eigen::VectorXd& x, Eigen::VectorXd& y) {  
      int n = x.size();  
      Eigen::Matrix3d A;  
      Eigen::Vector3d b;  
      A << n, x.sum(), (x.array() * x.array()).sum(),  
         x.sum(), (x.array() * x.array()).sum(), (x.array() * x.array() * x.array()).sum(),  
         (x.array() * x.array()).sum(), (x.array() * x.array() * x.array()).sum(), (x.array() * x.array() * x.array() * x.array()).sum();  
      b << y.sum(),  
         (x.array() * y.array()).sum(),  
         (x.array() * x.array() * y.array()).sum();  
      // slove A * coeff = b  
      Eigen::Vector3d coeff = A.colPivHouseholderQr().solve(b);  
      return coeff;  
    }
  
    const double S3T(Phase phase1, Phase phase2, double T, size_t i_unique){
      return tf.get_action(phase1, phase2, T, i_unique)/T;
    }
    
    const double dSdT(Phase phase1, Phase phase2, double T, size_t i_unique){
      Eigen::VectorXd x,y;
      for (int ii=0; ii <= np_dSdT; ii++ ){
        double Ti = T - h + ii * 2. * h / np_dSdT;
        double S = S3T(phase1,phase2,T,i_unique);
        if ( (not std::isnan(S)) and (not std::isnan(S)) ) { 
          x << Ti;
          y << S3T(phase1,phase2,T,i_unique);
        }
      }
      // TODO add check on num of valid x
      Eigen::Vector3d coeff = quadraticFit(x, y); // coeff[0] x^2 + coeff[1] x + coeff[2]
      return 2 * coeff[0] * T + coeff[1];
    }
    
    const double beta(Phase phase1, Phase phase2, double T, size_t i_unique){
      return T * H(T) * dSdT(phase1, phase2, T, i_unique);
    }
  
public:
  explicit GravWaveCalculator(TransitionFinder tf_) :
  tf(tf_) {
    std::vector<Transition> trans = tf.get_transitions();
    
    for (const auto& ti : trans){
      double TN = ti.TN;
      if (std::isnan(TN))
      {
        std::cout<<  "Nucleation temperature dose not exist. GW will not be calculated !";
      } else {
        std::vector<Eigen::VectorXd>  vacua = tf.get_vacua_at_T(ti.true_phase,ti.false_phase,TN,ti.key);

        transitions_params.push_back({alpha(vacua[0],vacua[1],TN), beta(ti.true_phase,ti.false_phase,TN,ti.key), ti.TN});     
      }
    }
  }
  virtual ~GravWaveCalculator() = default;
    
    void Set_frequency_list(double begin_log_frequency = -4, double end_log_frequency = 1, double num_frequency = 500);
    void Print_parameter();
    void Write_to_csv(const std::tuple<double, double, std::vector<double>, std::vector<double>>& data, const std::string& filename);
    double Kappa_sound_wave(double alpha);
    double GW_bubble_collision(double f, double alpha, double beta_H, double T_ref);
    double GW_sound_wave(double f, double alpha, double beta_H, double T_ref);
    double GW_turbulence(double f, double alpha, double beta_H, double T_ref);
    void GW_total_spectrum();

};

}  // namespace PhaseTracer
#endif //PHASETRACER_GRAVWAVECALCULATOR_H
