#ifndef PHASE_TRACER_GW_GRAVWAVECALCULATOR_H
#define PHASE_TRACER_GW_GRAVWAVECALCULATOR_H
#include <vector>
#include <tuple>
#include <fstream>

#include "transition_finder.hpp"


namespace PhaseTracer {


struct GravWaveTransition {
  double alpha;
  double beta_H;
  double T_ref;
};

class GravWaveCalculator {
private:
    double alpha;
    double beta_H;
    double vw;
    double T_ref;
    std::vector<double> frequency_list;
  PROPERTY(double, dof, 106.75)
  PROPERTY(double, epsilon, 0.1)
  PROPERTY(double, begin_log_frequency, -4)
  PROPERTY(double, end_log_frequency, 1)
  PROPERTY(double, num_frequency, 500)
  
  TransitionFinder tf;
  std::vector<GravWaveTransition> transitions;
  
public:
  explicit GravWaveCalculator(TransitionFinder tf_) :
  tf(tf_) {
    std::vector<Transition> trans = tf.get_transitions();
    for (const auto& ii : trans){
      double action = tf.get_action(ii.true_phase,ii.false_phase,ii.TN,ii.key);
      LOG(fatal) << "action(TN) = "<< action;
      action = tf.get_action(ii.true_phase,ii.false_phase,ii.TN-0.1,ii.key);
      LOG(fatal) << "action(TN-0.1) = "<< action;
      std::vector<Eigen::VectorXd>  vacua = tf.get_vacua_at_T(ii.true_phase,ii.false_phase,ii.TN,ii.key);
      LOG(fatal) << "vacuum 1 = "<< vacua[0];
//      double V1 =
      LOG(fatal) << "V(vacuum 1, TN) = "<< tf.pf.P.V(vacua[0],ii.TN);
      
      LOG(fatal) << "vacuum 2 = "<< vacua[1];
//      double V1 =
      LOG(fatal) << "V(vacuum 2, TN) = "<< tf.pf.P.V(vacua[1],ii.TN);
      
      transitions.push_back({0,0,ii.TN});
    }
    }
  virtual ~GravWaveCalculator() = default;
    
  
    void Set_parameters(double alpha_input, double beta_H_input, double vw_input, double T_ref_input);
    void Print_parameter();
    void Write_to_csv(const std::tuple<double, double, std::vector<double>, std::vector<double>>& data, const std::string& filename="GW_results.csv");
    double Kappa_sound_wave();
    double GW_bubble_collision(double f);
    double GW_sound_wave(double f);
    double GW_turbulence(double f);
    std::tuple<double, double, std::vector<double>,std::vector<double>> GW_total_spectrum();

};

}  // namespace PhaseTracer
#endif //PHASE_TRACER_GW_GRAVWAVECALCULATOR_H
