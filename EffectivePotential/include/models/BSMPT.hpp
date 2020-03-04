#ifndef POTENTIAL_BSMPT_HPP_INCLUDED
#define POTENTIAL_BSMPT_HPP_INCLUDED

/**
  Interface to a BSMPT model.
*/

#include <vector>
#include "potential.hpp"

namespace EffectivePotential {

template<typename BSMPTModel>
class BSMPTPotential : public Potential {
 public:
  BSMPTPotential(std::string parameters) {
    model.initModel(parameters);
    model_name = typeid(BSMPTModel).name();
  }
  double V(Eigen::VectorXd phi, double T) const override {
    auto model_copy = model;
    const std::vector<double> vMin(phi.data(), phi.data() + phi.rows() * phi.cols());
    std::vector<double> vIn;
    model_copy.MinimizeOrderVEV(vMin, vIn);
    return model_copy.VEff(vIn, std::max(0., T), 0);
  }

  size_t get_n_scalars() const override { return model.nVEV; }

  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    auto phi1 = phi;
    phi1[0] = - phi[0];
    phi1[1] = - phi[1];
    phi1[2] = - phi[2];
    phi1[3] = - phi[3];
    if (model_name.find("N2HDM") != std::string::npos) {
      auto phi2 = phi;
      phi2[4] = - phi[4];    
      return {phi1,phi2};
    }
    return {phi1};
  };

 private:
  BSMPTModel model;
  std::string model_name;
};

}  // namespace EffectivePotential

#endif
