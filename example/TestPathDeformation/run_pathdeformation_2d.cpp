/**
 The scale example BubbleProfiler
 ./run_BP_scale 1. 0.6 200.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "shooting.hpp"
#include "path_deformation.hpp"

class Example: public EffectivePotential::Potential {
public:
  Example(double c_, double fx_, double fy_) :
    c(c_), fx(fx_), fy(fy_) {}
  
  size_t get_n_scalars() const override {return 2;}
  
  double V(Eigen::VectorXd phi, double T) const override {
    double x=phi[0];
    double y=phi[1];
    double r1 = x*x+c*y*y;
    double r2 = c*pow(x-1.,2) + pow(y-1.,2);
    double r3 = fx*(0.25*pow(x,4) - pow(x,3)/3.);
    r3 += fy*(0.25*pow(y,4) - pow(y,3)/3.);
    return r1*r2 + r3;
  }
  
  virtual Eigen::VectorXd dV_dx(Eigen::VectorXd phi, double T) const override {
    double x=phi[0];
    double y=phi[1];
    double r1 = x*x+c*y*y;
    double r2 = c*pow(x-1.,2) + pow(y-1.,2);
    double dr1dx = 2*x;
    double dr1dy = 2*c*y;
    double dr2dx = 2*c*(x-1.);
    double dr2dy = 2*(y-1.);
    double dVdx = r1*dr2dx + dr1dx*r2 + fx*x*x*(x-1.);
    double dVdy = r1*dr2dy + dr1dy*r2 + fy*y*y*(y-1.);
    Eigen::VectorXd rval(2);
    rval << dVdx, dVdy;
    return rval;
  }
  
private:
  double c, fx, fy;
  
};

int main(int argc, char* argv[]) {

//  LOGGER(fatal);
  LOGGER(debug);
  
  Example p(5,0.,10);
  
  PhaseTracer::PathDeformation pd(p);
  
  std::vector<Eigen::VectorXd> path_pts;
  path_pts.push_back(Eigen::VectorXd(2));
  path_pts.push_back(Eigen::VectorXd(2));
  path_pts[0] << 1, 1;
  path_pts[1] << 0, 0;
  auto a = pd.full_tunneling(path_pts);

  LOG(debug)<< "Action = "<< std::setprecision(10) << a.action;
  LOG(debug)<< "fRatio = "<< std::setprecision(10) << a.fRatio;
  
  std::ofstream file1("profile1D_for_2d_example.txt");
  file1 << "R,Phi,dPhi\n";
  for (int i = 0; i < a.profile1D.R.size(); ++i) {
    file1 << std::setprecision(10) << a.profile1D.R[i] << "," << a.profile1D.Phi[i] << "," << a.profile1D.dPhi[i] << std::endl;
  }
  file1.close();

  std::ofstream file2("path_for_2d_example.txt");
  for (int i = 0; i < a.phi.size(); ++i) {
    file2 << std::setprecision(10) << a.phi[i][0] << "\t" << a.phi[i][1]  << std::endl;
  }
  file2.close();
  
  return 0;
}
