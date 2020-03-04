#include "models/Z2_scalar_singlet_model.hpp"
#include <iostream>

int main() {
    
    EffectivePotential::Z2ScalarSingletModel model;

    Eigen::VectorXd x(2);
    x << 100., 100.;
    double T = 100.;
    
    std::cout << "x = " << x << std::endl;
    std::cout << "T = " << T << std::endl;

    std::cout << "***********************************************************************" << std::endl;

    std::cout << "Total one-loop potential                = " << model.V(x,T) << std::endl;
    
    std::cout << "***********************************************************************" << std::endl;
    
    std::cout << "The derivative of the gradient of potential with respect to T = " << std::endl
              << model.d2V_dxdt(x, T) << std::endl;
    std::cout << "The Hessian matrix of potential = " << std::endl
              << model.d2V_dx2(x, T) << std::endl;

    std::cout << "***********************************************************************" << std::endl;
   
    return 0;
}

