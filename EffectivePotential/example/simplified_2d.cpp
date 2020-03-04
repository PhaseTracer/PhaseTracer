#include "models/2D_test_model.hpp"
#include <iostream>

int main() {
    
    EffectivePotential::TwoDimModel model;

    Eigen::VectorXd x(2);
    x << 100., 100.;
    double T = 100.;
    
    std::cout << "x = " << x << std::endl;
    std::cout << "T = " << T << std::endl;

    std::cout << "***********************************************************************" << std::endl;

    std::cout << "Tree level potential                    = " << model.V0(x) << std::endl;
    std::cout << "Daisy term                              = " << model.daisy(x, T) << std::endl;
    std::cout << "One-loop corrections                    = " << model.V1T(x, T) << std::endl;
    std::cout << "Finite-temperature one-loop corrections = " << model.V1(x) << std::endl;
    std::cout << "Total one-loop potential                = " << model.V(x,T) << std::endl;
    
    std::cout << "***********************************************************************" << std::endl;
    
    std::cout << "The derivative of the gradient of potential with respect to T = " << std::endl
              << model.d2V_dxdt(x, T) << std::endl;
    std::cout << "The Hessian matrix of potential = " << std::endl
              << model.d2V_dx2(x, T) << std::endl;

    std::cout << "***********************************************************************" << std::endl;
   
    return 0;
}

