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

