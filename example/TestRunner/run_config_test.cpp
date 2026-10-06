#include <iostream>

#include "phasetracer.hpp"

int main(int argc, char *argv[]) {


    auto config = PhaseTracer::Config();
    config.phase_finder.t_low = -1.0;

    auto status = config.validate();
    std::cout << status;

    return 0;
}
