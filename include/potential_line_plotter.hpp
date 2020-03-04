#ifndef POTENTIAL_LINE_PLOTTER_INCLUDED
#define POTENTIAL_LINE_PLOTTER_INCLUDED

#include <boost/filesystem.hpp>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "logger.hpp"
#include "potential.hpp"

namespace PhaseTracer {

void potential_line_plotter(EffectivePotential::Potential &P, double T, Eigen::VectorXd false_vacuum, Eigen::VectorXd true_vacuum, std::string prefix = "model") {

  std::ofstream output_file;
  const auto dat_name = prefix + "_potential_line.dat";
  output_file.open(dat_name);
  output_file << std::fixed << std::setprecision(20);

  const auto delta = true_vacuum - false_vacuum;
  const size_t n = 10000;

  for (unsigned int i = 0; i <= n; i++) {
    double frac = 1.5 * i / static_cast<double>(n) - 0.25;  // from -0.25 to 1.25
    output_file << frac << " " << P(false_vacuum + frac * delta, T) << std::endl;
  }

  output_file.close();


  const boost::filesystem::path this_file(__FILE__);
  const auto this_dir = this_file.parent_path();
  const boost::filesystem::path file("../make_plots/potential_line_plotter.gnu");
  const auto program = this_dir / file;

  const auto png_name = prefix + "_potential_line.png";
  std::string vars = "png_name = '" + png_name + "'; dat_name = '" + dat_name + "'; T = " + std::to_string(T) + ";";
  std::string command = "gnuplot -e \"" + vars + "\" " + program.string();

  LOG(debug) << "Executing " << command;
  const auto result = std::system(command.c_str());
}

void potential_line_plotter(EffectivePotential::Potential &P, Transition t, std::string prefix = "model") {
  potential_line_plotter(P, t.TC, t.false_vacuum, t.true_vacuum, prefix);
}

void potential_line_plotter(EffectivePotential::Potential &P, std::vector<Transition> t, std::string prefix = "model") {
    for (size_t i = 0; i < t.size(); ++i) {
        potential_line_plotter(P, t[i], prefix + "_" + std::to_string(i));
    }
}

}  // namespace PhaseTracer

#endif
