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

#ifndef PHASETRACER_POLYNOMIAL_FITTER_HPP_
#define PHASETRACER_POLYNOMIAL_FITTER_HPP_

#include <algorithm>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>

#include <eigen3/Eigen/Core>
#include <boost/cstdint.hpp>

namespace PhaseTracer {

class PolynomialFitterEigen {
public:
  void fit(const std::vector<double> &T_list_, const std::vector<double> &S_list_, double TC_, int degree) {
    TC = TC_;
    T_list = T_list_;
    S_list = S_list_;
    fit_flag = true;

    // Filter the nodes
    std::vector<double> T_select;
    std::vector<double> S_select;
    T_select.reserve(T_list.size());
    S_select.reserve(T_list.size());
    for (size_t i = 0; i < T_list.size(); ++i) {
      double S3T = S_list[i] / T_list[i];
      if (S3T > 1 and S3T < 1000) {
        T_select.push_back(T_list[i]);
        S_select.push_back(S_list[i]);
      }
    }
    T_select.shrink_to_fit();
    S_select.shrink_to_fit();

    // TODO: this is only for test
    MSE_pre = fit_(T_select, S_select, degree);
    LOG(debug) << "MSE of action fit (pre selection) = " << MSE_pre;

    select(T_select, S_select);
    if (T_select.size() < degree * 2)
      return;

    MSE = fit_(T_select, S_select, degree);
    LOG(debug) << "MSE of action fit = " << MSE;
    //    if (MSE < 10)
    success = true;
  }

  double fit_(const std::vector<double> T_select, const std::vector<double> S_select, int degree) {

    int n = T_select.size();
    std::vector<double> x = T_select;
    std::vector<double> y;
    y.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      y.push_back(S_select[i] * pow(TC - T_select[i], 2));
    }
    y.shrink_to_fit();

    Eigen::MatrixXd X(n, degree + 1);
    Eigen::VectorXd Y(n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= degree; j++) {
        X(i, j) = std::pow(x[i], j);
      }
      Y(i) = y[i];
    }
    coefficients = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

    fit_range_min = *std::min_element(T_select.begin(), T_select.end());
    fit_range_max = *std::max_element(T_select.begin(), T_select.end());

    // Calculate mean square error
    double MSE_ = 0;
    S_predict.clear();
    S_predict.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      S_predict.push_back(predict(T_select[i]));
      MSE_ += pow(S_select[i] / T_select[i] - S_predict[i] / pow(TC - T_select[i], 2) / T_select[i], 2);
    }
    S_predict.shrink_to_fit();
    MSE_ = MSE_ / n;
    return MSE_;
  }

  void select(std::vector<double> &x, std::vector<double> &y) {

    if (x.size() != y.size() || x.size() < 3) {
      throw std::runtime_error("Need at least 3 points for action fit");
    }
    if (x[0] < x.back()) {
      std::reverse(x.begin(), x.end());
      std::reverse(y.begin(), y.end());
    }
    std::vector<int> selected_indices;
    for (size_t i = 0; i < y.size() - 1; ++i) {
      if (y[i + 1] - y[i] < 0) {
        selected_indices.push_back(i);
        selected_indices.push_back(i + 1);
        break;
      }
    }
    if (selected_indices.empty()) {
      return;
    }

    double y_prev = y[selected_indices.back()];
    double diff_prev = y[selected_indices[1]] - y[selected_indices[0]];
    for (size_t i = selected_indices.back() + 1; i < y.size(); ++i) {
      double diff = y[i] - y_prev;
      if ((diff < 0) && (std::abs(diff) < std::abs(diff_prev * 2))) {
        selected_indices.push_back(i);
        diff_prev = diff;
        y_prev = y[i];
      }
    }
    std::vector<double> x_selected;
    std::vector<double> y_selected;
    for (int idx : selected_indices) {
      x_selected.push_back(x[idx]);
      y_selected.push_back(y[idx]);
    }
    x = std::move(x_selected);
    y = std::move(y_selected);
  }

  double predict(double x) const {
    if (x < fit_range_min or x > fit_range_max) {
      LOG(warning) << "Input value exceeds the action fitting function's domain.";
    }
    double result = 0.0;
    for (int i = 0; i < coefficients.size(); i++) {
      result += coefficients[i] * std::pow(x, i);
    }
    return result;
  }
  std::vector<double> predict(const std::vector<double> &x) const {
    std::vector<double> result;
    result.reserve(x.size());
    for (int i = 0; i < x.size(); i++) {
      result.push_back(predict(x[i]));
    }
    return result;
  }

  double derivative(double x) const {
    double result = 0.0;
    for (int i = 1; i < coefficients.size(); i++) {
      result += i * coefficients[i] * std::pow(x, i - 1);
    }
    return result;
  }
  double secondDerivative(double x) const {
    double result = 0.0;
    for (int i = 2; i < coefficients.size(); i++) {
      result += i * (i - 1) * coefficients[i] * std::pow(x, i - 2);
    }
    return result;
  }

  double findLocalMinimum(double initial_guess = 0.0, double tolerance = 1e-4, int max_iterations = 100) const {
    double x = initial_guess;
    for (int i = 0; i < max_iterations; i++) {
      double f = derivative(x);
      double f_prime = secondDerivative(x);
      if (std::abs(f_prime) < 1e-12 || f_prime < 0) {
        x = x - 0.1 * f;
        continue;
      }
      double delta = f / f_prime;
      x = x - delta;
      if (std::abs(delta) < tolerance && f_prime > 0) {
        break;
      }
    }
    if (secondDerivative(x) > 0) {
      return x;
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  const Eigen::VectorXd &getCoefficients() const {
    return coefficients;
  }

  const double get_MSE() const {
    return MSE;
  }

  const double get_MSE_pre() const {
    return MSE_pre;
  }

  const double get_S3T(double x) const {
    double S3T = predict(x) / pow(TC - x, 2) / x;
    return std::min(S3T, 1E30);
  }

  const double get_fit_flag() const {
    return fit_flag;
  }

  const bool get_success() const {
    return success;
  }

  const std::vector<double> get_T_list() const {
    return T_list;
  }

  const std::vector<double> get_S_list() const {
    return S_list;
  }

private:
  Eigen::VectorXd coefficients;
  double TC;
  std::vector<double> T_list;
  std::vector<double> S_list;
  std::vector<double> S_predict;
  double fit_range_min;
  double fit_range_max;
  double MSE = 1E10;
  double MSE_pre = 1E10; // TODO: this is only for test
  bool success = false;
  bool fit_flag = false;
};

} // namespace PhaseTracer

#endif // PHASETRACER_POLYNOMIAL_FITTER_HPP_
