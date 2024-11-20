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

#ifndef PHASETRACER_TIME_HPP_
#define PHASETRACER_TIME_HPP_

#include <chrono>
#include <iostream>

#define START_TIMER const int n_repeats = 1000; \
  auto t0 = std::chrono::high_resolution_clock::now(); \
  for (int i = 0; i < n_repeats; i++) {

#define STOP_TIMER } \
  auto t1 = std::chrono::high_resolution_clock::now(); \
  auto dt_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count(); \
  std::cout << "seconds per repeat = " << 1.e-6 * static_cast<double>(dt_microseconds) / n_repeats << std::endl;

#endif  // PHASETRACER_TIME_HPP_