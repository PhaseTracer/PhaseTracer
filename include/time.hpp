#include <chrono>
#include <iostream>

#define START_TIMER const int n_repeats = 1000; \
  auto t0 = std::chrono::high_resolution_clock::now(); \
  for (int i = 0; i < n_repeats; i++) {

#define STOP_TIMER } \
  auto t1 = std::chrono::high_resolution_clock::now(); \
  auto dt_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count(); \
  std::cout << "seconds per repeat = " << 1.e-6 * static_cast<double>(dt_microseconds) / n_repeats << std::endl;
