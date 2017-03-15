/// Causal Dynamical Triangulations in C++ using CGAL
///
/// Copyright © 2016-2017 Adam Getchell
///
/// Full run-through with default options used to calculate
/// optimal values for thermalization, etc. A simpler version
/// that encompasses the entire lifecycle. Also suitable for
/// scripting parallel runs.
///
/// \done Invoke Metropolis algorithm
/// \todo Print out graph of time-value vs. volume vs. pass number

/// @file cdt-opt.cpp
/// @brief Outputs values to determine optimizations
/// @author Adam Getchell

#include <utility>

#include "Measurements.h"
#include "Metropolis.h"
#include "Simulation.h"

int main() {
  std::cout << "cdt-opt started at " << currentDateTime() << std::endl;
  constexpr auto simplices  = static_cast<uintmax_t>(64000);
  constexpr auto timeslices = static_cast<uintmax_t>(16);
  constexpr auto alpha      = static_cast<long double>(0.6);
  constexpr auto k          = static_cast<long double>(1.1);
  constexpr auto lambda     = static_cast<long double>(0.01);
  constexpr auto passes     = static_cast<std::uintmax_t>(100);
  constexpr auto checkpoint = static_cast<std::uintmax_t>(10);

  // Initialize simulation
  Simulation my_simulation;

  // Initialize the Metropolis algorithm
  Metropolis my_algorithm(alpha, k, lambda, passes, checkpoint);

  // Make a triangulation
  SimplicialManifold universe(simplices, timeslices);

  // Queue up simulation with desired algorithm
  my_simulation.queue(
      [&my_algorithm](SimplicialManifold s) { return my_algorithm(s); });
  // Measure results
  my_simulation.queue(
      [](SimplicialManifold s) { return VolumePerTimeslice(s); });
  // my_simulation.queue(print_results())

  // Run it
  universe = my_simulation.start(std::move(universe));

  auto max_timevalue = universe.geometry->max_timevalue().get();

  /// \todo  Fix SimplicialManifold move and copy operators to save optionals
  if (max_timevalue < timeslices)
    std::cout << "You wanted " << timeslices << " timeslices, but only got "
              << max_timevalue << " ." << std::endl;
  return 0;
}
