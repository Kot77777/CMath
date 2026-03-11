#pragma once
#include <vector>
#include "State.hpp"

std::vector<State> gas_solver(double L, std::size_t NX, double tau_init,
                              const std::vector<State> &state_0,
                              double t_end, double gamma, double CFL_max);
