#include <fstream>

#include "cmath/PDE/gas_equation/gas_solver.hpp"
#include "cmath/PDE/gas_equation/State.hpp"
#include "gtest/gtest.h"

TEST(gas_equation, test) {
    constexpr std::size_t NX = 100;
    constexpr std::size_t NX_2 = NX / 2;
    constexpr double L = 10;
    constexpr double tau_init = 1e-7;
    constexpr double t_end = 0.015;
    constexpr double gamma = 5. / 3;
    constexpr double u_R = 0, u_L = 0;
    constexpr double rho_R = 1.3;
    constexpr double rho_L = 13;
    constexpr double P_R = 1 * 1e5;
    constexpr double P_L = 10 * 1e5;
    constexpr double e_R = P_R / rho_R / (gamma - 1);
    constexpr double e_L = P_L / rho_L / (gamma - 1);

    const State state_R{rho_R, u_R, e_R, gamma};
    const State state_L{rho_L, u_L, e_L, gamma};

    std::vector<State> state0(NX);
    for (std::size_t i = 0; i < NX_2; ++i) {
        state0[i] = state_L;
    }
    for (std::size_t i = NX_2; i < NX; ++i) {
        state0[i] = state_R;
    }

    const auto res = gas_solver(L, NX, tau_init, state0, t_end, gamma, 0.01);

    std::ofstream res_gas("gas_equation.csv");
    res_gas << "rho" << std::endl;
    for (const auto &i: res) {
        res_gas << i.rho_ << std::endl;
    }
}