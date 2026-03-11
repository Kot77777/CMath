#include "gas_solver.hpp"

#include <cmath>
#include "State.hpp"
#include "matrix.hpp"

std::vector<State> gas_solver(const double L, const std::size_t NX,
                              const double tau_init,
                              const std::vector<State> &state_0,
                              const double t_end, const double gamma, const double CFL_max) {
    std::vector<State> state_prev{state_0};
    const double h = 2 * L / (NX - 1);
    double t = 0, tau = tau_init, lambda_max = 0;
    while (t < t_end) {
        std::vector<State> state_next(NX);
        for (std::size_t i = 1; i < NX - 1; ++i) {
            const Eigen::Vector3d w_nl = state_prev[i].w_;
            const Eigen::Vector3d w_nl_prev = state_prev[i - 1].w_;
            const Eigen::Vector3d w_nl_next = state_prev[i + 1].w_;
            const matrix m{state_prev[i], gamma};

            const Eigen::Vector3d w_n_next_l = w_nl - tau * m.A_ * (w_nl_next - w_nl_prev) / (2 * h) + tau * (
                                                   m.omega_inv * m.lambda_abs * m.omega_) * (
                                                   w_nl_next - 2 * w_nl + w_nl_prev) / (2 * h);
            const State s{w_n_next_l, gamma};
            lambda_max = std::max(lambda_max, std::abs(s.u_ + s.c_));
            state_next[i] = s;
        }
        state_next[0] = state_next[1];
        state_next[NX - 1] = state_next[NX - 2];
        state_prev = state_next;
        if (tau * lambda_max / h > CFL_max) {
            tau = 0.95 * CFL_max * h / lambda_max;
        }
        lambda_max = 0;
        t += tau;
    }
    return state_prev;
}
