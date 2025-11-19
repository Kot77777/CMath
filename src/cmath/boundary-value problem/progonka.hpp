#pragma once
#include <vector>
#include <solution_SLAE/method_progonky.h>

template<typename Callable>
[[nodiscard]] std::tuple<std::vector<double>, std::vector<double> > progonka(const Callable &k, const auto &state_begin,
                                                                             const auto &state_end,
                                                                             const std::size_t N) {
    const auto [x_begin, y_begin] = state_begin;
    const auto [x_end, y_end] = state_end;
    const double h = (x_end - x_begin) / (N - 1);
    double x = x_begin;
    std::vector<double> alpha(N - 1), betta(N), gamma(N - 1), delta(N), x_vec(N);
    x_vec[0] = x;
    betta[0] = 1.;
    gamma[0] = 0.;
    delta[0] = y_begin;
    for (std::size_t i = 1; i < N - 1; ++i) {
        x += h;
        alpha[i - 1] = 1. - k.a(x) * h / 2;
        betta[i] = -2. + k.b(x) * h * h;
        gamma[i] = 1. + k.a(x) * h / 2;
        delta[i] = k.c(x) * h * h;
        x_vec[i] = x;
    }
    betta[N - 1] = 1.;
    alpha[N - 2] = 0.;
    delta[N - 1] = y_end;
    x_vec[N - 1] = x_end;

    const auto res = method_progonky(alpha, betta, gamma, delta);
    return {res, x_vec};
}
