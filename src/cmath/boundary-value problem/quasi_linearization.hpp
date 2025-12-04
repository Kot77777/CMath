#pragma once
#include <vector>
#include <cmath>
#include <solution_SLAE/method_progonky.h>

template<typename Callable>
[[nodiscard]] std::vector<double> progonka_for_ql(const Callable &k, const std::vector<double> &t,
                                                  const std::vector<double> &y,
                                                  const double step) {
    const std::size_t N = y.size();
    const double step2 = step * step;
    std::vector<double> alpha(N - 1), betta(N), gamma(N - 1), delta(N), x_vec(N);
    betta[0] = 1.;
    gamma[0] = 0.;
    delta[0] = 0.;
    for (std::size_t i = 1; i < N - 1; ++i) {
        const double y_der = (y[i + 1] - y[i - 1]) / 2 / step;
        const double y_der2 = (y[i + 1] - 2 * y[i] + y[i - 1]) / step2;

        alpha[i - 1] = 1. - k.df_dy_(t[i], y[i], y_der) * step / 2;
        betta[i] = -2. + k.df_dy(t[i], y[i], y_der) * step2;
        gamma[i] = 1. + k.df_dy_(t[i], y[i], y_der) * step / 2;

        delta[i] = k.rhs(t[i], y[i], y_der, y_der2) * step2;
    }
    betta[N - 1] = 1.;
    alpha[N - 2] = 0.;
    delta[N - 1] = 0.;

    const auto res = method_progonky(alpha, betta, gamma, delta);
    return res;
}

template<typename Callable, typename Callable_start>
[[nodiscard]] std::tuple<std::vector<double>, std::vector<double> > quasi_linearization(
    const Callable &k, const Callable_start &y0, const auto &state_begin,
    const auto &state_end, const std::size_t N, const double eps) {
    const auto [t_begin, y_begin] = state_begin;
    const auto [t_end, y_end] = state_end;

    const double step = (t_end - t_begin) / (N - 1);
    std::vector<double> t(N), y(N);

    double t_i = t_begin;
    for (std::size_t i = 0; i < N; ++i) {
        t[i] = t_i;
        y[i] = y0.value(t_i);
        t_i += step;
    }

    double tol = 0;
    std::vector<double> dy = progonka_for_ql(k, t, y, step);

    for (std::size_t i = 0; i < N; ++i) {
        y[i] = y[i] + dy[i];
        tol += dy[i];
    }
    while (std::sqrt(std::abs(tol)) > eps) {
        tol = 0;
        dy = progonka_for_ql(k, t, y, step);
        for (std::size_t i = 0; i < N; ++i) {
            y[i] = y[i] + dy[i];
            tol += dy[i];
        }
    }

    return {y, t};
}
