#include <fstream>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>

#include "cmath/boundary-value problem/progonka.hpp"

struct Callable {
    [[nodiscard]] double a(const double x) const {
        return x * x - 3;
    }

    [[nodiscard]] double b(const double x) const {
        return (x * x - 3) * std::cos(x);
    }

    [[nodiscard]] double c(const double x) const {
        const double x2 = x * x;
        const double x3 = x2 * x;
        const double x4 = x3 * x;
        const double e_x = std::exp(x);
        const double cos_x = std::cos(x);
        const double sin_x = std::sin(x);

        return 2 - 6 * x + 2 * x3 + (x2 - 3) * e_x * sin_x * (1 + cos_x) + cos_x * (e_x + (x2 - 1) + x4 - 3 * x2);
    }
};

TEST(progonka, progonka) {
    constexpr double pi = std::numbers::pi_v<double>;
    const Callable k{};
    constexpr std::tuple<double, double> state_begin{0., 0.};
    constexpr std::tuple<double, double> state_end{pi, pi * pi};
    constexpr std::size_t N = 1000;
    const auto [y_res, x_res] = progonka(k, state_begin, state_end, N);
    std::ofstream res_data("progonka.csv");
    res_data << "x,y" << std::endl;
    for (std::size_t i = 0; i < N; ++i) {
        res_data << x_res[i] << "," << y_res[i] << std::endl;
    }
}
