#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

#include "cmath/boundary-value problem/quasi_linearization.hpp"

struct Callable_start {
    [[nodiscard]] double value(const double t) const {
        return t * std::log(t);
    }
};

struct Callable {
    [[nodiscard]] double operator()(const double t, const double x0, const double x1) const {
        return std::sqrt(-std::exp(x1) * x0 + std::exp(1.) * x0 * x0 / std::log(t) + 1 / (t * t));
    }

    [[nodiscard]] double df_dy(const double t, const double y, const double y_) const {
        return -(2 * y * std::exp(1.) / std::log(t) - std::exp(y_)) / 2 / this->operator()(t, y, y_);
    }

    [[nodiscard]] double df_dy_(const double t, const double y, const double y_) const {
        return y * std::exp(y_) / 2 / this->operator()(t, y, y_);
    }

    [[nodiscard]] double rhs(const double t, const double y, const double y_, const double y_2) const {
        return this->operator()(t, y, y_) - y_2;
    }
};

TEST(test, test) {
    constexpr double pi = std::numbers::pi_v<double>;
    constexpr Callable_start y0{};
    constexpr Callable k{};
    const double e = std::exp(1.);
    const double e2 = e * e;
    const std::tuple<double, double> state_begin{e, e};
    const std::tuple<double, double> state_end{e2, 2 * e2};
    constexpr std::size_t N = 1e3;
    constexpr double eps = 1e-12;

    const auto [res, t] = quasi_linearization(k, y0, state_begin, state_end, N, eps);
    std::ofstream res_data("quasi_linearization.csv");
    res_data << "t,y" << std::endl;
    for (std::size_t i = 0; i < res.size(); i += 1) {
        res_data << t[i] << "," << res[i] << std::endl;
    }
}
