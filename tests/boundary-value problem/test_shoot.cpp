#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

#include "cmath/boundary-value problem/shoot.hpp"
#include "cmath/spline/spline.hpp"

struct RHS_bound1 {
    using VecType = double;

    struct State {
        VecType x;
        double t;
    };

    [[nodiscard]] double operator()(const double t, const double x0, const double x1) const {
        return -std::sqrt(-std::exp(x1) * x0 + std::exp(1.) * x0 * x0 / std::log(t) + 1 / (t * t));
    }

    [[nodiscard]] double derivative_(const double t, const double x0, const double x1) const {
        return -x0 * std::exp(x1) / 2 / this->operator()(t, x0, x1);
    }

    [[nodiscard]] double derivative(const double t, const double x0, const double x1) const {
        return (2 * x0 * std::exp(1.) / std::log(t) - std::exp(x1)) / 2 / this->operator()(t, x0, x1);
    }
};

TEST(shoot, shoot1) {
    const double e = std::exp(1.);
    const double e2 = e * e;
    const RHS_bound1 k;
    const RHS_bound1::State start{e, e};
    const RHS_bound1::State end{e2 * 2, e2};
    constexpr double step = 1e-3;
    constexpr double alpha0 = 1.;
    const auto res = shoot(k, start, end, alpha0, step, 1e-9);
    std::vector<double> y(res.size());
    std::vector<double> t(res.size());
    for (std::size_t i = 0; i < res.size(); ++i) {
        y[i] = res[i].x(0);
        t[i] = res[i].t;
    }
    const Spline sp{t, y};
    std::ofstream res_data("shoot_T4.csv");
    res_data << "y,x" << std::endl;
    for (const auto &i: res) {
        res_data << sp(i.t) << "," << i.t << std::endl;
    }
    std::cout << "x = 3.0: " << sp(3.) << '\n';
    std::cout << "x = 3.5: " << sp(3.5) << '\n';
    std::cout << "x = 4.0: " << sp(4.0) << '\n';
    std::cout << "x = 4.5: " << sp(4.5) << '\n';
    std::cout << "x = 5.0: " << sp(5.5) << '\n';
}

struct RHS_bound2 {
    using VecType = double;

    struct State {
        VecType x;
        double t;
    };

    [[nodiscard]] double operator()(const double t, const double x0, const double x1) const {
        return -t * std::sqrt(x0);
    }

    [[nodiscard]] double derivative_(const double t, const double x0, const double x1) const {
        return 0.;
    }

    [[nodiscard]] double derivative(const double t, const double x0, const double x1) const {
        return -t / 2 / std::sqrt(x0);
    }
};

TEST(shoot, shoot2) {
    const RHS_bound2 k;
    const RHS_bound2::State start{0., 0.};
    const RHS_bound2::State end{2., 1.};
    constexpr double step = 1e-3;
    constexpr double alpha0 = 1.;
    const auto res = shoot(k, start, end, alpha0, step, 1e-9);
    std::vector<double> y(res.size());
    std::vector<double> t(res.size());
    for (std::size_t i = 0; i < res.size(); ++i) {
        y[i] = res[i].x(0);
        t[i] = res[i].t;
    }
    const Spline sp{t, y};
    std::ofstream res_data("shoot_9_3.csv");
    res_data << "y,x" << std::endl;
    for (const auto &i: res) {
        res_data << sp(i.t) << "," << i.t << std::endl;
    }
}
