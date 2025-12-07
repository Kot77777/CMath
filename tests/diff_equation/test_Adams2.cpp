#include <fstream>
#include <gtest/gtest.h>
#include "cmath/diff_equation/Adams2.hpp"
#include <Eigen/Eigen>

struct Adams2_RHS {
    using VecType = double;

    struct State {
        VecType x;
        double t;
    };

    [[nodiscard]] VecType operator()(const double t, const VecType &x) const {
        return -2 * x;
    }
};

TEST(test, test) {
    constexpr Adams2_RHS rhs;
    constexpr Adams2_RHS::State state0{-2., 0.};
    constexpr double step = 1e-3;
    const auto res = Adams2(rhs, state0, 100., step);
    std::ofstream res_data("Adams2.csv");
    res_data << "x_th,x,t" << std::endl;
    for (const auto&[x, t]: res) {
        res_data << -2 * std::exp(-2 * t) << "," << x << "," << t << std::endl;
    }
}
