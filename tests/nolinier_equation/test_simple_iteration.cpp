#include <gtest/gtest.h>
#include <cmath>

#include "cmath/nolinier_equation/simple_iteration.hpp"

TEST(simple_iteration, test) {
    constexpr  double eps = 1e-3;

    const auto fun1 = [](const double x) {
        return std::exp(x * x - 0.5) / (2 * std::sqrt(2));
    };
    constexpr  double x1_0 = 0.25;
    const double x1 = simple_iteration(fun1, x1_0, eps);
    ASSERT_NEAR(x1, 0.226, eps);

    const auto fun2 = [](const double x) {
        return std::sqrt(std::log(2 * x * std::sqrt(2 * std::exp(1))));
    };
    constexpr  double x2_0 = 1.5;
    const double x2 = simple_iteration(fun2, x2_0, eps);
    ASSERT_NEAR(x2, 1.359, eps);

}