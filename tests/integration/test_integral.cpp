#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

#include "cmath/integration/integral.hpp"

TEST(integral, test1) { {
        const auto fun = [](const double x) {
            return 1.;
        };
        const Integral J(fun, 0, 2);
        ASSERT_EQ(J.rectangle(2.), 2);
    } {
        const auto fun = [](const double x) {
            return x;
        };
        const Integral J(fun, 0, 2);
        ASSERT_EQ(J.mean_dot(2.), 2);
    } {
        const auto fun = [](const double x) {
            return x;
        };
        const Integral J(fun, 0, 2);
        ASSERT_EQ(J.trapezoid(2.), 2);
    } {
        const auto fun = [](const double x) {
            return x * x * x;
        };
        const Integral J(fun, 0, 2);
        ASSERT_EQ(J.simpson(2.), 4);
    } {
        const auto fun = [](const double x) {
            return x * x * x;
        };
        const Integral J(fun, 0, 2);
        ASSERT_EQ(J.threeeights(2.), 4);
    }
}

TEST(integral, test2) {
    const auto fun = [](const double x) {
        return std::sin(100 * x) * std::exp(-x * x) * std::cos(2 * x);
    };

    const Integral J{fun, 0, 3};
    std::ofstream res("integral.txt");
    res.precision(10);

    const double rectangle = J.rectangle(1e-8);
    res << "Метод левых прямоугольников: " << rectangle << "\n";

    const double mean_dot = J.mean_dot(1e-7);
    res << "Метод прямоугольников с центральной точкой: " << mean_dot << "\n";

    const double trapezoid = J.trapezoid(1e-7);
    res << "Метод трапеции: " << trapezoid << "\n";

    const double simpson = J.simpson(1e-4);
    res << "Метод Симпсона: " << simpson << "\n";

    const double threeeights = J.threeeights(1e-4);
    res << "Правило '3/8': " << threeeights << "\n";
}
