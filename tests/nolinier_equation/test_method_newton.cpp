#include <gtest/gtest.h>
#include <cmath>
#include "cmath/nolinier_equation/method_newton.hpp"
#include "primitives/vector_from_vector.h"
#include "primitives/matrix_from_vector.h"

class Equation {
    double (*fun_)(double);

    double (*fun_deriv_)(double);

public:
    Equation(double (*fun)(double), double (*fun_deriv)(double)) : fun_(fun), fun_deriv_(fun_deriv) {
    }

    [[nodiscard]] double fun(const double x) const { return fun_(x); }
    [[nodiscard]] double reverse_fun_deriv(const double x) const { return 1. / fun_deriv_(x); }
    [[nodiscard]] double norm(const double x, const double x_next) const { return std::abs(x - x_next); }
};

template<typename T>
class SystemEquation {
    Vector<T> (*fun_)(const Vector<T> &);

    Matrix<T> J_rev;

public:
    SystemEquation(Vector<T> (*fun)(const Vector<T> &), const Matrix<T> &J_rev) : fun_(fun), J_rev(J_rev) {
    }

    [[nodiscard]] Vector<T> fun(const Vector<T> &x) const { return fun_(x); }
    [[nodiscard]] Matrix<T> reverse_fun_deriv(const Vector<T>& x) const { return J_rev; }
    [[nodiscard]] double norm(const Vector<T>& x, const Vector<T>& x_next) const { return std::abs((x - x_next).norm()); }
};

TEST(method_newton, equation) {
    constexpr double eps = 1e-3;

    const auto fun = [](const double x) {
        return x * std::exp(-x * x) - std::exp(-0.5) / (2 * std::sqrt(2));
    };
    const auto fun_deriv = [](const double x) {
        const double x2 = x * x;
        return std::exp(-x2) * (1 - 2 * x2);
    };

    const Equation eq{fun, fun_deriv};
    constexpr double x1_0 = 0.2;
    const double x1 = method_newton(eq, x1_0, eps);
    ASSERT_NEAR(x1, 0.226, eps);

    constexpr double x2_0 = 1.5;
    const double x2 = method_newton(eq, x2_0, eps);
    ASSERT_NEAR(x2, 1.359, eps);
}

TEST(method_newton, system) {
    constexpr double eps = 1e-6;

    const auto fun = [](const Vector<double>& u) {
        const double x = u(0);
        const double y = u(1);
        return Vector<double>{{(x * x + y * y - 1), y - std::tan(x)}, 2};
    };
    const auto fun_deriv = [](const Vector<double>& u) {
        const double x = u(0);
        const double y = u(1);
        const double cos_x = std::cos(x);
        const double cos_x2 = cos_x * cos_x;

        const double det = 2 *(x + y / cos_x2);
        const Matrix<double> matrix{{1, -2 * y, 1 / cos_x2, 2 * x}, 2, 2};

        return Matrix<double>{matrix * (1 / det)};
    };

    const Vector<double> x1_0{{0.5, 0.5}, 2};
    const SystemEquation<double> eq_system1{fun, fun_deriv(x1_0)};
    const Vector<double> res1 = method_newton(eq_system1, x1_0, eps);
    ASSERT_NEAR(res1(0), 0.649889, eps);
    ASSERT_NEAR(res1(1), 0.760029, eps);

    const Vector<double> x2_0{{-0.5, -0.5}, 2};
    const SystemEquation<double> eq_system{fun, fun_deriv(x2_0)};
    const Vector<double> res2 = method_newton(eq_system, x2_0, eps);
    ASSERT_NEAR(res2(0), -0.649889, eps);
    ASSERT_NEAR(res2(1), -0.760029, eps);
}