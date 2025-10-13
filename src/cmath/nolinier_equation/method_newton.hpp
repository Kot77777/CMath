#pragma once

#include <cmath>
#include <utility>

template<typename T, typename SomeEquation>
[[nodiscard]] T method_newton(const SomeEquation& some, const T& x0, double eps) {
    T x = x0, x_next = x0 - some.reverse_fun_deriv(x0) * some.fun(x0);
    while (some.norm(x, x_next) > eps) {
        x = std::exchange(x_next, x_next - some.reverse_fun_deriv(x_next) * some.fun(x_next));
    }
    return x_next;
}