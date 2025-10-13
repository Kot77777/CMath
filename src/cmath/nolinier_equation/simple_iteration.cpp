#include "simple_iteration.hpp"
#include <cmath>
#include <utility>

double simple_iteration(double (*fun)(double), const double x0, const double eps) {
    double x = x0, x_next = fun(x0);
    while (std::abs(x - x_next) > eps) {
        x = std::exchange(x_next, fun(x_next));
    }
    return x_next;
}
