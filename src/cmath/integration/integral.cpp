#include "cmath/integration/integral.hpp"

double Integral::rectangle(const double h) const {
    double res = 0;
    for (double x = x_start; x < x_end; x += h) {
        res += fun(x);
    }
    return res * h;
}

double Integral::mean_dot(const double h) const {
    double res = 0;
    for (double x = x_start; x < x_end; x += h) {
        res += fun(x + h / 2);
    }
    return res * h;
}

double Integral::trapezoid(const double h) const {
    double res = 0;
    for (double x = x_start; x < x_end; x += h) {
        res += (fun(x) + fun(x + h)) / 2;
    }
    return res * h;
}

double Integral::simpson(const double h) const {
    double res = 0;
    for (double x = x_start; x < x_end; x += h) {
        res += (fun(x) + 4 * fun(x + h / 2) + fun(x + h)) / 6;
    }
    return res * h;
}

double Integral::threeeights(const double h) const {
    double res = 0;
    for (double x = x_start; x < x_end; x += h) {
        const double h_3 = h / 3;
        res += (fun(x) + 3 * fun(x + h_3) + 3 * fun(x + 2 * h_3) + fun(x + h)) / 8;
    }
    return res * h;
}
