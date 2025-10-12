#pragma once

class Integral {
    double (*fun)(double);
    double x_start;
    double x_end;

public:
    Integral(const auto& fun, const double x_start, const double x_end) : fun(fun), x_start(x_start), x_end(x_end){}

    [[nodiscard]] double rectangle(double h) const;
    [[nodiscard]] double mean_dot(double h) const;
    [[nodiscard]] double trapezoid(double h) const;
    [[nodiscard]] double simpson(double h) const;
    [[nodiscard]] double threeeights(double h) const;
};
