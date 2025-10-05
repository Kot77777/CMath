#pragma once
#include <vector>

class Spline {
    std::vector<double> x_;
    std::vector<double> f_;
    std::size_t N;

    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;

    void find_polynomials();

public:
    Spline(const std::vector<double> &x, const std::vector<double> &f);

    [[nodiscard]] double operator()(double x) const;
};
