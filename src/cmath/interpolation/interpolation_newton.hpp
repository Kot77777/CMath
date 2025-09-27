#pragma once
#include <vector>

class Interpolation_Newton {
    std::vector<double> x_;
    std::vector<double> f_;
    std::vector<double> res_up;
    std::vector<double> res_down;
    std::size_t len_arr;

    void update(const std::vector<double>& f);

public:
    Interpolation_Newton(const std::vector<double>& x, const std::vector<double>& f);

    double operator()(double x) const;
};