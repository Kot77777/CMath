#pragma once
#include <vector>

#include "Grid.hpp"

template<typename Template, typename U0>
std::vector<std::vector<double> > diffusion_solver(const Grid &grid, const Template &temp, const U0 &u0) {
    std::vector<std::vector<double> > res;
    std::vector<double> u0_(grid.MX() + 1);
    for (std::size_t i = 0; i <= grid.MX(); ++i) {
        u0_[i] = u0(grid.value_X(i));
    }
    res.push_back(u0_);
    for (double t = grid.tau(); t <= grid.T(); t += grid.tau()) {
        res.push_back(temp(res.back(), grid));
    }
    return res;
}
