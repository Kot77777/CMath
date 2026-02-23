#include "Corner.hpp"

std::vector<double> Corner::operator()(const std::vector<double> &u, const Grid &grid) const {
    std::vector<double> res(u.size());
    for (std::size_t i = 1; i < res.size(); ++i) {
        res[i] = u[i] * (1 - grid.CFL()) + grid.CFL() * u[i - 1];
    }
    res[0] = res.back();
    return res;
}
