#include "Grid.hpp"

Grid::Grid(const double L, const double T, const double h, const double CFL) : L_(L), T_(T), h_(h), CFL_(CFL),
                                                                               tau_(CFL * h) {
}

std::size_t Grid::MX() const {
    return static_cast<std::size_t>(L_ / h_);
}

double Grid::value_X(const std::size_t i) const {
    return i * h_;
}

double Grid::T() const { return T_; }

double Grid::tau() const { return tau_; }

double Grid::CFL() const { return CFL_; }
