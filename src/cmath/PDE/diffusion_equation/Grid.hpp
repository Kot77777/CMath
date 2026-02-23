#pragma once
#include <cstddef>

class Grid {
    double L_;
    double T_;
    double h_;
    double CFL_;
    double tau_;

public:
    Grid(double L, double T, double h, double CFL);

    [[nodiscard]] std::size_t MX() const;
    [[nodiscard]] double value_X(std::size_t i) const;
    [[nodiscard]] double T() const;
    [[nodiscard]] double tau() const;
    [[nodiscard]] double CFL() const;
};
