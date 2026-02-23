#pragma once
#include <vector>

#include "cmath/PDE/diffusion_equation/Grid.hpp"

struct Corner {
    [[nodiscard]] std::vector<double> operator()(const std::vector<double> &u, const Grid& grid) const;
};
