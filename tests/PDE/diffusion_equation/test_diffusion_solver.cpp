#include <cmath>
#include <format>
#include <fstream>
#include <gtest/gtest.h>

#include "cmath/PDE/diffusion_equation/diffusion_solver.hpp"
#include "cmath/PDE/diffusion_equation/Grid.hpp"
#include "cmath/PDE/diffusion_equation/templates/Corner.hpp"
#include "cmath/PDE/diffusion_equation/templates/Lax-Wendroff.hpp"

TEST(diffusion_solver, corner_template) {
    const double pi = std::numbers::pi_v<double>;
    const double L = 20.;
    const double T = 50.;
    const double h = 0.5;
    const double CFL = 0.6;

    const Grid grid{L, T, h, CFL};
    constexpr Corner corner_template;
    const auto u0 = [&](const double x) {
        return std::sin(4 * pi * x / L);
    };

    const auto res = diffusion_solver(grid, corner_template, u0);

    std::ofstream data(
        "/home/kostya/Repositories/python_lesson/PDE/diffusion_solver/data/diffusion_solver_Corner_temp_cfl_0.6.csv"
    );

    for (const auto & re : res) {
        for (std::size_t i = 0; i < grid.MX(); ++i) {
            data << re[i] << ",";
        }
        data << re.back() << std::endl;
    }
}

TEST(diffusion_solver, Lax_Wendroff_template) {
    const double pi = std::numbers::pi_v<double>;
    const double L = 20.;
    const double T = 50.;
    const double h = 0.5;
    const double CFL = 1.01;

    const Grid grid{L, T, h, CFL};
    constexpr Lax_Wendroff LW_template;
    const auto u0 = [&](const double x) {
        return std::sin(4 * pi * x / L);
    };

    const auto res = diffusion_solver(grid, LW_template, u0);

    std::ofstream data(
        "/home/kostya/Repositories/python_lesson/PDE/diffusion_solver/data/diffusion_solver_LW_temp_cfl_1.01.csv"
    );

    for (const auto & re : res) {
        for (std::size_t i = 0; i < grid.MX(); ++i) {
            data << re[i] << ",";
        }
        data << re.back() << std::endl;
    }
}
