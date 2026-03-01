#include "Lax-Wendroff.hpp"

std::vector<double> Lax_Wendroff::operator()(const std::vector<double> &u, const Grid &grid) const {
    std::vector<double> res(u.size());
    const auto LW = [&grid](const double u_next, const double u_curr, const double u_prev) {
        return u_curr - grid.CFL() / 2 * (u_next - u_prev) + grid.CFL() * grid.CFL() / 2 * (
                   u_next - 2 * u_curr + u_prev);
    };
    const std::size_t N = u.size();
    for (std::size_t i = 1; i < N - 1; ++i) {
        res[i] = LW(u[i + 1], u[i], u[i - 1]);
    }
    res[N - 1] = LW(u[0], u[N - 1], u[N - 2]);
    res[0] = LW(u[1], u[0], u[N - 1]);
    // res[0] = res[N - 1];
    return res;
}
