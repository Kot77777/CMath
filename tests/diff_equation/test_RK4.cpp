#include <fstream>
#include "Eigen/Eigen"
#include "cmath/diff_equation/RK4.hpp"

class RK4_RHS {
    double mu_;

public:
    RK4_RHS(const double mu) : mu_(mu) {
    }

    using VecType = Eigen::Vector2d;

    struct State {
        VecType x;
        double t;
    };

    [[nodiscard]] VecType operator()(const double t, const VecType &x) const {
        const double x0 = x(0), x1 = x(1);
        const double dx0_dt = x1;
        const double dx1_dt = mu_ * (1 - x1 * x1) * x1 - x0;

        return VecType{dx0_dt, dx1_dt};
    }
};

int main() {
    constexpr double mu = 1e3;
    const RK4_RHS rhs{mu};
    const RK4_RHS::State state0{{0, 0.001}, 0};
    constexpr double step = 0.001;
    constexpr double endTime = 1000.;

    const auto res = RK4(rhs, state0, endTime, step);
    std::ofstream res_data("RK4.csv");
    res_data << "x,t" << std::endl;
    for (const auto&[x, t]: res) {
        res_data << x(0) << "," << t << std::endl;
    }
}