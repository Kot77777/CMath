#include <fstream>
#include <iomanip>
#include <gtest/gtest.h>
#include <cmath/diff_equation/Dormand_Prince5.hpp>
#include "Eigen/Eigen"
#include <stdfloat>

#include "cmath/diff_equation/RK4.hpp"

class RHS {
    long double mu_;
    long double eta_;

public:
    RHS(const double mu) : mu_(mu), eta_(1 - mu) {
    }

    using Type = long double;
    using VecType = Eigen::Vector<Type, 4>;

    struct State {
        VecType x;
        Type t;
    };

    [[nodiscard]] VecType operator()(const Type t, const VecType &x) const {
        const Type x0 = x(0), x1 = x(1), x2 = x(2), x3 = x(3);
        const Type a = (x0 + mu_) * (x0 + mu_) + x2 * x2;
        const Type b = (x0 - eta_) * (x0 - eta_) + x2 * x2;
        const Type A = std::sqrt(a * a * a);
        const Type B = std::sqrt(b * b * b);

        const Type dx0_dt = x1;
        const Type dx1_dt = x0 + 2 * x3 - eta_ * (x0 + mu_) / A - mu_ * (x0 - eta_) / B;
        const Type dx2_dt = x3;
        const Type dx3_dt = x2 - 2 * x1 - eta_ * (x2 / A) - mu_ * (x2 / B);

        return {dx0_dt, dx1_dt, dx2_dt, dx3_dt};
    }
};

TEST(DP5, DP5) {
    const RHS rhs{0.012277471L};
    const RHS::State state0{{0.994, 0., 0., -2.00158510637908252240537862224}, 0.};
    const double T = 17.0652165601579625588917206249;
    const double endTime = 5 * T, step = 1e-3;
    constexpr double eps = 1e-14;

    const auto res = DP5(rhs, state0, endTime, step, eps);
    std::ofstream res_data("DP5.csv");
    res_data << std::setprecision(30);
    res_data << "x,u,y,nu,t" << std::endl;
    for (const auto &[x, t]: res) {
        res_data << x(0) << "," << x(1) << "," << x(2) << "," << x(3) << "," << t << std::endl;
    }
}
