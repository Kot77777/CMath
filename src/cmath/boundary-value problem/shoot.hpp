#pragma once
#include <vector>

#include "Eigen/Eigen"
#include "cmath/diff_equation/RK4.hpp"

template<typename Callable>
class RHS {
    const Callable &f;

public:
    RHS(const Callable &f_) : f(f_) {
    }

    using VecType = Eigen::Vector4d;

    struct State {
        VecType x;
        double t;
    };

    [[nodiscard]] VecType operator()(const double t, const VecType &x) const {
        const double x0 = x(0), x1 = x(1), x2 = x(2), x3 = x(3);
        const double dx0_dt = x1;
        const double dx1_dt = -f(t, x0, x1);
        const double dx2_dt = x3;
        const double dx3_dt = -f.derivative_(t, x0, x1) * x3 - f.derivative(t, x0, x1) * x2;

        return VecType{dx0_dt, dx1_dt, dx2_dt, dx3_dt};
    }
};

template<typename Callable>
[[nodiscard]] std::vector<typename RHS<Callable>::State> shoot(const Callable &f,
                                                               const typename Callable::State &startVec,
                                                               const typename Callable::State &endVec,
                                                               const double alpha0,
                                                               const double step, const double eps) {
    const RHS rhs{f};
    const typename RHS<Callable>::State state0{{startVec.x, alpha0, 0., 1.}, startVec.t};
    std::vector<typename RHS<Callable>::State> res = RK4(rhs, state0, endVec.t, step);
    double r = res.back().x(0) - endVec.x;
    double alpha = alpha0;
    while (r > eps) {
        alpha = alpha - r / res.back().x(2);
        const typename RHS<Callable>::State state{{startVec.x, alpha, 0., 1.}, startVec.t};
        res = RK4(rhs, state, endVec.t, step);
        r = res.back().x(0) - endVec.x;
    }
    return res;
}
