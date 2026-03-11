#pragma once
#include "eigen3/Eigen/Eigen"

struct State {
    double rho_{};
    double u_{};
    double e_{};
    double p_{};
    double c_{};
    Eigen::Vector3d w_;

    State() = default;
    State(const double rho, const double u, const double e, const double gamma) : rho_(rho), u_(u), e_(e),
        p_((gamma - 1) * rho * e), c_(gamma * (gamma - 1) * e), w_({rho_, rho_ * u_, rho_ * e_}) {
    }

    State(const Eigen::Vector3d& w, const double gamma) : w_(w) {
        rho_ = w(0);
        u_ = w(1) / rho_;
        e_ = w(2) / rho_;
        p_ = (gamma - 1) * rho_ * e_;
        c_ = std::sqrt(gamma * (gamma - 1) * e_);
    }
};
