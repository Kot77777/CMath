#include "matrix.hpp"

matrix::matrix(const State &s, const double gamma) {
    omega_ = Eigen::Matrix3d{
        {-s.u_ * s.c_, s.c_, gamma - 1},
        {-s.c_ * s.c_, 0., gamma - 1},
        {s.u_ * s.c_, -s.c_, gamma - 1}
    };

    const double c2 = 2 * s.c_ * s.c_;
    omega_inv = Eigen::Matrix3d{
        {1. / c2, -2 / c2, 1 / c2},
        {(s.u_ + s.c_) / c2, -2 * s.u_ / c2, (s.u_ - s.c_) / c2},
        {1. / 2 / (gamma - 1), 0., 1. / 2 / (gamma - 1)}
    };

    lambda = Eigen::Matrix3d{
        {s.u_ + s.c_, 0., 0.},
        {0., s.u_, 0.},
        {0., 0., s.u_ - s.c_}
    };

    lambda_abs = Eigen::Matrix3d{
        {std::abs(s.u_ + s.c_), 0., 0.},
        {0., std::abs(s.u_), 0.},
        {0., 0., std::abs(s.u_ - s.c_)}
    };

    A_ = Eigen::Matrix3d{
        {0., 1., 0.},
        {-s.u_ * s.u_, 2 * s.u_, gamma - 1},
        {-gamma * s.u_ * s.e_, gamma * s.e_, s.u_}
    };
    // A_ = omega_inv * lambda * omega_;
}
