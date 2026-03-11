#pragma once
#include "State.hpp"

struct matrix {
    Eigen::Matrix3d omega_;
    Eigen::Matrix3d omega_inv;
    Eigen::Matrix3d lambda;
    Eigen::Matrix3d lambda_abs;
    Eigen::Matrix3d A_;

    matrix(const State& s, double gamma);
};