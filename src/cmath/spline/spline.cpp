#include "spline.hpp"
#include "solution_SLAE/method_progonky.h"

Spline::Spline(const std::vector<double> &x, const std::vector<double> &f) : x_(x), f_(f), N(x.size() - 1) {
    b_.resize(N);
    c_.resize(N);
    d_.resize(N);
    find_polynomials();
}

void Spline::find_polynomials() {
    const auto vec_transform = [](const auto &rule, const std::size_t N) {
        std::vector<double> res(N);
        for (std::size_t i = 0; i < N; ++i) {
            res[i] = rule(i);
        }
        return res;
    };

    const std::vector<double> h_ = vec_transform([this](const std::size_t i) {
                                                            return x_[i + 1] - x_[i];}, N);
    const std::vector<double> diag(N - 1, 2);

    const std::vector<double> diag_up = vec_transform([&h_](const std::size_t i) {
                                                             return h_[i + 1] / (h_[i + 1] + h_[i]);}, N - 2);

    const std::vector<double> diag_down = vec_transform([&h_](const std::size_t i) {
                                                             return h_[i + 1] / (h_[i + 2] + h_[i + 1]);}, N - 2);

    const std::vector<double> U_1 = vec_transform([&h_, this](const std::size_t i) {
                                                             return (f_[i + 1] - f_[i]) / h_[i];}, N);

    const std::vector<double> U_2 = vec_transform([&U_1, &h_](const std::size_t i) {
                                                             return 6 * (U_1[i + 1] - U_1[i]) / (h_[i + 1] + h_[i]);}, N - 1);

    const std::vector<double> res = method_progonky(diag_down, diag, diag_up, U_2);
    std::copy(res.begin(), res.end(), c_.begin());

    b_[0] = c_[0] * h_[0] / 3 + U_1[0];
    d_[0] = c_[0] / h_[0];
    for (std::size_t i = 1; i <= N - 1; ++i) {
        b_[i] = (c_[i] + c_[i - 1] / 2) * h_[i] / 3 + U_1[i];
        d_[i] = (c_[i] - c_[i - 1]) / h_[i];
    }
}

double Spline::operator()(const double x) const {
    const auto find_i = [&]() {
        auto i = static_cast<std::size_t>(x == x_[0]);
        while (x > x_[i]) { i++;}
        return i;
    };

    const auto i = x > x_.back() ? N : x < x_[0] ? 1 : find_i();

    const double dif = x - x_[i];
    const double dif2 = dif * dif;
    const double dif3 = dif * dif2;

    return f_[i] + b_[i - 1] * dif + c_[i - 1] / 2 * dif2 + d_[i - 1] / 6 * dif3;
}
