#include "interpolation_newton.hpp"

Interpolation_Newton::Interpolation_Newton(const std::vector<double> &x, const std::vector<double> &f) : x_(x), f_(f) {
    len_arr = x.size();
    res_down.resize(len_arr);
    res_up.resize(len_arr);
    update(f_);
}

void Interpolation_Newton::update(const std::vector<double>& f) {
    if (const std::size_t N = f.size(); N != 0) {
        const std::size_t p = len_arr - N;
        res_up[p] = f[0];
        res_down[p] = f[N - 1];
        std::vector<double> cache(N - 1);
        for (std::size_t i = 0; i < N - 1; ++i) {
            cache[i] = (f[i+1] - f[i]) / (x_[i + p + 1] - x_[i]);
        }
        update(cache);
    }
}

[[nodiscard]] double Interpolation_Newton::operator()(const double x) const {
    double res = res_up[0];
    for (std::size_t i = 1; i < len_arr; ++i) {
        double p = res_up[i];
        for (std::size_t j = 0; j < i; ++j) {
            p *= x - x_[j];
        }
        res += p;
    }
    return res;
}


