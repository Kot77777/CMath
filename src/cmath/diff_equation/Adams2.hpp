#pragma once
#include <vector>
#include "cmath/diff_equation/RK4.hpp"

template<typename Callable>
[[nodiscard]] std::vector<typename Callable::State> Adams2(const Callable &f, const typename Callable::State &state0,
                                                           const double endTime,
                                                           const double step) {
    std::vector<typename Callable::State> res{state0};
    res.push_back({RK4_step<Callable>(f, res.back(), step), res.back().t + step});
    std::size_t i = 1;
    while (res.back().t < endTime) {
        const double b1 = f(res[i].t, res[i].x);
        const double b2 = f(res[i-1].t, res[i-1].x);
        res.push_back({res[i].x + step * (1.5 * b1 - 0.5 * b2), res.back().t + step});
        i += 1;
    }
    return res;
}
