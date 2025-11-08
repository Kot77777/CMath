#pragma once
#include <vector>

template<typename Callable>
typename Callable::VecType RK4_step(const Callable &f, const typename Callable::State &state, const double step) {
    const auto k1 = f(state.t, state.x);
    const auto k2 = f(state.t + 0.5 * step, state.x + 0.5 * step * k1);
    const auto k3 = f(state.t + 0.5 * step, state.x + 0.5 * step * k2);
    const auto k4 = f(state.t + step, state.x + step * k3);

    const typename Callable::VecType x_next = state.x + step / 6 * (k1 + k4 + 2 * (k2 + k3));
    return x_next;
}

template<typename Callable>
std::vector<typename Callable::State> RK4(const Callable &f, const typename Callable::State &state0,
                                          const double endTime,
                                          const double step) {
    std::vector<typename Callable::State> res{state0};
    while (res.back().t < endTime) {
        res.push_back({RK4_step<Callable>(f,res.back(), step), res.back().t + step});
    }
    return res;
}
