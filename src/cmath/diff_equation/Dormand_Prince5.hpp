#pragma once
#include <cmath>
#include <vector>

template<typename Callable>
[[nodiscard]] std::tuple<typename Callable::VecType, typename Callable::Type> DP5_step(
    const Callable &f, const typename Callable::State &state,
    const typename Callable::Type step) {
    const auto k1 = f(state.t, state.x);
    const auto k2 = f(state.t + 0.2 * step, state.x + 0.2 * step * k1);
    const auto k3 = f(state.t + 0.3 * step, state.x + step * (3. / 40 * k1 + 9. / 40 * k2));
    const auto k4 = f(state.t + 0.8 * step, state.x + step * (44. / 45 * k1 - 56. / 15 * k2 + 32. / 9 * k3));
    const auto k5 = f(state.t + 8. / 9 * step,
                      state.x + step * (19372. / 6561 * k1 - 25360. / 2187 * k2 + 64448. / 6561 * k3 - 212. / 729 *
                                        k4));
    const auto k6 = f(state.t + step,
                      state.x + step * (9017. / 3168 * k1 - 355. / 33 * k2 + 46732. / 5247 * k3 + 49. / 176 * k4 - 5103.
                                        / 18656 * k5));
    const auto k7 = f(state.t + step,
                      state.x + step * (35. / 384 * k1 + 500. / 1113 * k3 + 125. / 192 * k4 - 2187. / 6784 * k5 + 11. /
                                        84 * k6));

    const typename Callable::VecType x_next =
            state.x + step * (35. / 384 * k1 + 500. / 1113 * k3 + 125. / 192 * k4 - 2187. / 6784 * k5 + 11. / 84 * k6);

    const typename Callable::VecType x_next_ =
            state.x + step * (5179. / 57600 * k1 + 7571. / 16695 * k3 + 393. / 640 * k4 - 92097. / 339200 * k5 + 187. /
                              2100 * k6 + 1. / 40 * k7);

    const typename Callable::Type eps = (x_next_ - x_next).norm();

    return {x_next, eps};
}

template<typename Callable>
[[nodiscard]] std::vector<typename Callable::State> DP5(const Callable &f, const typename Callable::State &state0,
                                                        const typename Callable::Type endTime,
                                                        const typename Callable::Type step, const double tol) {
    std::vector<typename Callable::State> res{state0};
    const typename Callable::Type step_min = 1e-6, step_max = 1.;
    typename Callable::Type step_next = step;
    while (res.back().t < endTime) {
        const auto [res_i, eps] = DP5_step(f, res.back(), step_next);
        if (eps <= tol) {
            res.push_back({res_i, res.back().t + step_next});
            step_next = std::clamp(0.7 * step_next * std::pow(tol / eps, 0.2), step_min, step_max);
        }
        else {
            step_next = std::clamp(0.7 * step_next * std::pow(tol / eps, 0.2), step_min, step_max);
        }

    }
    return res;
}
