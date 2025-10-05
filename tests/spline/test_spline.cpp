#include <gtest/gtest.h>
#include <fstream>
#include "cmath/spline/spline.hpp"

TEST(spline, test1) {
    const auto check = [](const double res, const auto exp) {
        ASSERT_NEAR(res, exp, 1e-15);
    }; {
        const std::vector<double> x{0, 1, 2};
        const std::vector<double> f{0, 1, 2};
        const Spline sp{x, f};
        check(sp(0), 0);
        check(sp(1), 1);
        check(sp(2), 2);
        check(sp(0.5), 0.5);
        check(sp(1.5), 1.5);
        check(sp(3), 3);
    } {
        const std::vector<double> x{1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
        const std::vector<double> f{
            92228496, 106021537, 123202624, 132164569, 151325798,
            179323175, 203211926, 226545805, 248709873, 281421906
        };
        const Spline sp{x, f};
        std::cout << "2010 год:" << static_cast<int>(sp(2010));
        std::ofstream res("data_spline.csv");
        res << "year" << "," << "number" << "\n";
        for (double i = 1910; i <= 2010; i += 1) {
            res << i << "," << sp(i) << "\n";
        }
    }
}
