#include <gtest/gtest.h>
#include <fstream>
#include "cmath/interpolation/interpolation_newton.hpp"

TEST(interpolation_newton, create) {
    const auto check = [](const double res, const double exp) {
        ASSERT_EQ(res, exp);
    }; {
        const Interpolation_Newton pol{{0, 1, 2, 3}, {0, 1, 4, 9}};
        check(pol(1.), 1.);
        check(pol(2.), 4.);
        check(pol(4.), 16.);
    } {
        const Interpolation_Newton pol{{0, 1, 2, 3}, {0, 1, 8, 27}};
        check(pol(1.), 1.);
        check(pol(2.), 8.);
        check(pol(4.), 64.);
    } {
        const Interpolation_Newton pol{{0, 1, 2, 3}, {0, 1, 16, 81}};
        check(pol(1.), 1.);
        check(pol(2.), 16.);
        check(pol(4.), 232.);
    } {
        const Interpolation_Newton pol{
            {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000},
            {
                92228496, 106021537, 123202624, 132164569, 151325798,
                179323175, 203211926, 226545805, 248709873, 281421906
            }
        };
        std::cout << "2010 год:" << static_cast<int>(pol(2010));
        std::ofstream res("data_interpolation.csv");
        res << "year" << "," << "number" << "\n";
        for (double i = 1910; i <= 2010; i += 1) {
            res << i << "," << pol(i) << "\n";
        }
    }
}
