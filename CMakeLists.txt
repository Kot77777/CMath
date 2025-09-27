CMAKE_MINIMUM_REQUIRED(VERSION 3.20)
project(CMath)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CMATH_WITH_TESTS "Build tests" OFF)

add_subdirectory(src)
add_subdirectory(third_party)

if (CMATH_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
