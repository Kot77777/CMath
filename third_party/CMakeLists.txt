include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(
        System_of_equations
        GIT_REPOSITORY https://github.com/Kot77777/System_of_equations.git
        GIT_TAG main
)
FetchContent_MakeAvailable(System_of_equations)