file(GLOB_RECURSE SOURCES cmath/*.cpp cmath/*.hpp)

add_library(CMath STATIC ${SOURCES})
target_link_libraries(CMath PUBLIC SLAE_project)
target_compile_features(CMath PUBLIC cxx_std_23)

target_include_directories(CMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})