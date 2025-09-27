file(GLOB_RECURSE SOURCES cmath/*.cpp cmath/*.hpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_23)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})