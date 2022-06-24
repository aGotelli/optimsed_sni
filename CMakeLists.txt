cmake_minimum_required(VERSION 3.16)

project(optimized_sni LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)


include_directories(include)

find_package(Eigen3 REQUIRED)

find_package(benchmark REQUIRED)




add_library(ode45_numerical_integration  include/rod.h
                                         src/rod.cpp    )
target_link_libraries(ode45_numerical_integration Eigen3::Eigen)



add_library(spectral_numerical_integration  include/chebyshev_differentiation.h
                                            include/integrator_base.h
                                            include/quaternion_integrator.h
                                            src/quaternion_integrator.cpp
                                            include/position_integrator.h
                                            src/position_integrator.cpp         )
target_link_libraries(spectral_numerical_integration Eigen3::Eigen)


add_executable(integrate main.cpp)
target_link_libraries(integrate spectral_numerical_integration)

add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark spectral_numerical_integration ode45_numerical_integration benchmark::benchmark)
