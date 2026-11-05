cmake_minimum_required(VERSION 3.10)
project(CertifiedSearchVerifier CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Core library
add_library(verifier_core src/lj_potential.cpp)
target_include_directories(verifier_core PUBLIC src)

# Testing
enable_testing()
add_executable(test_lj tests/test_lj_potential.cpp)
target_link_libraries(test_lj verifier_core)
add_test(NAME LJMathTest COMMAND test_lj)