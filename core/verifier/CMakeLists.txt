cmake_minimum_required(VERSION 3.10)
project(CertifiedSearchVerifier CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMP for parallel batch verification
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: version=${OpenMP_CXX_VERSION}, flags=${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP not found; batch verifiers will compile without parallelism")
endif()

# Core library
add_library(verifier_core src/lj_potential.cpp)
target_include_directories(verifier_core PUBLIC src)

# Interval arithmetic spike
find_library(FLINT_LIB NAMES flint PATHS /usr/lib /usr/lib/x86_64-linux-gnu)
find_library(MPFI_LIB NAMES mpfi PATHS /usr/lib /usr/lib/x86_64-linux-gnu)
find_library(MPFR_LIB NAMES mpfr PATHS /usr/lib /usr/lib/x86_64-linux-gnu)
find_library(GMP_LIB NAMES gmp PATHS /usr/lib /usr/lib/x86_64-linux-gnu)

if(FLINT_LIB AND MPFI_LIB AND MPFR_LIB AND GMP_LIB)
    add_executable(spike_mpfi spike/interval_mpfi.cpp)
    target_link_libraries(spike_mpfi ${MPFI_LIB} ${MPFR_LIB} ${GMP_LIB})
    
    add_executable(spike_arb spike/interval_arb.cpp)
    target_link_libraries(spike_arb ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
    
    add_executable(krawczyk_lj13 src/krawczyk_lj13.cpp)
    target_link_libraries(krawczyk_lj13 verifier_core ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
    
    add_executable(krawczyk_lj38 src/krawczyk_lj38.cpp)
    target_link_libraries(krawczyk_lj38 verifier_core ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
    
    # Batch verifier (original, unoptimized baseline)
    add_executable(krawczyk_lj38_batch src/krawczyk_lj38_batch.cpp)
    target_link_libraries(krawczyk_lj38_batch verifier_core ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
    if(OpenMP_CXX_FOUND)
        target_link_libraries(krawczyk_lj38_batch OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(krawczyk_lj38_batch stdc++fs)
    
    # Optimized batch verifier (thread pooling, async I/O, capped threads)
    add_executable(krawczyk_lj38_batch_opt src/krawczyk_lj38_batch_opt.cpp)
    target_link_libraries(krawczyk_lj38_batch_opt verifier_core ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
    if(OpenMP_CXX_FOUND)
        target_link_libraries(krawczyk_lj38_batch_opt OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(krawczyk_lj38_batch_opt stdc++fs)
    
    message(STATUS "Interval libraries found: FLINT=${FLINT_LIB}, MPFI=${MPFI_LIB}")
else()
    message(WARNING "Interval libraries not found, skipping spike targets")
endif()

# Testing
enable_testing()
add_executable(test_lj tests/test_lj_potential.cpp)
target_link_libraries(test_lj verifier_core)
add_test(NAME LJMathTest COMMAND test_lj)