cmake_minimum_required(VERSION 3.18)

# Configure in a separate build directory for each compiler. CMake cannot use
# icpc and g++ in the same build tree:
#   cmake -S REOS_UnitTests -B build/unit-gcc -DCMAKE_CXX_COMPILER=g++
#   cmake -S REOS_UnitTests -B build/unit-icpc -DCMAKE_CXX_COMPILER=icpc
project(RF_EMT_UnitTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # Tests use compiler extensions and x86 intrinsics.

option(RF_EMT_REGISTER_TESTS "Register executables with CTest (some need specific CPU features)" OFF)

get_filename_component(RF_EMT_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(RF_EMT_SRC_ROOT "${RF_EMT_ROOT}/src")

# Unit tests include headers directly from many module directories.
file(GLOB_RECURSE _rf_emt_headers CONFIGURE_DEPENDS
    "${RF_EMT_SRC_ROOT}/*.h"
    "${RF_EMT_SRC_ROOT}/*.hpp"
    "${RF_EMT_SRC_ROOT}/*.inl"
)
set(_rf_emt_include_dirs "${RF_EMT_SRC_ROOT}")
foreach(_header IN LISTS _rf_emt_headers)
    get_filename_component(_include_dir "${_header}" DIRECTORY)
    list(APPEND _rf_emt_include_dirs "${_include_dir}")
endforeach()
list(REMOVE_DUPLICATES _rf_emt_include_dirs)

# Look up implementation files that match directly included project headers,
# e.g. GMS_trapezoid_waveform.h -> src/modulation/GMS_trapezoid_waveform.cpp.
file(GLOB_RECURSE _rf_emt_cpp_sources CONFIGURE_DEPENDS "${RF_EMT_SRC_ROOT}/*.cpp")
set(_rf_emt_cpp_source_names)
foreach(_source IN LISTS _rf_emt_cpp_sources)
    get_filename_component(_source_name "${_source}" NAME)
    list(APPEND _rf_emt_cpp_source_names "${_source_name}")
endforeach()

file(GLOB _rf_emt_test_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/unit_test_*.cpp")
if(NOT _rf_emt_test_sources)
    message(FATAL_ERROR "No unit_test_*.cpp sources found in ${CMAKE_CURRENT_LIST_DIR}")
endif()

find_package(Threads REQUIRED)
find_package(OpenMP QUIET COMPONENTS CXX)

foreach(_test_source IN LISTS _rf_emt_test_sources)
    get_filename_component(_test_name "${_test_source}" NAME_WE)
    add_executable("${_test_name}" "${_test_source}")
    target_include_directories("${_test_name}" PRIVATE ${_rf_emt_include_dirs})
    target_link_libraries("${_test_name}" PRIVATE Threads::Threads)

    # Add matching implementation units for project headers included directly
    # by this test. Any additional/transitive sources can be supplied with
    # RF_EMT_TEST_SUPPORT_SOURCES below.
    file(STRINGS "${_test_source}" _include_lines REGEX "^[ \t]*#[ \t]*include[ \t]*[<\"]")
    set(_test_support_sources)
    foreach(_include_line IN LISTS _include_lines)
        if(_include_line MATCHES "#[ \t]*include[ \t]*[<\"]([^>\"]+)[>\"]")
            get_filename_component(_include_name "${CMAKE_MATCH_1}" NAME_WE)
            set(_matching_source_name "${_include_name}.cpp")
            if(_matching_source_name IN_LIST _rf_emt_cpp_source_names)
                foreach(_candidate IN LISTS _rf_emt_cpp_sources)
                    get_filename_component(_candidate_name "${_candidate}" NAME)
                    if(_candidate_name STREQUAL _matching_source_name)
                        list(APPEND _test_support_sources "${_candidate}")
                    endif()
                endforeach()
            endif()
        endif()
    endforeach()
    if(RF_EMT_TEST_SUPPORT_SOURCES)
        list(APPEND _test_support_sources ${RF_EMT_TEST_SUPPORT_SOURCES})
    endif()
    if(_test_support_sources)
        list(REMOVE_DUPLICATES _test_support_sources)
        target_sources("${_test_name}" PRIVATE ${_test_support_sources})
    endif()

    if(OpenMP_CXX_FOUND)
        target_link_libraries("${_test_name}" PRIVATE OpenMP::OpenMP_CXX)
    endif()

    # Apply ISA flags only to tests whose filenames indicate the required SIMD
    # path. These binaries must only be run on CPUs supporting that ISA.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options("${_test_name}" PRIVATE -O2)
        if(_test_name MATCHES "avx512|zmm|512")
            target_compile_options("${_test_name}" PRIVATE -mavx512f -mavx512dq -mavx512bw -mavx512vl -mfma)
        elseif(_test_name MATCHES "avx2")
            target_compile_options("${_test_name}" PRIVATE -mavx2 -mfma)
        elseif(_test_name MATCHES "avx|ymm")
            target_compile_options("${_test_name}" PRIVATE -mavx -mfma)
        elseif(_test_name MATCHES "sse|xmm")
            target_compile_options("${_test_name}" PRIVATE -msse4.2)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Intel classic (icpc)
        target_compile_options("${_test_name}" PRIVATE -O2)
        if(_test_name MATCHES "avx512|zmm|512")
            target_compile_options("${_test_name}" PRIVATE -xCORE-AVX512 -mavx512f)
        elseif(_test_name MATCHES "avx2")
            target_compile_options("${_test_name}" PRIVATE -xCORE-AVX2)
        elseif(_test_name MATCHES "avx|ymm")
            target_compile_options("${_test_name}" PRIVATE -xAVX)
        elseif(_test_name MATCHES "sse|xmm")
            target_compile_options("${_test_name}" PRIVATE -xSSE4.2)
        endif()
    else()
        message(WARNING "No SIMD flags selected for ${_test_name}: unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()

    if(RF_EMT_REGISTER_TESTS)
        include(CTest)
        add_test(NAME "${_test_name}" COMMAND "$<TARGET_FILE:${_test_name}>")
    endif()
endforeach()

message(STATUS "Configured ${CMAKE_PROJECT_NAME} with ${CMAKE_CXX_COMPILER_ID}: ${CMAKE_CXX_COMPILER}")
message(STATUS "Discovered ${_rf_emt_test_sources} unit-test source files")
