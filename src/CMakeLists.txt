cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(RF-EMT 
    VERSION 1.0.0 
    LANGUAGES CXX
)

# ---------------------------------------------------------------------------
# 1. User Options
# ---------------------------------------------------------------------------
option(RF_EMT_ENABLE_CUDA "Enable NVIDIA CUDA GPGPU Acceleration" ON)
option(RF_EMT_ENABLE_AVX512 "Enable Intel AVX-512 SIMD Vectorization" ON)
option(RF_EMT_BUILD_EXAMPLES "Build standalone examples and tools" ON)

# Force modern C++ standards
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# 2. CUDA Language Initialization (Dynamic)
# ---------------------------------------------------------------------------
if(RF_EMT_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        message(STATUS "CUDA Support Enabled. Compiler: ${CMAKE_CUDA_COMPILER}")
        
        # Target modern compute architectures (adjust based on deployment hardware)
        # 80 = Ampere (A100/RTX 30xx), 89 = Ada Lovelace (RTX 40xx), 90 = Hopper
        set(CMAKE_CUDA_ARCHITECTURES "80;89;90")
        
        # Ensure proper separation of device code tracking
        set(CMAKE_CUDA_SEPARATE_COMPILATION ON)
    else()
        message(WARNING "CUDA compiler not found. Falling back to CPU-only SIMD implementation.")
        set(RF_EMT_ENABLE_CUDA OFF)
    endif()
endif()

# ---------------------------------------------------------------------------
# 3. Compiler-Specific Flags (SIMD / Intel Intrinsics Alignment)
# ---------------------------------------------------------------------------
if(MSVC)
    add_compile_options(/W4 /WX /Zc:__cplusplus)
    if(RF_EMT_ENABLE_AVX512)
        add_compile_options(/arch:AVX512)
    else()
        add_compile_options(/arch:AVX2)
    endif()
else() # GCC / Clang
    add_compile_options(-Wall -Wextra -Wpedantic -O3)
    if(RF_EMT_ENABLE_AVX512)
        # Enable full array of AVX512 instructions required for analytical RCS/Antenna modeling
        add_compile_options(-mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl -mfma)
    else()
        add_compile_options(-mavx2 -mfma -msse4.2)
    endif()
endif()

# ---------------------------------------------------------------------------
# 4. Define the Primary Computational Backend Target
# ---------------------------------------------------------------------------
# Collect files natively. Group logic conditionally depending on compilation path
file(GLOB_RECURSE RF_EMT_HOST_SOURCES 
    "src/*.cpp" 
    "src/*.hpp" 
    "src/*.h"
)

if(RF_EMT_ENABLE_CUDA)
    file(GLOB_RECURSE RF_EMT_CUDA_SOURCES 
        "src/*.cu" 
        "src/*.cuh"
    )
endif()

# Define core toolkit library
add_library(rf_emt_backend STATIC 
    ${RF_EMT_HOST_SOURCES}
    ${RF_EMT_CUDA_SOURCES}
)

# Set Target Include Directories (Abides by modern "Target-based" CMake rules)
target_include_directories(rf_emt_backend 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kernels
)

# ---------------------------------------------------------------------------
# 5. Core Feature Macro Definitions
# ---------------------------------------------------------------------------
if(RF_EMT_ENABLE_CUDA)
    target_compile_definitions(rf_emt_backend PUBLIC RF_EMT_USE_CUDA)
endif()

if(RF_EMT_ENABLE_AVX512)
    target_compile_definitions(rf_emt_backend PUBLIC RF_EMT_USE_AVX512)
endif()

# ---------------------------------------------------------------------------
# 6. Linking Third-Party Modules (Standard Math/Threading)
# ---------------------------------------------------------------------------
find_package(Threads REQUIRED)
target_link_libraries(rf_emt_backend PUBLIC Threads::Threads)

# ---------------------------------------------------------------------------
# 7. Executables / Prototyping Drivers
# ---------------------------------------------------------------------------
if(RF_EMT_BUILD_EXAMPLES)
    # Checks for an explicit main entry point simulation script
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/apps/main.cpp")
        add_executable(rf_emt_sim apps/main.cpp)
        target_link_libraries(rf_emt_sim PRIVATE rf_emt_backend)
        message(STATUS "Created executable: rf_emt_sim")
    endif()
endif()
