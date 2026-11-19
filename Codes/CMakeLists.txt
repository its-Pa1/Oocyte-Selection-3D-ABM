# ============================================================================
# Require minimal version of cmake
# ============================================================================
cmake_minimum_required(VERSION 3.10)


# ============================================================================
# Set project name and languages
# ============================================================================
project(simuscale
        VERSION 1.2.0
        LANGUAGES CXX C)


# ============================================================================
# Set Apple architecture automatically
# ============================================================================
if(APPLE)
    execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE SYSTEM_ARCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(SYSTEM_ARCH STREQUAL "arm64")
        set(CMAKE_OSX_ARCHITECTURES "arm64")
        message(STATUS "Building for Apple Silicon (arm64)")
    else()
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
        message(STATUS "Building for Intel Mac (x86_64)")
    endif()
endif()


# ============================================================================
# Standard directories
# ============================================================================
include(GNUInstallDirs)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# ============================================================================
# C++ Standard (global)
# ============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# ============================================================================
# BUILD_VIEW option
# ON  by default on macOS and Windows (local machines)
# ON  by default on Linux desktops (DISPLAY or WAYLAND_DISPLAY is set)
# OFF by default on headless Linux servers (no display detected)
# Can always be overridden manually: cmake -DBUILD_VIEW=ON/OFF ..
# ============================================================================
if(APPLE OR WIN32)
    set(DEFAULT_BUILD_VIEW ON)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(DEFINED ENV{DISPLAY} OR DEFINED ENV{WAYLAND_DISPLAY})
        set(DEFAULT_BUILD_VIEW ON)
        message(STATUS "Display detected (Linux desktop) → BUILD_VIEW defaults ON")
    else()
        set(DEFAULT_BUILD_VIEW OFF)
        message(STATUS "No display detected (Linux server/HPC) → BUILD_VIEW defaults OFF")
    endif()
else()
    set(DEFAULT_BUILD_VIEW ON)
endif()

option(BUILD_VIEW "Build the view/GUI component" ${DEFAULT_BUILD_VIEW})
message(STATUS "BUILD_VIEW = ${BUILD_VIEW}")


# ============================================================================
# Windows (MSVC) flags
# ============================================================================
if(WIN32)
    message(STATUS "Configuring for Windows (MSVC)")

    add_compile_definitions(
        _USE_MATH_DEFINES
        _CRT_SECURE_NO_WARNINGS
        NOMINMAX
    )

    add_compile_options(/wd4244 /wd4267 /wd4305 /wd4146)
    add_compile_options(/permissive-)

    find_package(GSL REQUIRED)
endif()


# ============================================================================
# GSL Detection
# ============================================================================
if(APPLE)
    execute_process(
        COMMAND brew --prefix gsl
        OUTPUT_VARIABLE GSL_HOMEBREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE BREW_RESULT
    )

    if(BREW_RESULT EQUAL 0 AND EXISTS ${GSL_HOMEBREW_PREFIX})
        message(STATUS "Using Homebrew GSL: ${GSL_HOMEBREW_PREFIX}")

        add_library(GSL::gslcblas SHARED IMPORTED)
        set_target_properties(GSL::gslcblas PROPERTIES
            IMPORTED_LOCATION "${GSL_HOMEBREW_PREFIX}/lib/libgslcblas.dylib"
            INTERFACE_INCLUDE_DIRECTORIES "${GSL_HOMEBREW_PREFIX}/include"
        )

        add_library(GSL::gsl SHARED IMPORTED)
        set_target_properties(GSL::gsl PROPERTIES
            IMPORTED_LOCATION "${GSL_HOMEBREW_PREFIX}/lib/libgsl.dylib"
            INTERFACE_INCLUDE_DIRECTORIES "${GSL_HOMEBREW_PREFIX}/include"
            INTERFACE_LINK_LIBRARIES GSL::gslcblas
        )

    else()
        message(STATUS "Falling back to system GSL (macOS)")
        find_package(GSL REQUIRED)
    endif()

elseif(UNIX AND NOT APPLE)
    message(STATUS "Configuring for Linux")

    if(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}")
        set(CONDA_ENV $ENV{CONDA_PREFIX})
        message(STATUS "Using Conda environment: ${CONDA_ENV}")

        set(GSL_INCLUDE_DIR "${CONDA_ENV}/include")
        set(GSL_LIBRARY "${CONDA_ENV}/lib/libgsl.so")
        set(GSL_CBLAS_LIBRARY "${CONDA_ENV}/lib/libgslcblas.so")
        set(GSL_INCLUDE_DIRS "${CONDA_ENV}/include")
        set(GSL_LIBRARIES "${CONDA_ENV}/lib/libgsl.so" "${CONDA_ENV}/lib/libgslcblas.so")

        list(APPEND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CONDA_ENV}/lib")

    else()
        message(STATUS "Conda not active → using system GSL")
        find_package(GSL REQUIRED)

        if(NOT GSL_FOUND)
            message(FATAL_ERROR "GSL not found. Install with: sudo apt install libgsl-dev")
        endif()
    endif()

    if(NOT TARGET GSL::gsl)
        add_library(GSL::gslcblas SHARED IMPORTED)
        set_target_properties(GSL::gslcblas PROPERTIES
            IMPORTED_LOCATION "${GSL_CBLAS_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIRS}"
        )

        add_library(GSL::gsl SHARED IMPORTED)
        set_target_properties(GSL::gsl PROPERTIES
            IMPORTED_LOCATION "${GSL_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES GSL::gslcblas
        )
    endif()

else()
    find_package(GSL REQUIRED)
endif()


# ============================================================================
# Verify GSL
# ============================================================================
if(NOT TARGET GSL::gsl)
    message(FATAL_ERROR "GSL not found!")
else()
    message(STATUS "GSL configured successfully")
endif()


# ============================================================================
# Output directory
# ============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


# ============================================================================
# Compilation flags (Linux/macOS debug only)
# ============================================================================
if(NOT WIN32)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -O0")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -O0")
endif()


# ============================================================================
# Pass BUILD_VIEW to source code as a compile definition
# ============================================================================
if(BUILD_VIEW)
    add_compile_definitions(BUILD_VIEW=1)
else()
    add_compile_definitions(BUILD_VIEW=0)
endif()


# ============================================================================
# Subdirectories
# ============================================================================
add_subdirectory(src)


# ============================================================================
# CPack (optional packaging)
# ============================================================================
include(CPack)

