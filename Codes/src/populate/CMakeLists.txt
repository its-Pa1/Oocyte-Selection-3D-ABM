cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
# ============================================================================
# Find packages
# ============================================================================
find_package(OpenGL REQUIRED)

# Manually find freeglut installed by Homebrew
find_library(FREEGLUT_LIBRARY glut PATHS /opt/homebrew/lib)
include_directories(/opt/homebrew/include)

# ============================================================================
# Add target populate add its dependencies
# ============================================================================
add_executable(populate populate.cpp)
# We use C++11
target_compile_options(populate PRIVATE "-std=c++11")

# ============================================================================
# Use OpenGL and freeglut
# ============================================================================
target_link_libraries(populate ${OPENGL_LIBRARIES} ${FREEGLUT_LIBRARY})

# ============================================================================
# Set behaviour on make install
# ============================================================================
#install(TARGETS populate DESTINATION ${CMAKE_INSTALL_BINDIR}/simuscale)