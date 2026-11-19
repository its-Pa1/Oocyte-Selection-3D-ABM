cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
project(view)

# ============================================================================
# Find packages
# ============================================================================
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(PNG REQUIRED)
find_package(OpenCV REQUIRED COMPONENTS core imgproc imgcodecs videoio)

# ============================================================================
# Add targets
# ============================================================================
add_executable(view
    view.cpp
    make_colormap.cpp make_colormap.h
    read_display_parameters.cpp read_display_parameters.h
)
add_executable(colormap colormap.c)

# ============================================================================
# Windows-specific fixes
# ============================================================================
if(WIN32)
    target_compile_options(view PRIVATE "/std:c++20" "/D_CRT_SECURE_NO_WARNINGS")
    target_compile_options(colormap PRIVATE "/D_CRT_SECURE_NO_WARNINGS")

    # getopt-win32 via vcpkg - hardcoded path as fallback
    set(VCPKG_INC "C:/vcpkg/installed/x64-windows/include")
    set(VCPKG_LIB "C:/vcpkg/installed/x64-windows/lib")

    target_include_directories(view PRIVATE ${VCPKG_INC})
    target_include_directories(colormap PRIVATE ${VCPKG_INC})

    find_library(GETOPT_LIBRARY NAMES getopt PATHS ${VCPKG_LIB} NO_DEFAULT_PATH)
    if(GETOPT_LIBRARY)
        target_link_libraries(view PRIVATE ${GETOPT_LIBRARY})
        target_link_libraries(colormap PRIVATE ${GETOPT_LIBRARY})
    else()
        message(WARNING "getopt library not found in ${VCPKG_LIB}")
    endif()
else()
    target_compile_options(view PRIVATE "-std=c++11")
endif()

# ============================================================================
# Include directories
# ============================================================================
target_include_directories(view PRIVATE
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    ${PNG_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

# ============================================================================
# Link libraries
# ============================================================================
target_link_libraries(view PRIVATE
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    PNG::PNG
    ${OpenCV_LIBRARIES}
)

# ============================================================================
# Install
# ============================================================================
install(TARGETS view DESTINATION ${CMAKE_INSTALL_BINDIR})