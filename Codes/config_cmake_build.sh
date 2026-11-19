#!/bin/bash
# configure.sh - Complete configuration and build script

set -e  # Exit on any error

# Create build directory if it doesn't exist
if [ ! -d "build" ]; then
    echo "Creating build directory..."
    mkdir build
fi

# Change to build directory
cd build

echo "Configuring project with CMake..."
cmake ..

echo "Building project..."
make # Use all available CPU cores

echo "Build completed successfully, calling from script!"