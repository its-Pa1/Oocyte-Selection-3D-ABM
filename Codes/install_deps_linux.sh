#!/bin/bash
# install_deps_linux.sh — Install all dependencies for building simuscale view on Linux

set -e  # exit on any error

echo "=================================================="
echo " SiMuScale View — Dependency Installer (Linux)"
echo "=================================================="

# ── Check we are on Linux ──────────────────────────────────────────────────
if [[ "$(uname)" != "Linux" ]]; then
    echo "❌ This script is for Linux only. Exiting."
    exit 1
fi

# ── Check for apt ──────────────────────────────────────────────────────────
if ! command -v apt-get &> /dev/null; then
    echo "❌ apt-get not found. Are you on a Debian/Ubuntu system?"
    exit 1
fi

# ── System packages ────────────────────────────────────────────────────────
echo ""
echo "📦 Installing system packages..."
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    freeglut3-dev \
    libpng-dev \
    libopencv-dev \
    libgsl-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev

# ── Verify key libraries are findable ─────────────────────────────────────
echo ""
echo "🔍 Verifying installations..."

MISSING=0

check_lib() {
    if ldconfig -p | grep -q "$1"; then
        echo "  ✅ $1 found"
    else
        echo "  ❌ $1 NOT found"
        MISSING=$((MISSING + 1))
    fi
}

check_lib "libGL.so"
check_lib "libGLU.so"
check_lib "libglut.so"
check_lib "libpng"
check_lib "libopencv"
check_lib "libgsl.so"
check_lib "libgslcblas.so"

# ── OpenCV cmake path ──────────────────────────────────────────────────────
OPENCV_CMAKE_DIR=$(find /usr/lib -name "OpenCVConfig.cmake" 2>/dev/null | head -1 | xargs dirname 2>/dev/null)
if [[ -n "$OPENCV_CMAKE_DIR" ]]; then
    echo "  ✅ OpenCV cmake config found at: $OPENCV_CMAKE_DIR"
else
    echo "  ⚠️  OpenCV cmake config not found — you may need -DOpenCV_DIR=... manually"
fi

# ── Summary ────────────────────────────────────────────────────────────────
echo ""
if [[ $MISSING -eq 0 ]]; then
    echo "✅ All dependencies installed successfully!"
    echo ""
    echo "To build, run:"
    echo "  cd src/view"
    echo "  mkdir -p build && cd build"
    if [[ -n "$OPENCV_CMAKE_DIR" ]]; then
        echo "  cmake .. -DOpenCV_DIR=${OPENCV_CMAKE_DIR}"
    else
        echo "  cmake .."
    fi
    echo "  make view -j\$(nproc)"
else
    echo "⚠️  $MISSING library/libraries not found after install."
    echo "   Try running: sudo ldconfig"
    echo "   Then re-run this script."
fi
