#!/bin/bash
# =============================================================================
# install_deps_mac.sh — macOS dependency installer for SimuScale
# =============================================================================

set -e  # Exit immediately on error
set -o pipefail

# ── Colors for output ─────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_step()  { echo -e "\n${BLUE}▶ $1${NC}"; }
print_ok()    { echo -e "${GREEN}✔ $1${NC}"; }
print_warn()  { echo -e "${YELLOW}⚠ $1${NC}"; }
print_error() { echo -e "${RED}✖ $1${NC}"; }

echo ""
echo "=============================================="
echo "   SimuScale — macOS Dependency Installer"
echo "=============================================="
echo ""

# ── Step 1: Check for Homebrew ────────────────────────────────────────────────
print_step "Checking for Homebrew..."

if ! command -v brew &>/dev/null; then
    print_warn "Homebrew not found. Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    # Add Homebrew to PATH for both Apple Silicon and Intel Macs
    if [[ -x "/opt/homebrew/bin/brew" ]]; then
        echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> "$HOME/.zprofile"
        eval "$(/opt/homebrew/bin/brew shellenv)"
    elif [[ -x "/usr/local/bin/brew" ]]; then
        echo 'eval "$(/usr/local/bin/brew shellenv)"' >> "$HOME/.zprofile"
        eval "$(/usr/local/bin/brew shellenv)"
    fi

    print_ok "Homebrew installed successfully."
else
    print_ok "Homebrew already installed: $(brew --version | head -1)"
fi

# ── Step 2: Update Homebrew ───────────────────────────────────────────────────
print_step "Updating Homebrew..."
brew update
print_ok "Homebrew updated."

# ── Step 3: Install core build dependencies ───────────────────────────────────
print_step "Installing core build dependencies (cmake, gsl, open-mpi)..."
brew install cmake gsl open-mpi
print_ok "Core build dependencies installed."

# ── Step 4: Install image/graphics libraries ─────────────────────────────────
print_step "Installing image and graphics libraries (opencv, libpng, freeglut)..."
brew install opencv libpng freeglut
print_ok "Image and graphics libraries installed."

# ── Step 5: Verify all required tools ────────────────────────────────────────
print_step "Verifying installations..."

MISSING=()

check_tool() {
    if command -v "$1" &>/dev/null || brew list "$2" &>/dev/null 2>&1; then
        print_ok "$2 — OK"
    else
        print_error "$2 — NOT FOUND"
        MISSING+=("$2")
    fi
}

check_tool cmake     cmake
check_tool mpirun    open-mpi
brew list gsl      &>/dev/null && print_ok "gsl — OK"      || { print_error "gsl — NOT FOUND";      MISSING+=(gsl); }
brew list opencv   &>/dev/null && print_ok "opencv — OK"   || { print_error "opencv — NOT FOUND";   MISSING+=(opencv); }
brew list libpng   &>/dev/null && print_ok "libpng — OK"   || { print_error "libpng — NOT FOUND";   MISSING+=(libpng); }
brew list freeglut &>/dev/null && print_ok "freeglut — OK" || { print_error "freeglut — NOT FOUND"; MISSING+=(freeglut); }

# ── Step 6: Print OpenCV cmake path (needed for build) ───────────────────────
print_step "Detecting OpenCV CMake path..."
OPENCV_CMAKE_DIR=$(brew --prefix opencv)/lib/cmake/opencv4
if [ -d "$OPENCV_CMAKE_DIR" ]; then
    print_ok "OpenCV CMake dir: $OPENCV_CMAKE_DIR"
else
    print_warn "OpenCV CMake directory not found at expected path: $OPENCV_CMAKE_DIR"
fi

# ── Step 7: Summary ───────────────────────────────────────────────────────────
echo ""
echo "=============================================="
if [ ${#MISSING[@]} -eq 0 ]; then
    echo -e "${GREEN}✔ All dependencies installed successfully!${NC}"
    echo ""
    echo "You can now build SimuScale with:"
    echo ""
    echo "  mkdir -p build && cd build"
    echo "  cmake -DCMAKE_BUILD_TYPE=Release -DOpenCV_DIR=\$(brew --prefix opencv)/lib/cmake/opencv4 .."
    echo "  cmake --build ."
    echo ""
    echo "And build the viewer with:"
    echo ""
    echo "  cd src/view"
    echo "  mkdir -p build && cd build"
    echo "  cmake .. -DOpenCV_DIR=\$(brew --prefix opencv)/lib/cmake/opencv4"
    echo "  make view -j\$(sysctl -n hw.logicalcpu)"
else
    echo -e "${RED}✖ The following packages failed to install:${NC}"
    for pkg in "${MISSING[@]}"; do
        echo "   - $pkg"
    done
    echo ""
    echo "Try installing them manually with:"
    echo "  brew install ${MISSING[*]}"
fi
echo "=============================================="
echo ""
