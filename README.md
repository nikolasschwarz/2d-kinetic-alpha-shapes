# Kinetic Data Structures

## Prerequisites

- **CMake** 3.14 or higher
- **C++20 compatible compiler** (MSVC recommended for Windows)
- **vcpkg** - C++ package manager
- **Visual Studio 2022** (includes MSVC compiler and build tools)
- **Ninja** (optional) - Only needed if using the `vcpkg-x64-ninja` preset

## Setup Instructions

### 1. Set Up vcpkg

If you haven't already installed vcpkg:

```powershell
# Clone vcpkg
git clone https://github.com/Microsoft/vcpkg.git path\to\vcpkg

# Bootstrap vcpkg
cd path\to\vcpkg
.\bootstrap-vcpkg.bat
```

**Set the VCPKG_ROOT environment variable:**
```powershell
# PowerShell (current session)
$env:VCPKG_ROOT = "path\to\vcpkg"

# To make it permanent (run as Administrator)
[System.Environment]::SetEnvironmentVariable("VCPKG_ROOT", "path\to\vcpkg", "Machine")
```

Or set it via System Properties → Environment Variables.

### 2. Install Dependencies

Install required packages using vcpkg:

```powershell
vcpkg install glm:x64-windows
vcpkg install cgal:x64-windows
```

### 3. Build the Project

**Using CMake Presets (Recommended):**

```powershell
# Configure (uses Visual Studio 2022 generator with MSVC)
cmake --preset vcpkg-x64

# Build
cmake --build build --config Debug
```

**Or manually:**

```powershell
# Configure
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -G "Visual Studio 17 2022" -A x64

# Build
cmake --build build --config Debug
```

**Using Ninja (Optional):**

If you prefer Ninja for faster builds, install it first (see below), then use:

```powershell
# Configure with Ninja preset
cmake --preset vcpkg-x64-ninja

# Build
cmake --build build
```

**Install Ninja (only if using Ninja preset):**

- **Using winget:** `winget install Ninja-build.Ninja`
- **Using Chocolatey:** `choco install ninja`
- **Using Scoop:** `scoop install ninja`
- **Manual:** Download from [ninja-build.org](https://ninja-build.org/)

**In VS Code:**
1. Open the Command Palette (Ctrl+Shift+P)
2. Select "CMake: Select Configure Preset"
3. Choose "vcpkg (x64)"
4. Build using the CMake extension

## Project Structure

- `kinDS/` - Main library source code
- `eigen/` - Eigen library (included as subdirectory)
- `main.cpp` - Application entry point
- `CMakePresets.json` - CMake configuration presets
- `CMakeLists.txt` - Main CMake build configuration

## Notes

- The project uses `CMakePresets.json` to automatically configure vcpkg integration
- The default preset (`vcpkg-x64`) uses Visual Studio 2022 generator with MSVC compiler, which matches vcpkg's CGAL build
- The preset uses the `VCPKG_ROOT` environment variable to locate vcpkg
- If `VCPKG_ROOT` is not set, update `CMakePresets.json` with your vcpkg path
- Ninja is optional and only needed if you want to use the `vcpkg-x64-ninja` preset for faster incremental builds
