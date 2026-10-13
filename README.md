# 2D Kinetic Alpha Shapes

A C++20 library for 2D kinetic alpha shapes, including kinetic Delaunay triangulation and related algorithms.

## Prerequisites

- **CMake** 3.14 or higher
- **C++20 compatible compiler** (MSVC recommended for Windows)
- **vcpkg** - C++ package manager
- **Visual Studio 2022** (includes MSVC compiler and build tools)
- **Ninja** - Required for the default CMake preset (`vcpkg-x64`); install via winget/choco/scoop (see below)

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

Or set it via System Properties â†’ Environment Variables.

### 2. Install Dependencies

Install required packages using vcpkg:

```powershell
vcpkg install glm:x64-windows
vcpkg install cgal:x64-windows
```

### 3. Build the Project

**Using CMake Presets (Recommended):**

Use **Developer PowerShell for VS 2022** (or run `vcvars64.bat`) so `cl.exe` is on `PATH`, then:

```powershell
# Configure (Ninja + MSVC via vcpkg toolchain; Release by default)
cmake --preset vcpkg-x64

# Build the library and tests (single-config: no --config flag)
cmake --build build --target kinDS kinDS-tests

# Debug build
cmake --preset vcpkg-x64-debug
cmake --build build --target kinDS kinDS-tests
```

The demo executable will be built to `build/bin/kinDS-demo.exe` (or `build/bin/kinDS-demo` on Unix).

**Or manually (Ninja):**

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="cmake/vcpkg-ninja-msvc.toolchain.cmake"
cmake --build build --target kinDS kinDS-tests
```

**CGAL + compiler note:** vcpkg's CGAL 5.6.x does not build with **Clang** on Windows (invalid `operator bool` in CGAL headers). The Ninja presets use `cmake/vcpkg-ninja-msvc.toolchain.cmake` to force **MSVC** (`cl.exe`). If you still see Clang in compile commands, delete `build/` and reconfigure.

**Visual Studio generator (optional, multi-config):**

```powershell
cmake --preset vcpkg-x64-vs
cmake --build build-vs --config Release --target kinDS kinDS-tests
```

**Install Ninja:**

- **Using winget:** `winget install Ninja-build.Ninja`
- **Using Chocolatey:** `choco install ninja`
- **Using Scoop:** `scoop install ninja`
- **Manual:** Download from [ninja-build.org](https://ninja-build.org/)

**In VS Code / Cursor:**

The repo includes `.vscode/settings.json` so **Delete Cache and Reconfigure** uses the `vcpkg-x64` preset by default (Ninja + MSVC + vcpkg). Debug launch is wired to `build/bin/kinDS-demo.exe` (Ninja's flat output path). Switch the status-bar configure preset to `vcpkg-x64-debug` when you need a Debug build with symbols.

1. Open the kinDS folder as the workspace root.
2. Right-click `CMakeLists.txt` → **Delete Cache and Reconfigure**.
3. Select **kinDS-demo** as the launch target, then **Debug** (or F5).

**In Visual Studio 2022:**

Enable **CMake Presets** under **Tools → Options → CMake → General**, then reopen the folder.

## Project Structure

- `kinDS/` - Main library source code (headers and implementation)
- `eigen/` - Eigen library (included as subdirectory)
- `main.cpp` - Demo application showing library usage
- `CMakePresets.json` - CMake configuration presets
- `CMakeLists.txt` - Main CMake build configuration

## Using as a Submodule

This project can be used as a Git submodule in other CMake projects:

1. **Add as a submodule:**
   ```bash
   git submodule add <repository-url> path/to/kinDS
   git submodule update --init path/to/kinDS
   ```

2. **In your CMakeLists.txt:**
   ```cmake
   # Add the kinDS subdirectory
   add_subdirectory(path/to/kinDS)
   
   # Link to the library in your target
   target_link_libraries(your-target PRIVATE kinDS)
   ```

3. **In your code:**
   ```cpp
   #include "kinDS/KineticDelaunay.hpp"
   #include "kinDS/Polynomial.hpp"
   // ... use kinDS namespace
   ```

**Note on Eigen dependency:**
- If your project already provides Eigen (via `add_subdirectory` or `find_package`), kinDS will automatically use your Eigen installation
- If Eigen is not found, kinDS will use its bundled Eigen submodule (make sure to initialize it with `git submodule update --init path/to/kinDS`)
- The library automatically handles this detection, so no additional configuration is needed

The library will automatically handle its dependencies (glm, CGAL, Eigen) when used as a submodule.

## Notes

- The project is structured as a library (`kinDS` target) with an optional demo executable (`kinDS-demo`)
- The demo executable is only built when this is the main project (not when used as a submodule)
- The project uses `CMakePresets.json` to automatically configure vcpkg integration
- The default preset (`vcpkg-x64`) uses **Ninja** with the MSVC toolchain (Developer shell / `vcvars64`) and vcpkg; use `vcpkg-x64-vs` for the Visual Studio multi-config generator
- If `VCPKG_ROOT` is not set, update `CMakePresets.json` with your vcpkg path
- Set `VCPKG_ROOT` before configuring; update the vcpkg path in `CMakePresets.json` if yours is not `C:/src/vcpkg` (or copy `CMakeUserPresets.json.example` to `CMakeUserPresets.json` and override `VCPKG_ROOT` there)
- CGAL is optional - the library will build without it, but some features will be disabled
