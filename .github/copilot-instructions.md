# Gkeyll AI Agent Instructions

## Project Overview

**Gkeyll** is a multi-scale, multi-physics simulation framework for plasma physics, space physics, general relativity, and high-energy astrophysics. Written primarily in **C (C99/C11)** with minimal external dependencies, it features a modular hierarchical architecture and an optional Lua scripting layer for simulation input. Always generate code that is machine-precision correct, numerically stable, and optimized for clarity and maintainability before micro‑optimizations. Prefer small, well‑named functions, explicit memory management, and exhaustive error checking; preserve existing patterns, abstractions, and coding style when extending the code.

**Key Philosophy**: Minimize dependencies, keep code simple and dependency-free, follow a "Minimalist Approach to Software" (see https://www.ammar-hakim.org/sj/pn/pn0/pn0-minimalism.html).

## Architecture Layers

The codebase is organized in strict hierarchical layers, each depending only on the layer below:

```
core → moments → vlasov → gyrokinetic
                      ↘ pkpm
```

- **`core/`**: Essential infrastructure (grids, parallelism, I/O, DG operations, Lua tools)
- **`moments/`**: Finite-volume solvers for hyperbolic PDEs (multi-fluid, Maxwell, Einstein equations)
- **`vlasov/`**: Modal DG solvers for Vlasov-Maxwell/Poisson and Hamiltonian systems
- **`gyrokinetic/`**: Modal DG solvers for gyrokinetic equations with tokamak/mirror geometries
- **`pkpm/`**: Parallel-Kinetic-Perpendicular-Moment system solvers
- **`gkeyll/`**: Lua scripting executable and main entry point

Each layer can be built independently: `make core`, `make moments`, `make vlasov`, `make gyrokinetic`, `make pkpm`.

## Critical Naming Conventions

### Struct Naming Pattern
All public structs use the `gkyl_` prefix:
```c
struct gkyl_array
struct gkyl_gyrokinetic_app
struct gkyl_rect_grid
```

### Private Implementation Pattern
Private implementation details use `_priv.h` headers:
```c
#include <gkyl_array_rio_priv.h>      // Private implementation
#include <gkyl_gyrokinetic_priv.h>    // Private app details
```

**Rule**: Never expose private headers in public API. Only include `_priv.h` in `.c` implementation files.

### Directory Structure Per Layer
```
<layer>/
  ├── zero/          # Core implementation (C files)
  ├── apps/          # Application-level code
  ├── unit/          # Unit tests
  ├── creg/          # C regression tests
  ├── luareg/        # Lua regression tests
  └── Makefile-<layer>
```

## Build System & Workflows

### Configuration & Building
```bash
# 1. Install dependencies (optional if cluster has them)
cd install-deps && ./mkdeps.sh --build-openblas=yes --build-superlu=yes

# 2. Configure (use machine files when available)
./machines/configure.<machine>.sh
# OR manually: ./configure --prefix=$HOME/gkylsoft --app=gyrokinetic

# 3. Compile
make install -j 8              # Full build
make gyrokinetic -j 8         # Single layer
make gkeyll -j 8              # Lua executable only
```

### Key Build Flags
- `BUILD_APP=<layer>`: Build specific layer (moments, vlasov, gyrokinetic, pkpm)
- `CC=nvcc`: Enable CUDA support (changes `BUILD_DIR` to `cuda-build/`)
- `CUDA_ARCH=70`: V100 GPUs (use 80 for A100)
- `USE_MPI=1`: Enable MPI support
- `PREFIX`: Installation directory (default: `$HOME/gkylsoft`)

### Testing Commands
```bash
# Build & run unit tests
make core-unit -j 8
./build/core/unit/ctest_array

# Build & run regression tests  
make gyrokinetic-regression -j 8
./build/gyrokinetic/creg/rt_gk_mirror_boltz_elc_1x2v_p1 -s1

# Run all unit tests
make check -j 8

# Memory leak detection (CRITICAL before commits)
valgrind --leak-check=full ./build/core/unit/ctest_array
compute-sanitizer --tool memcheck --leak-check full <executable>  # GPU
```

**All code MUST be valgrind-clean and compute-sanitizer-clean.**

## Running Simulations

### Lua Input Files (Preferred)
```bash
gkylsoft/gkeyll/bin/gkeyll simulation.lua
```

Lua apps follow this pattern:
```lua
local Gyrokinetic = require "App.Gyrokinetic"

gkApp = Gyrokinetic.App.new {
  tEnd = 1e-6,
  nFrame = 10,
  -- ... configuration
}

gkApp:run()
```

### C Input Files (For Advanced Features)
```bash
# Compile C input file
cp gkylsoft/gkeyll/Makefile .
make  # Creates executable from <name>.c

# Run (or use cgkeyll for CPU-only)
./<executable>
gkylsoft/gkeyll/bin/cgkeyll input.c  # Compile + run
```

C apps use runtime flags: `-h` for help, `-s<N>` for N steps, `-g` for GPU.

## Code Style & Standards

### Memory Management
- Use `gkyl_malloc()` / `gkyl_free()` from `<gkyl_alloc.h>`
- Reference counting pattern: `_acquire()` / `_release()` methods
- **Zero tolerance** for memory leaks - always run valgrind

### Function Naming
```c
// Public API: gkyl_<module>_<action>
struct gkyl_array* gkyl_array_new(...)
void gkyl_gyrokinetic_app_write(...)

// Updaters follow: gkyl_<category>_updater_<name>
struct gkyl_dg_updater_collisions* gkyl_dg_updater_lbo_gyrokinetic_new(...)
```

### Header Guards
Use `#pragma once` (project standard).

### Documentation
Reference https://en.cppreference.com/w/ for C standards. Follow existing patterns in similar files.

## Common Pitfalls

1. **Don't use `ker/` generated kernels directly** - these are auto-generated. Modify templates instead.

2. **Layer dependencies are strict** - never include higher-layer headers in lower layers:
   ```c
   // WRONG in core/: #include <gkyl_vlasov.h>
   // RIGHT in vlasov/: #include <gkyl_array.h> from core
   ```

3. **Configure before building** - machine files set critical paths:
   ```bash
   # Always check configure script for required modules on clusters
   module load cudatoolkit/12.4  # From configure.stellar-amd.sh
   ```

4. **GPU vs CPU builds are different**:
   - `make CC=cc`: builds in `build/`
   - `make CC=nvcc`: builds in `cuda-build/`

5. **Lua function context** - When adding Lua-callable functions, use `lua_func_ctx` pattern (see `moment_lw.c`, `vlasov_lw.c` for examples).

## File Output Patterns

Simulations output `.gkyl` binary files: `<sim_name>-<field>_<frame>.gkyl`

Example: `gk_mirror-ion_M0_0.gkyl` (ion density moment at frame 0)

## GPU Development

When `CC=nvcc`:
- `.cu` files compiled with CUDA
- Use `-arch=sm_${CUDA_ARCH}` for compute capability
- Include CUDA math libs: `-lcublas -lcusparse -lcusolver`
- Multi-GPU: decompose with `-e <dir>` flag (direction 0, 1, or 2)

## Testing Philosophy

- **Unit tests** (`unit/`): Test individual functions/modules
- **C regression tests** (`creg/`): Full C-based simulations  
- **Lua regression tests** (`luareg/`): Lua-based simulation tests

Regression tests support `-s1` (1 step), `-x4 -y4 -z4` (grid resolution), etc.

## Key Abbreviations

- `gkyl` = Gkeyll (common in code despite project name)
- `dg` = Discontinuous Galerkin
- `mom` = moments
- `lbo` = Lenard-Bernstein Operator
- `creg` = C regression tests
- `luareg` = Lua regression tests

## Getting Help

- Documentation: http://gkeyll.rtfd.io
- Build issues: Check `machines/configure.<machine>.sh` for examples
- Memory errors: **Always** run valgrind before submitting PRs

// ...existing code...

## Key Abbreviations

- `gkyl` = Gkeyll (common in code despite project name)
- `dg` = Discontinuous Galerkin
- `mom` = moments
- `lbo` = Lenard-Bernstein Operator
- `creg` = C regression tests
- `luareg` = Lua regression tests

## Architectural Design Principles

Do not use forward declarations. All local functions to files should not be declared in headers or only in the private header.

### App vs Updater Responsibilities

**Apps (`apps/`)** should be thin scripting layers that orchestrate updaters:
- No intensive algorithms - delegate to `zero/` updaters
- Think of apps as "glue code" connecting updaters
- Do not update or access members of updaters directly

**Updaters (`zero/`)** contain all computational logic:
- Self-contained modules that own their data
- Compute derived quantities internally during construction
- Provide getter functions for data that apps need
- Handle all memory allocation for their internal state

```c
// WRONG - App doing computation
void app_init(...) {
  // Intensive computation in app
  for (int i = 0; i < n; i++) {
    compute_expensive_thing(...);
  }
}

// RIGHT - App delegates to updater
void app_init(...) {
  // Updater does computation internally
  struct gkyl_my_updater *up = gkyl_my_updater_new(inp);
  // App retrieves results via getter
  const double *result = gkyl_my_updater_get_result(up);
}
```

### Static Helper Functions and Function Pointers

Helper functions are defined as `static` at the top of the file, then assigned to function pointers in structs during initialization. This enables runtime polymorphism without virtual functions.

**Core Principle: Move Conditionals to Initialization**

Evaluate `if` statements once during `_new()`, not repeatedly in hot paths. When runtime behavior varies based on configuration, assign different function pointers at init time rather than checking conditions on every call.

```c
// WRONG - Conditional evaluated every call (potentially millions of times)
void my_struct_advance(struct my_struct *s, double tm)
{
  if (s->feature_enabled) {
    do_feature_work(s, tm);
  }
  // ... rest of advance
}

// CORRECT - Conditional evaluated once at init, function pointer called directly
static void advance_with_feature(struct my_struct *s, double tm)
{
  do_feature_work(s, tm);
  // ... rest of advance
}

static void advance_without_feature(struct my_struct *s, double tm)
{
  // ... rest of advance (no feature work)
}

// In _new(), decide once which function to use
struct my_struct* my_struct_new(struct my_struct_inp *inp)
{
  struct my_struct *s = gkyl_malloc(sizeof(*s));
  
  if (inp->feature_enabled) {
    s->advance_func = advance_with_feature;
  } else {
    s->advance_func = advance_without_feature;
  }
  return s;
}

// Caller just invokes the function pointer - no conditional
void my_struct_advance(struct my_struct *s, double tm)
{
  s->advance_func(s, tm);
}
```
