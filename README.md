# xmath: High-Performance C++ 20 Math Library

**xmath** is a modern, header-only C++ library designed for efficient 3D math operations in games, graphics, simulations, and scientific applications. Built for C++20, it offers vectors, matrices, quaternions, angles, trigonometry, and general numeric functions with a focus on performance, type safety, and ease of use. Whether you're building a game engine or need reliable linear algebra, xmath delivers without bloat or dependencies.

## Key Features
- **Vectors & Matrices**: 2D/3D/4D vectors (`fvec2`, `fvec3`, `fvec4`) and 3x3/4x4 matrices (`fmat3`, `fmat4`) with SIMD/CPU variants for speed vs. memory trade-offs.
- **Quaternions**: `fquat` for gimbal-free rotations, with interpolation (Lerp, Slerp, Squad) and conversions to/from Euler/matrix.
- **Angles & Trig**: Type-safe `radian` and `degree` wrappers, Euler angles (`radian3`), and full trig suite (Sin, Cos, Tan, etc.) with normalization and minimal difference utils.
- **General Math**: Exponential, log, power, rounding, clamping, lerp, quadratic solver, and more—all templated and constexpr where possible.
- **Strong Typing**: CRTP-based wrappers prevent unit errors (e.g., mixing radians/degrees).
- **Performance Optimizations**: Optional SSE SIMD, mutable chaining for fluent APIs, inline/noexcept functions, and debug assertions for finite values.
- **Graphics Helpers**: Look-at, billboard, perspective/ortho projections, from-to rotations, vector transforms.
- **No Dependencies**: Pure C++20, header-only—drop in and go.
- **MIT License**: Free to use, modify, and distribute.
- **doubles, floats, and half-floats**: WIP

## Why Choose xmath?
xmath stands out for developers seeking a lightweight, powerful math toolkit:

- **Zero Overhead, No Dependencies**: Header-only design means no build setup or external libs—just include and compile. Perfect for cross-platform projects.
- **Performance Tailored**: SIMD variants accelerate hot paths (e.g., vector ops ~4x faster on SSE), while CPU modes save memory. Toggle per-type for flexibility.
- **Type Safety & Usability**: Strong typing (e.g., `radian` vs. raw float) prevents bugs; mutable chaining (e.g., `mat.Translate(t).Rotate(q)`) makes code readable.
- **Game/Graphics Focus**: Built-in utilities like billboard matrices, look-at, and safe normalization make it ideal for engines—fewer custom hacks needed.
- **Modern C++**: Leverages C++20 features (constexpr, concepts) for compile-time safety and efficiency.
- **Comprehensive Documentation**: Each class/header has detailed docs with examples, tables, and tutorials (generated for clarity).
- **MIT Licensed**: Open-source freedom without restrictions—use in commercial projects hassle-free.
- **Unit-Tested**: Test functions to make sure everything is working...
- **Sub-headers**: To only get what you need.

If you're tired of bloated libraries or raw float errors, xmath streamlines your math code while boosting performance. Join developers building faster, safer apps—try it today!

## Comparison to Other Libraries
xmath bridges the gap between graphics-focused and general-purpose math libs, offering unique balances:

- **vs. GLM (OpenGL Mathematics)**: GLM is excellent for GLSL-like graphics math (header-only, no deps, SIMD in some funcs). xmath adds SIMD/CPU toggles per type, strong typing for angles (e.g., `radian`), mutable chaining, and extras like billboard/look-at—from GLM's search results, xmath is more customizable for CPU-specific needs while matching graphics utilities.
- **vs. Eigen**: Eigen shines in linear algebra (dense/sparse matrices, solvers, broad CPU support). xmath is lighter (no solvers, focused on 3D transforms), header-only with no deps, and game-optimized (e.g., safe quaternion interp, projection helpers)—Eigen's features suit scientific computing, but xmath is faster to integrate for games without overhead.

xmath's unique edge: Hybrid SIMD/CPU, type-safe angles, and chaining APIs make it versatile yet simple, outperforming in targeted use cases like real-time rendering.

## Getting Started
1. **Include Headers**: Copy the `xmath` folder to your project and `#include "xmath.h"` (or specific headers like `xmath_fvector.h`).
2. **Basic Usage**:
   ```cpp
   #include "xmath.h"

   int main() {
       xmath::fvec3 pos(1, 2, 3);
       xmath::fmat4 mat = xmath::fmat4::fromTranslation(pos).Rotate(xmath::fquat::fromAxisAngle(xmath::fvec3::fromUnitX(), 90.0_xdeg));
       xmath::fvec3 transformed = mat * xmath::fvec3(0,1,0);

       return 0;
   }
   ```
3. **Build**: Compile with C++20 enabled (e.g., `-std=c++20`). No linking needed.
4. **Examples & Docs**: Check headers for per-class docs; examples in tutorials above.

Contribute or report issues on GitHub (link if available).

