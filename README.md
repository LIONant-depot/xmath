# xmath: Fast C++ Math Library for Games & Graphics

**xmath** is a simple, header-only C++ library for 3D math in games, graphics, and simulations. Made for C++20, it includes vectors, matrices, quaternions, angles, trig functions, and more. Focuses on speed, safety, and easy use—no extras or setup needed.

## Key Features
- **Vectors & Matrices**: 2D/3D/4D vectors (`fvec2`, `fvec3`, `fvec4`) and 3x3/4x4 matrices (`fmat3`, `fmat4`). Choose SIMD or CPU mode for speed or small size.
- **Quaternions**: `fquat` for smooth rotations without gimbal lock. Includes blending (Lerp, Slerp) and easy swaps with Euler angles or matrices.
- **Angles & Trig**: Safe `radian` and `degree` types, Euler angles (`radian3`), and trig tools (Sin, Cos, etc.) with wrapping and shortest-path helpers.
- **General Math**: Exp, log, power, round, clamp, lerp, quadratic solver—all fast and ready for templates.
- **Safe Types**: Stops mistakes like mixing radians and degrees.
- **Speed Boosts**: Optional SSE SIMD, chainable methods (like `mat.Translate(t).Rotate(q)`), inline functions.
- **Graphics Tools**: Look-at, billboard, perspective/ortho views, from-to rotations.
- **No Dependencies**: Just C++20 headers—include and go.
- **MIT License**: Free for any use.
- **Extras (WIP)**: Support for double, float, half-float.
- **Tested & Modular**: Unit tests ensure reliability; sub-headers let you pick what you need.

## Why Use xmath?
xmath helps you write better code faster:

- **No Setup Hassle**: Header-only—no builds or external libs. Works anywhere C++20 runs.
- **Custom Speed**: SIMD for quick math, CPU for low memory. Switch as needed.
- **Fewer Bugs**: Types like `radian` catch errors early; chaining makes code clean.
- **Game-Ready**: Tools for cameras, rotations, projections—perfect for engines.
- **Modern & Efficient**: C++20 features like constexpr for safe, fast code.
- **Great Docs**: Each part has clear guides, examples, and tables.
- **Open & Free**: MIT license means no limits—use it commercially.

Tired of big libraries or float mix-ups? xmath keeps things simple and quick. Start using it and build better apps!

## How xmath Compares
xmath fits between graphics and general math libs, with smart trade-offs:

- **vs. GLM**: GLM is great for GLSL-style math (header-only, some SIMD). xmath adds type safety for angles, chaining, and SIMD/CPU choices—more flexible for mixed needs.
- **vs. Eigen**: Eigen is strong for advanced algebra (matrices, solvers). xmath is lighter, no deps, and tuned for games (safe rotations, projections)—easier for quick integrates.

xmath shines with its mix of speed, safety, and simplicity—ideal for real-time work.

## Get Started
1. **Add to Project**: Copy the `xmath` folder and include headers like `#include "xmath.h"`.
2. **Simple Example**:
   ```cpp
   #include "xmath.h"

   int main() {
       xmath::fvec3 pos(1, 2, 3);
       xmath::fmat4 mat = xmath::fmat4::fromTranslation(pos).Rotate(xmath::fquat::fromAxisAngle(xmath::fvec3::fromUnitX(), 90.0_xdeg));
       xmath::fvec3 transformed = mat * xmath::fvec3(0, 1, 0);

       return 0;
   }
   ```
3. **Build**: Use C++20 (e.g., `-std=c++20`). No linking required.
4. **More Help**: Headers have docs and examples.

Join the fun—try xmath now and speed up your project! 🚀
