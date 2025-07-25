# `xmath_trigonometry.h` Documentation

## Overview

The `xmath_trigonometry.h` header in the `xmath` namespace provides type-safe angle representations and trigonometric functions optimized for mathematical computations in applications like game development, graphics, and physics simulations. It defines:
- **Angle Types**: `radian` (float radians), `degree` (float degrees), `dradian` (double radians), `ddegree` (double degrees) as strongly-typed wrappers.
- **Trigonometric Functions**: Standard functions (`Sin`, `Cos`, `Tan`, `Asin`, `Acos`, `ATan`, `ATan2`) and utilities (`SinCos`, `ModAngle`, `ModAngle2`, `MinAngleDiff`, `LerpAngle`).
- **Constants**: Convenient aliases for π and related values (`pi_v`, `pi2_v`, `pi_over2_v`, `pi_over4_v`).
- **Conversions**: Templates `DegToRad` and `RadToDeg`, plus user-defined literals (e.g., `90.0_xdeg`, `1.0_xrad`).

This header is header-only, requiring no external dependencies beyond the C++ standard library (`<numbers>` for π). It emphasizes type safety, performance, and usability, avoiding SIMD for portability but providing inline functions for efficiency.

### Key Features
- **Strong Typing**: `radian` and `degree` prevent mixing units; conversions via `getRadians()`, `getDegrees()`, or constructors.
- **Precision Variants**: Single (`radian`, `degree`) and double (`dradian`, `ddegree`) precision for flexibility.
- **Trigonometric Efficiency**: `SinCos` computes both sine and cosine in one call; all functions are `inline noexcept`.
- **Angle Normalization**: `ModAngle` ([0, 2π)), `ModAngle2` ([-π, π)), `MinAngleDiff` for shortest angular paths.
- **User-Friendly Literals**: `90.0_xdeg`, `1.0_xrad` for readable angle definitions.
- **Constants**: Predefined `radian`/`dradian` values for π, 2π, π/2, π/4, sourced from `std::numbers::pi_v`.

### Comparison to Similar Libraries
- Unlike GLM's trigonometric functions (which operate on raw floats or `glm::vecN`), `xmath` uses `radian`/`degree` for type safety, similar to `radian3`'s approach but for scalar angles.
- Compared to Eigen's math utilities, `xmath_trigonometry.h` is lighter, avoiding matrix dependencies, and provides game-focused features like `MinAngleDiff` and literals. Eigen's functions are more integrated with its matrix types.
- Engine analogs: Similar to Unity's `Mathf` (but type-safe with `radian`) or Unreal's `FMath`, but standalone and more explicit with angle types.

### Usage Notes
- **Type Safety**: Use `radian(float)` or `1.0_xrad` for radians; `degree(float)` or `90.0_xdeg` for degrees. Avoid raw floats to prevent unit errors.
- **Conversions**: Use `getDegrees()`, `getRadians()`, `DegToRad`, or `RadToDeg` for unit switching.
- **Constants**: Prefer `xmath::pi_v` (type `radian`) or `xmath::pi64_v` (`dradian`) over raw π.
- **Performance**: All functions are `inline noexcept`; `SinCos` is faster for combined sine/cosine calculations.
- **Precision**: Use `dradian`/`ddegree` for high-precision tasks (e.g., scientific computing); `radian`/`degree` for games.

### Quick Start Example
```cpp
#include "xmath_trigonometry.h"

int main() {
    xmath::radian angle = 90.0_xdeg.getRadians();  // Convert 90° to π/2 radians
    float s = xmath::Sin(angle);  // sin(π/2) = 1.0
    float c = xmath::Cos(angle);  // cos(π/2) ≈ 0.0
    xmath::radian target = xmath::pi_v;  // 180°
    xmath::radian diff = xmath::MinAngleDiff(angle, target);  // -π/2

    std::cout << "Sine: " << s << ", Diff: " << float(diff) << " rad" << std::endl;
    return 0;
}
```

## Angle Types

### `degree_t<T>` (Alias: `degree`, `ddegree`)
A templated strong-typing wrapper for degree angles, inheriting from `strong_typing_numerics_t`. Aliases:
- `degree`: `degree_t<float>` (single precision).
- `ddegree`: `degree_t<double>` (double precision).

| Member/Method | Description | Example |
|---------------|-------------|---------|
| `degree_t(T value)` | Construct from raw value (degrees). | `degree(90.0f)` |
| `degree_t(radian_t<G> rad)` | Convert from radians. | `degree(radian(1.5708f))` // ~90° |
| `radian_t<T> getRadians() const` | Convert to radians. | `90.0_xdeg.getRadians()` // π/2 |
| Arithmetic Operators | `+`, `-`, `*`, `/`, `+=`, `-=`, `*=`, `/=` | `degree a = 90.0_xdeg + 45.0_xdeg;` |
| Comparison | `<=>` (three-way comparison) | `90.0_xdeg == 90.0_xdeg` |

**Literals**:
- `_xdeg`: Creates `degree`, e.g., `90.0_xdeg`, `45_xdeg`.

### `radian_t<T>` (Alias: `radian`, `dradian`)
A templated strong-typing wrapper for radian angles. Aliases:
- `radian`: `radian_t<float>` (single precision).
- `dradian`: `radian_t<double>` (double precision).

| Member/Method | Description | Example |
|---------------|-------------|---------|
| `radian_t(T value)` | Construct from raw value (radians). | `radian(1.5708f)` // π/2 |
| `radian_t(degree_t<G> deg)` | Convert from degrees. | `radian(90.0_xdeg)` // π/2 |
| `degree_t<T> getDegrees() const` | Convert to degrees. | `xmath::pi_v.getDegrees()` // 90° |
| Arithmetic Operators | `+`, `-`, `*`, `/`, `+=`, `-=`, `*=`, `/=` | `radian a = 1.0_xrad + xmath::pi_v;` |
| Comparison | `<=>` (three-way comparison) | `xmath::pi_v == 3.1416_xrad` |

**Literals**:
- `_xrad`: Creates `radian`, e.g., `1.0_xrad`, `3.1416_xrad`.
- `_xrad64`: Creates `dradian`, e.g., `3.1416_xrad64`.

## Constants
| Constant | Type | Value | Description |
|----------|------------|-------|-------------|
| `pi_v` | `radian` | π ≈ 3.1416 | π radians (single precision). |
| `pi64_v` | `dradian` | π ≈ 3.141592653589793 | π radians (double precision). |
| `pi2_v` | `radian` | 2π ≈ 6.2832 | 2π radians. |
| `pi_over2_v` | `radian` | π/2 ≈ 1.5708 | π/2 radians (90°). |
| `pi_over4_v` | `radian` | π/4 ≈ 0.7854 | π/4 radians (45°). |

**Notes**: Sourced from `std::numbers::pi_v`; use `pi_v` for `radian`, `pi64_v` for `dradian`.

## Conversion Functions
| Template | Description | Example |
|----------|----------------|---------|
| `DegToRad<T>(T Deg)` | Convert degrees to radians. | `xmath::DegToRad(90.0f)` // π/2 |
| `RadToDeg<T>(T Rad)` | Convert radians to degrees. | `xmath::RadToDeg(xmath::pi_v)` // 90° |

**Notes**: Templated for `float` or `double`; `constexpr noexcept`.

## Trigonometric Functions
All functions are templated for `std::floating_point` types (`float`, `double`), `inline noexcept`, and operate on `radian_t<T>` inputs where applicable.

| Function | Returns | Description | Example |
|----------|---------|----------------|---------|
| `Sin(radian_t<T> x)` | `T` | Sine of angle. | `xmath::Sin(xmath::pi_over2_v)` // 1.0 |
| `Cos(radian_t<T> x)` | `T` | Cosine of angle. | `xmath::Cos(xmath::pi_v)` // -1.0 |
| `Tan(radian_t<T> x)` | `T` | Tangent of angle. | `xmath::Tan(xmath::pi_over4_v)` // 1.0 |
| `Asin(T x)` | `radian_t<T>` | Arc sine ([-1, 1] → [-π/2, π/2]). | `xmath::Asin(1.0f)` // π/2 |
| `Acos(T x)` | `radian_t<T>` | Arc cosine ([-1, 1] → [0, π]). | `xmath::Acos(-1.0f)` // π |
| `ATan(T x)` | `radian_t<T>` | Arc tangent (→ [-π/2, π/2]). | `xmath::ATan(1.0f)` // π/4 |
| `ATan2(T y, T x)` | `radian_t<T>` | Arc tangent with quadrant correction (→ [-π, π]). | `xmath::ATan2(1.0f, 1.0f)` // π/4 |
| `SinCos(radian_t<T> Angle, T& S, T& C)` | `void` | Compute sine and cosine together. | `float s, c; xmath::SinCos(xmath::pi_v, s, c);` // s=0, c=-1 |
| `ModAngle(radian_t<T> Angle)` | `radian_t<T>` | Wrap to [0, 2π). | `xmath::ModAngle(radian(7.0f))` // ~0.7168 |
| `ModAngle2(radian_t<T> Angle)` | `radian_t<T>` | Wrap to [-π, π). | `xmath::ModAngle2(radian(7.0f))` // ~-2.5664 |
| `MinAngleDiff(radian_t<T> Angle1, radian_t<T> Angle2)` | `radian_t<T>` | Shortest angle between two angles. | `xmath::MinAngleDiff(xmath::pi_v, 0.0_xrad)` // -π |
| `LerpAngle(T t, radian_t<T> Angle1, radian_t<T> Angle2)` | `radian_t<T>` | Linear interpolation with shortest path. | `xmath::LerpAngle(0.5f, 0.0_xrad, xmath::pi_v)` // π/2 |

**Notes**:
- `Sin`, `Cos`, `Tan` take `radian_t<T>`; return raw `T` (e.g., `float`).
- `Asin`, `Acos`, `ATan`, `ATan2` take raw `T` (range [-1, 1] for `Asin`/`Acos`); return `radian_t<T>`.
- `SinCos` is optimized for cases needing both values.
- `LerpAngle` uses `MinAngleDiff` internally for shortest-path interpolation.

## Tutorials and Examples

### Basic Trigonometric Calculations
Compute sine and cosine for a game character's rotation:
```cpp
xmath::radian angle = 45.0_xdeg.getRadians();  // π/4
float s = xmath::Sin(angle);  // sin(π/4) ≈ 0.707
float c = xmath::Cos(angle);  // cos(π/4) ≈ 0.707
// Faster alternative
float s2, c2;
xmath::SinCos(angle, s2, c2);  // Same results
std::cout << "Sin: " << s << ", Cos: " << c << std::endl;
```

### Angle Normalization
Wrap angles for a rotating object:
```cpp
xmath::radian angle = radian(7.0f);  // > 2π
xmath::radian wrapped1 = xmath::ModAngle(angle);  // ~0.7168 [0, 2π)
xmath::radian wrapped2 = xmath::ModAngle2(angle);  // ~-2.5664 [-π, π)
std::cout << "Wrapped [0, 2π): " << float(wrapped1) << std::endl;
```

### Minimal Angle Difference
Smooth AI character turning:
```cpp
xmath::radian current = radian(3.0f) * xmath::pi_v;  // 3π
xmath::radian target = xmath::pi_over2_v;  // π/2
xmath::radian diff = xmath::MinAngleDiff(current, target);  // -π/2
current += diff * 0.1f;  // Turn 10% per frame
current = xmath::ModAngle2(current);  // Stay in [-π, π)
```

### Angle Interpolation
Interpolate camera rotation:
```cpp
xmath::radian start = 0.0_xrad;
xmath::radian end = xmath::pi_v;
for (float t = 0; t <= 1; t += 0.1f) {
    xmath::radian angle = xmath::LerpAngle(t, start, end);  // Shortest path
    float s = xmath::Sin(angle);  // Use in transform
}
```

### Degree and Radian Conversion
Handle user input in degrees:
```cpp
xmath::degree input = 180.0_xdeg;
xmath::radian rad = input.getRadians();  // π
// Or use template
float rad2 = xmath::DegToRad(180.0f);  // π
xmath::degree deg = rad.getDegrees();  // 180°
if (deg == 180.0_xdeg) { /* Match */ }
```

### Double Precision for Scientific Computing
High-precision calculations:
```cpp
xmath::dradian angle = 3.1416_xrad64;  // π
double s = xmath::Sin(angle);  // sin(π) ≈ 0
double c = xmath::Cos(angle);  // cos(π) ≈ -1
xmath::dradian atan = xmath::ATan2(1.0, 0.0);  // π/2
std::cout << "ATan2: " << double(atan) << std::endl;
```

### Error Handling
- **Invalid Inputs**: Ensure `Asin`/`Acos` inputs are in [-1, 1]:
  ```cpp
  if (std::abs(value) <= 1.0f) {
      xmath::radian angle = xmath::Asin(value);
  }
  ```
- **Zero Handling**: `ATan2` handles `y=0`, `x=0` correctly:
  ```cpp
  xmath::radian angle = xmath::ATan2(0.0f, 0.0f);  // Defined (0 or platform-specific)
  ```
- **Wrapping**: Use `ModAngle` or `ModAngle2` post-arithmetic:
  ```cpp
  xmath::radian sum = angle1 + angle2;
  sum = xmath::ModAngle2(sum);  // Keep in [-π, π)
  ```

## Performance Tips
- **SinCos**: Use `SinCos` when both sine and cosine are needed:
  ```cpp
  float s, c;
  xmath::SinCos(angle, s, c);  // Single call, faster
  ```
- **Inlining**: All functions are `inline noexcept`, minimizing overhead.
- **Precision**: Use `radian`/`degree` for games; `dradian`/`ddegree` for high-precision tasks.
- **Constants**: Cache `xmath::pi_v` or `xmath::pi2_v` in loops to avoid repeated access.
- **Avoid Redundant Conversions**: Convert units once:
  ```cpp
  xmath::radian rad = 90.0_xdeg.getRadians();
  // Use rad repeatedly
  ```

## Related Classes (Brief Context)
- **`radian3`**: Uses `radian` for 3D Euler angles (ZXY order).
- **`fquat`**: Quaternions for gimbal-free rotations, likely uses `radian` for construction.
- **`fmat3`**, **`fmat4`**: Matrices for transforms, may accept `radian` angles.
- **`fvec3`**, **`fvec4`**: Vectors that may use trigonometric results for rotations.

---
