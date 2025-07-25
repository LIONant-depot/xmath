# `xmath_functions.h` Documentation

## Overview

The `xmath_functions.h` header in the `xmath` namespace provides a collection of mathematical constants and functions for common numerical operations, 
focusing on floating-point arithmetic. It includes exponential, logarithmic, power, modulo, rounding, clamping, interpolation, and validation functions, 
all designed for performance in applications like game development, graphics, and scientific computing. The functions are templated where possible for 
flexibility with `float` or `double`, and many are `constexpr` or `inline noexcept` for compile-time evaluation and efficiency.

This header is header-only and uses `std::numbers` for constants. It complements other `xmath` utilities (e.g., trigonometry in `xmath_trigonometry.h`), 
providing general-purpose math without dependencies on vectors or matrices.

### Key Features
- **Constants**: Predefined mathematical values like e, ln(2), sqrt(2), and a float tolerance (0.001f).
- **Exponential and Logarithmic**: `Exp`, `Pow`, `Log`, `Log2`, `Log10`.
- **Modulo and Rounding**: `FMod`, `ModFX`, `Round`, `Ceil`, `Floor`, `Trunc`, `LRound`.
- **Clamping and Interpolation**: `Clamp`, `Lerp`, `Range`, `isInRange`.
- **Validation**: `FEqual`, `FLess`, `FGreater`, `isValid`, `isFinite`.
- **Miscellaneous**: `Sqr`, `Sqrt`, `InvSqrt`, `Abs`, `Sign`, `CopySign`, `FSel`, `i2f`, `f2i`, `SolvedQuadraticRoots`.
- **Type Safety**: Templated for floating-point and integer types; uses concepts like `std::floating_point`.
- **Performance**: Many functions are `constexpr` for compile-time use; all `noexcept` for safety.

### Comparison to Similar Libraries
- Unlike GLM's math functions (which are vector-focused), `xmath_functions.h` is scalar-oriented with extras like `ModFX` (modulo with quotient output) and `SolvedQuadraticRoots`.
- Compared to Eigen's utilities, it's lighter and C++20-focused (e.g., concepts for templates), with game-friendly features like `flt_tol_v` 
     for approximate comparisons. Eigen integrates better with matrices, but this header is standalone.
- Engine analogs: Similar to Unity's `Mathf` or Unreal's `FMath`, but more templated and constexpr-capable for modern C++.

### Usage Notes
- **Tolerance**: Use `flt_tol_v` (0.001f) for floating-point comparisons (`FEqual`, `FLess`, `FGreater`).
- **Precision**: Templated for `float`/`double`; use `double` for high accuracy (e.g., `Exp<double>(x)`).
- **Signed/Unsigned**: Functions like `Abs` and `Sign` work with signed types; check constraints.
- **Quadratic Solver**: `SolvedQuadraticRoots` returns true if real roots found, setting `root1`/`root2`.
- **Performance**: Inline and constexpr; suitable for hot paths.

### Quick Start Example
```cpp
#include "xmath_functions.h"

int main() {
    float x = 2.0f;
    float exp = xmath::Exp(x);  // e^2 ≈ 7.389
    float pow = xmath::Pow(x, 3.0f);  // 8.0
    float sqrt = xmath::Sqrt(x);  // ≈1.414
    float clamped = xmath::Clamp(3.5f, 0.0f, 2.0f);  // 2.0
    bool equal = xmath::FEqual(1.0005f, 1.0f);  // true (within tol)

    float root1, root2;
    if (xmath::SolvedQuadraticRoots(root1, root2, 1.0f, -3.0f, 2.0f)) {  // x^2 - 3x + 2 = 0
        // root1=1.0, root2=2.0
    }

    return 0;
}
```

## Constants

| Constant | Type | Value | Description |
|----------|------|-------|-------------|
| `e_v` | `float` | e ≈ 2.71828 | Base of natural logarithm. |
| `log2e_v` | `float` | log₂(e) ≈ 1.4427 | Log base 2 of e. |
| `log10e_v` | `float` | log₁₀(e) ≈ 0.434294 | Log base 10 of e. |
| `ln2_v` | `float` | ln(2) ≈ 0.693147 | Natural log of 2. |
| `ln10_v` | `float` | ln(10) ≈ 2.30259 | Natural log of 10. |
| `sqrt2_v` | `float` | √2 ≈ 1.41421 | Square root of 2. |
| `flt_tol_v` | `float` | 0.001f | Default tolerance for float comparisons. |

**Notes**: Sourced from `std::numbers`; all `constexpr`.

## Functions

### Exponential and Power Functions

| Function | Returns | Description | Example |
|----------|---------|-------------|---------|
| `Exp<T>(T x)` | `T` | e^x. | `xmath::Exp(1.0f)` // e ≈ 2.718 |
| `Pow<T>(T a, T b)` | `T` | a^b. | `xmath::Pow(2.0f, 3.0f)` // 8.0 |

**Notes**: Templated for `std::floating_point T`; `inline noexcept`.

### Modulo Functions

| Function | Returns | Description | Example |
|----------|---------|-------------|---------|
| `FMod<T>(T x, T y)` | `T` | x mod y (floating-point). | `xmath::FMod(5.3f, 2.0f)` // 1.3 |
| `ModFX<T>(T x, T& y)` | `T` | x mod y, sets y to quotient. | `float q; xmath::ModFX(5.3f, q)` // 1.3, q=2.0 |

**Notes**: Templated for `std::floating_point T`; `inline noexcept`.

### Logarithmic Functions

| Function | Returns | Description | Example |
|----------|---------|-------------|---------|
| `Log<T>(T x)` | `T` | Natural log (ln x). | `xmath::Log(xmath::e_v)` // 1.0 |
| `Log2<T>(T x)` | `T` | Log base 2. | `xmath::Log2(8.0f)` // 3.0 |
| `Log10<T>(T x)` | `T` | Log base 10. | `xmath::Log10(100.0f)` // 2.0 |

**Notes**: Templated for `std::floating_point T`; `inline noexcept`.

### Conversion Functions

| Function | Returns | Description | Example |
|----------|---------|-------------|---------|
| `i2f<T>(T i)` | `T` | Int to float (constexpr). | `xmath::i2f(5)` // 5.0f |
| `f2i<T>(T f)` | `std::int32_t` | Float to int (constexpr). | `xmath::f2i(5.9f)` // 5 |

**Notes**: Templated for `std::floating_point T`; `constexpr noexcept`.

### Miscellaneous Functions

| Function | Returns | Description | Example |
|----------|---------|-------------|---------|
| `FSel<T>(T a, T b, T c)` | `T` | Select b if a >= 0, else c. | `xmath::FSel(1.0f, 2.0f, 3.0f)` // 2.0 |
| `Sqr<T>(T x)` | `T` | x². | `xmath::Sqr(3.0f)` // 9.0 |
| `Sqrt<T>(T x)` | `T` | √x. | `xmath::Sqrt(4.0f)` // 2.0 |
| `InvSqrt<T>(T x)` | `T` | 1/√x. | `xmath::InvSqrt(4.0f)` // 0.5 |
| `Min<T1, T2>(T1 a, T2 b)` | `decltype(a + b)` | Minimum of a and b. | `xmath::Min(1.0f, 2)` // 1.0f |
| `Max<T1, T2>(T1 a, T2 b)` | `decltype(a + b)` | Maximum of a and b. | `xmath::Max(1.0f, 2)` // 2 |
| `FEqual<T>(T f0, T f1, T tol = flt_tol_v)` | `bool` | f0 ≈ f1 within tol. | `xmath::FEqual(1.0005f, 1.0f)` // true |
| `FLess<T>(T f0, T f1, T tol = flt_tol_v)` | `bool` | f0 < f1 with tol. | `xmath::FLess(0.999f, 1.0f)` // true |
| `FGreater<T>(T f0, T f1, T tol = flt_tol_v)` | `bool` | f0 > f1 with tol. | `xmath::FGreater(1.001f, 1.0f)` // true |
| `Sign<T>(T x)` | `bool` | true if x >= 0. | `xmath::Sign(-1.0f)` // false |
| `LRound<T>(T x)` | `std::int32_t` | Round to nearest int. | `xmath::LRound(2.6f)` // 3 |
| `Round<T>(T a, T b)` | `T` | Round a to nearest multiple of b. | `xmath::Round(5.3f, 2.0f)` // 6.0 |
| `Round<T>(T x)` | `T` | Round to nearest integer. | `xmath::Round(2.6f)` // 3.0 |
| `Ceil<T>(T x)` | `T` | Ceiling. | `xmath::Ceil(2.3f)` // 3.0 |
| `Floor<T>(T x)` | `T` | Floor. | `xmath::Floor(2.7f)` // 2.0 |
| `isInRange<T>(T x, T min, T max)` | `bool` | x in [min, max]. | `xmath::isInRange(5, 1, 10)` // true |
| `Range<T>(T x, T min, T max)` | `T` | Clamp x to [min, max]. | `xmath::Range(15, 1, 10)` // 10 |
| `Abs<T>(T x)` | `T` | Absolute value. | `xmath::Abs(-5.0f)` // 5.0 |
| `Clamp<T>(const T& value, const T& low, const T& high)` | `T` | Clamp value to [low, high]. | `xmath::Clamp(15.0f, 1.0f, 10.0f)` // 10.0 |
| `Lerp<T>(float t, T a, T b)` | `T` | Linear interpolation. | `xmath::Lerp(0.5f, 0.0f, 10.0f)` // 5.0 |
| `isValid<T>(T x)` | `bool` | x is finite. | `xmath::isValid(INFINITY)` // false |
| `Trunc<T>(T x)` | `T` | Truncate toward zero. | `xmath::Trunc(2.7f)` // 2.0 |
| `isFinite<T>(T x)` | `bool` | x is finite (no NaN/Inf). | `xmath::isFinite(NAN)` // false |
| `CopySign<T>(T x, T y)` | `T` | x with sign of y. | `xmath::CopySign(5.0f, -1.0f)` // -5.0 |

**Notes**: Templated where applicable; `constexpr` for compile-time use.

## Quadratic Solver

| Function | Returns | Description | Example |
|----------|---------|-------------|---------|
| `SolvedQuadraticRoots(float& root1, float& root2, float a, float b, float c)` | `bool` | Solve ax² + bx + c = 0; true if real roots. | `float r1, r2; xmath::SolvedQuadraticRoots(r1, r2, 1, -3, 2)` // true, r1=1, r2=2 |

**Notes**: `inline noexcept`; sets roots if discriminant >=0.

## Tutorials and Examples

### Exponential and Power in Physics
Calculate growth:
```cpp
float growth = xmath::Exp(0.05f * 10.0f);  // e^{0.5} ≈ 1.6487
float power = xmath::Pow(1.05f, 10.0f);  // 1.05^10 ≈ 1.6289
```

### Modulo and Rounding for Games
Wrap positions:
```cpp
float pos = xmath::FMod(5.3f, 2.0f);  // 1.3
float q;
float mod = xmath::ModFX(5.3f, q);  // 1.3, q=2.0
float rounded = xmath::Round(5.3f, 2.0f);  // 6.0
```

### Clamping and Interpolation for UI
Lerp animation:
```cpp
float val = xmath::Lerp(0.5f, 0.0f, 10.0f);  // 5.0
float clamped = xmath::Clamp(15.0f, 0.0f, 10.0f);  // 10.0
bool inRange = xmath::isInRange(5.0f, 0.0f, 10.0f);  // true
```

### Floating-Point Comparisons
Approximate equality:
```cpp
bool eq = xmath::FEqual(1.0005f, 1.0f);  // true
bool less = xmath::FLess(0.999f, 1.0f);  // true
bool finite = xmath::isFinite(1.0f / 0.0f);  // false
```

### Quadratic Equation Solver
Solve for intersections:
```cpp
float r1, r2;
if (xmath::SolvedQuadraticRoots(r1, r2, 1.0f, 5.0f, 6.0f)) {  // x² + 5x + 6 = 0
    // r1=-2, r2=-3
} else {
    // Complex roots
}
```

### Error Handling
- **Division by Zero**: Functions like `Pow`, `InvSqrt` may return Inf/NaN; check with `isFinite`.
- **Invalid Inputs**: Use `isValid` or `isFinite` post-calculation:
  ```cpp
  float res = xmath::Log(-1.0f);  // NaN
  if (!xmath::isFinite(res)) { /* Handle */ }
  ```
- **Tolerance**: Customize tol in comparisons:
  ```cpp
  bool eq = xmath::FEqual(a, b, 0.0001f);  // Tighter tol
  ```

## Performance Tips
- **Constexpr**: Use in templates or constexpr contexts:
  ```cpp
  constexpr float sq = xmath::Sqr(2.0f);  // 4.0 at compile-time
  ```
- **Inline**: All functions inline for zero overhead.
- **Tolerance**: Use `flt_tol_v` for default comparisons to avoid floating-point pitfalls.
- **Quadratic**: Efficient for real-time (e.g., raytracing); check return bool for roots.

## Related Classes (Brief Context)
- **`xmath_trigonometry.h`**: Uses similar templates for angles (`radian`, `degree`).
- **`radian3`**: Applies these functions for 3D rotations.
- **`fquat`**, **`fmat3`**, **`fmat4`**: May use these for transforms.

---