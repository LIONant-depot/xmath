# `fvec4` Class Documentation

## Overview

The `fvec4` class represents a 4D vector using single-precision floating-point numbers (`float`), optimized for SIMD operations via SSE 
(using `floatx4`, likely `__m128`). It is part of the `xmath` namespace and is aligned to 16 bytes for performance. This class is ideal for graphics 
(e.g., homogeneous coordinates in 3D transformations, RGBA colors), simulations, and high-dimensional math where SIMD acceleration matters.

Unlike `fvec3_t`, `fvec4` is not templated and always uses SIMD padding. The fourth component (`w`) is often used for perspective division in graphics or as a weight in other contexts.

### Key Features
- **SIMD Optimization**: Leverages SSE for faster dot products, normalizations, and component-wise ops.
- **Uninitialized Default Construction**: Prioritizes speed; initialize explicitly.
- **Mutable vs. Immutable**: In-place methods (e.g., `Normalize()`) chainable; copies suffixed `_Copy`.
- **Homogeneous Support**: Methods like `HomogeneousCopy()` for dividing by `w`.
- **Extensive Swizzling**: Full HLSL-style permutations for all dimensions (scalars to `fvec4`).
- **Validation and Safety**: `isFinite()`, safe normalization, range checks.
- **Interoperability**: Constructs from `fvec2`, `fvec3`, spans, and converts to strings/arrays.

Header-only; requires `xmath_vector.h`. Performance-focused—use `const` for immutability.

### Comparison to Similar Classes
- Similar to GLM's `glm::vec4` but with SIMD built-in, fast inverse sqrt variants, and more swizzles.
- Analogous to Eigen's `Eigen::Vector4f` but game/graphics-oriented with homogenize, reflection, and grid snap; SIMD rivals Eigen's but with easier swizzling.
- Engine equivalents: Unreal's `FVector4`, but standalone and more math functions (e.g., `SmoothStep()`).

### Performance Notes
- Aligned for SSE; use aligned allocators if storing in arrays.
- `InvSqrtFast()` uses approximate instructions for speed (e.g., _mm_rsqrt_ps), less precise than `InvSqrt()`.
- For colors (RGBA), treat as [r,g,b,a]; for positions, [x,y,z,1].

### Quick Start Example
```cpp
#include <xmath_vector.h>  // Assuming includes

int main() {
    xmath::fvec4 v(1.0f, 2.0f, 3.0f, 1.0f);  // Homogeneous point
    v.Normalize();  // Normalize (w remains)
    fvec4 homog = v.HomogeneousCopy();  // Divide by w if !=1
    float dot = v.Dot(fvec4::fromOne());  // Sum of components

    std::cout << v.ToString() << std::endl;  // "[x, y, z, w]"
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `fvec4()` | Default (uninitialized). | `fvec4 v;` |
| `fvec4(fvec4&&)` | Move (noexcept). | Implicit. |
| `fvec4(const fvec4&)` | Copy (noexcept). | `fvec4 v2(v1);` |
| `fvec4(float x, float y, float z, float w)` | From components. | `fvec4(1,2,3,4);` |
| `fvec4(float value)` | All components `value`. | `fvec4(0.5f);` |
| `fvec4(const fvec3& other, float w)` | Extend 3D with w. | `fvec4(vec3, 1);` |
| `fvec4(float x, const fvec3& other)` | Prepend x to 3D. | `fvec4(1, vec3);` |
| `explicit fvec4(const floatx4& reg)` | From SSE register. | Low-level. |
| `fvec4(const fvec2& xy, const fvec2& zw)` | From two 2D vectors. | `fvec4(xy, zw);` |
| `fvec4(std::span<float> Span)` | From span (exactly 4 elements). | `std::array<float,4> arr{1,2,3,4}; fvec4 v(arr);` |

## Assignment and Conversion Operators

| Operator/Method | Description | Example |
|-----------------|-------------|---------|
| `fvec4& operator=(const fvec4&)` | Assignment. | `v1 = v2;` |
| `operator std::array<double,4>() const` | To doubles array. | `std::array<double,4> arr = v;` |
| `operator std::string() const` | To string. | `std::string s = v;` |
| `std::string ToString() const` | Explicit string. | `v.ToString();` // e.g., "[1, 2, 3, 4]" |
| `friend std::ostream& operator<<(std::ostream& os, const fvec4& vec)` | Stream. | `std::cout << v;` |

## Static Properties

Convenience vectors:

| Method | Returns | Description |
|--------|---------|-------------|
| `fromZero()` | [0,0,0,0] | Zero. |
| `fromOne()` | [1,1,1,1] | One. |
| `fromUnitX()` | [1,0,0,0] | Unit X. |
| `fromUnitY()` | [0,1,0,0] | Unit Y. |
| `fromUnitZ()` | [0,0,1,0] | Unit Z. |
| `fromUnitW()` | [0,0,0,1] | Unit W. |
| `fromUp()` | [0,1,0,0] | Up. |
| `fromDown()` | [0,-1,0,0] | Down. |
| `fromLeft()` | [-1,0,0,0] | Left. |
| `fromRight()` | [1,0,0,0] | Right. |
| `fromForward()` | [0,0,1,0] | Forward. |
| `fromBack()` | [0,0,-1,0] | Back. |
| `fromRandomUnitVector()` | Random unit (length 1). | Not constexpr. |

## Static Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Dot(const fvec4& a, const fvec4& b)` | `float` | Dot product. |
| `Min(const fvec4& a, const fvec4& b)` | `fvec4` | Component min. |
| `Max(const fvec4& a, const fvec4& b)` | `fvec4` | Component max. |
| `Lerp(const fvec4& a, const fvec4& b, float t)` | `fvec4` | Linear interp. |
| `Distance(const fvec4& a, const fvec4& b)` | `float` | Euclidean dist. |

Instance overloads exist (e.g., `this->Dot(a)`).

## Instance Methods - Basic Operations

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `Length() const` | `float` | Magnitude. | No |
| `LengthSq() const` | `float` | Squared mag. | No |
| `LimitLengthCopy(float MaxLength) const` | `fvec4` | Cap length copy. | No |
| `NormalizeCopy() const` | `fvec4` | Normalized copy. | No |
| `Normalize()` | `fvec4&` | Normalize in-place. | Yes |
| `HomogeneousCopy() const` | `fvec4` | Divide x,y,z by w copy. | No |
| `Homogenize()` | `fvec4&` | Homogenize in-place (w=1 after). | Yes |
| `NormalizeSafeCopy() const` | `fvec4` | Safe normalize copy. | No |
| `NormalizeSafe()` | `fvec4&` | Safe in-place. | Yes |
| `isFinite() const` | `bool` | All finite. | No |
| `isInRange(float min, float max) const` | `bool` | Components in range. | No |
| `Reflection(const fvec4& normal) const` | `fvec4` | Reflect over normal. | No |
| `DistanceSquare(const fvec4& v) const` | `float` | Squared dist. | No |
| `AngleBetween(const fvec4& v) const` | `radian` | Angle between. | No |
| `GridSnap(float gridX, float gridY, float gridZ, float gridW)` | `fvec4&` | Snap to grid. | Yes |

## Instance Methods - Component-Wise Math

Copy/mutable pairs for per-component ops:

| Method (Copy/Mutable) | Description |
|-----------------------|-------------|
| `AbsCopy() / Abs()` | Absolute. |
| `OneOverCopy() / OneOver()` | Reciprocal. |
| `SqrtCopy() / Sqrt()` | Sqrt. |
| `InvSqrtFastCopy() / InvSqrtFast()` | Fast approx inv sqrt. |
| `InvSqrtCopy() / InvSqrt()` | Precise inv sqrt. |
| `SignCopy() / Sign()` | Sign. |
| `FloorCopy() / Floor()` | Floor. |
| `CeilCopy() / Ceil()` | Ceil. |
| `FractCopy() / Fract()` | Fractional. |
| `RoundCopy() / Round()` | Round. |
| `TruncCopy() / Trunc()` | Trunc. |
| `ModCopy(float) / Mod(float)` | Modulo. |
| `ClampCopy(float, float) / Clamp(float, float)` | Scalar clamp. |
| `ClampCopy(fvec4, fvec4) / Clamp(fvec4, fvec4)` | Vector clamp. |
| `Step(float) const` | Step. |
| `SmoothStep(float, float) const` | Smooth step. |
| `LogCopy() / Log()` | Ln. |
| `Log2Copy() / Log2()` | Log2. |
| `PowCopy(float) / Pow(float)` | Power. |
| `SinCopy() / Sin()` | Sin. |
| `CosCopy() / Cos()` | Cos. |
| `TanCopy() / Tan()` | Tan. |
| `AsinCopy() / Asin()` | Asin. |
| `AcosCopy() / Acos()` | Acos. |
| `AtanCopy() / Atan()` | Atan. |
| `Atan2Copy(fvec4) / Atan2(fvec4)` | Atan2 (this=y). |

## Swizzle Methods

HLSL-style component rearrangement (all const noexcept):

- Scalars: `x()`, `y()`, `z()`, `w()`
- `fvec2`: `xx()`, `xy()`, `xz()`, `xw()`, `yx()`, ..., `ww()` (16 total)
- `fvec3`: `xxx()`, `xxy()`, ..., `www()` (64 total)
- `fvec4`: `xxxx()`, `xxxy()`, ..., `wwww()` (256 total, all permutations)

**Example**:
```cpp
fvec4 v(1,2,3,4);
fvec2 rg = v.xy();  // (1,2) for RGBA
fvec3 rgb = v.xyz();  // (1,2,3)
fvec4 rgba = v.xyzw();  // (1,2,3,4)
fvec4 argb = v.wxyz();  // (4,1,2,3)
```

Swizzles enable efficient shader-like code without temporaries.

## Operator Overloads

| Operator | Description | Example |
|----------|-------------|---------|
| `+ (const fvec4&)` | Add. | `v1 + v2` |
| `- (const fvec4&)` | Subtract. | |
| `* (float)` | Scalar mul. | `v * 2` |
| `/ (float)` | Scalar div. | |
| `+= (const fvec4&)` | Add-assign. | |
| `-= (const fvec4&)` | Sub-assign. | |
| `*= (float)` | Mul-assign. | |
| `/= (float)` | Div-assign. | |
| `== (const fvec4&)` | Equal. | |
| `!= (const fvec4&)` | Not equal. | |
| `[] (std::int32_t) const` | Access. | `v[2]` (z) |
| `[] (std::int32_t)` | Mutable. | `v[3] = 1;` |
| `friend * (float, const fvec4&)` | Left scalar. | `2 * v` |
| `friend - (const fvec4&)` | Negate. | `-v` |

## Tutorials and Examples

### Basic Operations and Graphics
Homogeneous transform:
```cpp
fvec4 pos(1, 2, 3, 1);  // Point
// Assume matrix mul (not in class, but external)
pos.Homogenize();  // If w changed, divide
fvec4 normal(0,1,0,0);  // Direction (w=0)
fvec4 reflected = pos.Reflection(normal);
```

### Component-Wise for Colors
RGBA manipulation:
```cpp
fvec4 color(0.2f, 0.4f, 0.6f, 1.0f);  // Blue-ish
color.Clamp(0,1);  // Ensure range
color *= 1.5f;  // Brighten
color.Pow(2.2f);  // Gamma correct
fvec4 grayscale = color.xxxw();  // (avg, avg, avg, a) but compute avg first
```

### Advanced: Normalization and Angles
Direction vectors:
```cpp
fvec4 dir = fvec4::fromRandomUnitVector();
if (!dir.NormalizeSafe()) { /* Handle zero */ }
radian angle = dir.AngleBetween(fvec4::fromUnitX());
dir.GridSnap(0.1f, 0.1f, 0.1f, 1.0f);  // Snap xyz, keep w
```

### SIMD-Specific Tips
- Batch operations: Process arrays of `fvec4` for SSE benefits.
- Fast vs. Precise: Use `InvSqrtFast()` in inner loops for approx (error ~0.001); `InvSqrt()` for accuracy.

### Error Handling
- `NormalizeSafe()` returns zero vector if length zero.
- Check `isFinite()` post-computation to catch NaN/Inf.
- For homogeneous, ensure w != 0 before `Homogenize()`.

## Related Classes
- `fvec2`, `fvec3`: Lower dimensions.
- `floatx4`: Internal SIMD type.
- Matrices/Quaternions in `xmath` for transformations.

---