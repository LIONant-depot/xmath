# `fvec3_t` Class Documentation

## Overview

The `fvec3_t` class is a templated 3D vector implementation in the `xmath` namespace, supporting single-precision floating-point 
(`float`) components. It comes in two specializations:
- **`fvec3` (alias for `fvec3_t<true>`)**: SIMD-optimized using SSE instructions 
     (via `floatx4`, assumed to be `__m128`). Aligned to 16 bytes, pads with a dummy `m_W` component 
     (ignore its value). This variant is faster for vector operations but consumes extra memory 
     (4 floats instead of 3) and requires aligned allocation.
- **`fvec3d` (alias for `fvec3_t<false>`)**: Standard CPU version with exactly 3 floats, no SIMD or padding. 
     Portable but slower for bulk operations.

Designed for performance-critical applications like 3D graphics, physics simulations, and game engines. 
It extends `fvec2`'s features to 3D, adding cross products, Euler angles, and plane projections.

### Key Features
- **SIMD vs. CPU Variants**: Use `fvec3` for speed in loops; `fvec3d` for memory efficiency or non-SSE platforms. 
     Conversions between them are supported.
- **Uninitialized Default**: No initialization for speed; always set values explicitly.
- **Mutable/Immutable Ops**: In-place methods (e.g., `Normalize()`) return `*this` for chaining; copies suffixed `_Copy`.
- **Constexpr/Inline**: Many operations are `constexpr` or `inline` for optimization.
- **Validation**: `isFinite()`, `isNormalized()`, assertions in debug.
- **Swizzling**: Extensive HLSL-style swizzles for component permutation.
- **Geometry Focus**: 3D-specific like rotation axes, pitch/yaw, reflections, and line/plane interactions.
- **Interoperability**: Converts to/from `std::array<double, 3>`, `std::span<float>`, `fvec2`, strings.

Requires `xmath_fvector.h`. Emphasizes performance—use `const` for safety.

### Comparison to Similar Classes
- Like GLM's `glm::vec3` but with SIMD specialization, safe ops, and more geometry (e.g., `PitchYaw()`).
- Similar to Eigen's `Eigen::Vector3f` with better game-oriented features (e.g., `MoveTowards()`) and swizzles; SIMD variant rivals Eigen's dense matrix speed.
- Engine analogs: Unreal's `FVector`, Godot's `Vector3`, but standalone and templated.

### Performance Considerations
- **SIMD (`fvec3`)**: Leverages SSE for dot/cross/normalize (~4x faster on compatible hardware). 
     Ensure aligned memory (e.g., via `alignas(16)` or allocators).
- **CPU (`fvec3d`)**: Scalar ops; use for embedded systems or when alignment is an issue.
- Benchmarks: In loops, SIMD excels for batched vectors; test with your workload.
- Tip: Convert `fvec3d` to `fvec3` for hot paths: `fvec3 simd_vec(cpu_vec);`.

### Quick Start Example
```cpp
#include <xmath_fvector.h>  // Assumes necessary includes

int main() {
    xmath::fvec3 v1(1.0f, 2.0f, 3.0f);  // SIMD variant
    xmath::fvec3d v2(4.0f, 5.0f, 6.0f); // CPU variant
    xmath::fvec3 cross = v1.Cross(v2);  // Cross product
    float dot = v1.Dot(v2);  // 32.0f
    v1.Normalize();  // In-place normalize

    std::cout << v1.ToString() << std::endl;  // e.g., "[0.267, 0.535, 0.802]"
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `fvec3_t()` | Default (uninitialized). | `fvec3 v;` |
| `fvec3_t(float x, float y, float z)` | From components. | `fvec3(1, 2, 3);` |
| `fvec3_t(float value)` | All components `value`. | `fvec3(5.0f);` // [5, 5, 5] |
| `fvec3_t(radian pitch, radian yaw)` | From Euler angles (direction vector). | `fvec3(pitch, yaw);` |
| `explicit fvec3_t(const floatx4& reg)` (SIMD only) | From SSE register. | For low-level use. |
| `fvec3_t(const fvec3_t<!T_USE_SIMD_V>& other)` | Convert between SIMD/CPU. | `fvec3 simd(v2);` (from CPU) |
| `fvec3_t(const fvec2& other, float z)` | Extend 2D vector. | `fvec3(vec2, 0);` |
| `fvec3_t(float x, const fvec2& other)` | Prepend x to 2D. | `fvec3(0, vec2);` |
| `fvec3_t(std::span<float> Span)` | From span (exactly 3 elements). | `std::array<float, 3> arr{1,2,3}; fvec3 v(arr);` |
| `fvec3_t(const std::array<double,3>& Conversion)` | From doubles. | |

## Assignment and Conversion Operators

| Operator/Method | Description | Example |
|-----------------|-------------|---------|
| `operator std::array<double,3>() const` | To array of doubles. | `std::array<double,3> arr = v;` |
| `operator std::string() const` | To string (e.g., "[x, y, z]"). | `std::string s = v;` |
| `std::string ToString() const` | Explicit string. | `v.ToString();` |
| `template<bool V> friend std::ostream& operator<<(std::ostream& os, const fvec3_t<V>& vec)` | Stream output. | `std::cout << v;` |

## Static Properties

Factory methods for common vectors:

| Method | Returns | Description |
|--------|---------|-------------|
| `fromZero()` | [0, 0, 0] | Zero vector. |
| `fromOne()` | [1, 1, 1] | One vector. |
| `fromUp()` | [0, 1, 0] | Up (Y-axis). |
| `fromDown()` | [0, -1, 0] | Down. |
| `fromLeft()` | [-1, 0, 0] | Left (negative X). |
| `fromRight()` | [1, 0, 0] | Right (X-axis). |
| `fromForward()` | [0, 0, 1] | Forward (Z-axis). |
| `fromBack()` | [0, 0, -1] | Back. |
| `fromRandomUnitVector()` | Random unit vector. | Not `consteval`; uses RNG. |

## Static Methods

Operate without instance:

| Method | Returns | Description |
|--------|---------|-------------|
| `Dot(const fvec3_t& a, const fvec3_t& b)` | `float` | Dot product. |
| `Cross(const fvec3_t& a, const fvec3_t& b)` | `fvec3_t` | Cross product. |
| `Min(const fvec3_t& a, const fvec3_t& b)` | `fvec3_t` | Component min. |
| `Max(const fvec3_t& a, const fvec3_t& b)` | `fvec3_t` | Component max. |
| `Lerp(const fvec3_t& a, const fvec3_t& b, float t)` | `fvec3_t` | Linear interp. |
| `Distance(const fvec3_t& a, const fvec3_t& b)` | `float` | Euclidean distance. |

Instance versions available (e.g., `this->Dot(a)`).

## Instance Methods - Basic Operations

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `Length() const` | `float` | Magnitude. | No |
| `LengthSq() const` | `float` | Squared magnitude. | No |
| `NormalizeCopy() const` | `fvec3_t` | Normalized copy. | No |
| `Normalize()` | `fvec3_t&` | Normalize in-place. | Yes |
| `NormalizeSafeCopy() const` | `fvec3_t` | Safe (zero-handling). | No |
| `NormalizeSafe()` | `fvec3_t&` | Safe in-place. | Yes |
| `LimitLengthCopy(float MaxLength) const` | `fvec3_t` | Cap length. | No |
| `isFinite() const` | `bool` | All finite? | No |
| `isInRange(float min, float max) const` | `bool` | Components in range. | No |
| `Equals(const fvec3_t& other, float tolerance) const` | `bool` | Approx equal. | No |
| `DistanceSquare(const fvec3_t& v) const` | `float` | Squared distance. | No |

## Instance Methods - Component-Wise Math

Similar to `fvec2`, applied per component. Copy/mutable pairs:

| Method (Copy/Mutable) | Description |
|-----------------------|-------------|
| `AbsCopy() / Abs()` | Absolute. |
| `OneOverCopy() / OneOver()` | Reciprocal. |
| `SqrtCopy() / Sqrt()` | Square root. |
| `InvSqrtCopy() / InvSqrt()` | Inverse sqrt. |
| `SignCopy() / Sign()` | Sign. |
| `FloorCopy() / Floor()` | Floor. |
| `CeilCopy() / Ceil()` | Ceil. |
| `FractCopy() / Fract()` | Fractional. |
| `RoundCopy() / Round()` | Round. |
| `TruncCopy() / Trunc()` | Trunc. |
| `ModCopy(float) / Mod(float)` | Modulo. |
| `ClampCopy(float, float) / Clamp(float, float)` | Scalar clamp. |
| `ClampCopy(fvec3_t, fvec3_t) / Clamp(fvec3_t, fvec3_t)` | Vector clamp. |
| `Step(float) const` | Step function. |
| `SmoothStep(float, float) const` | Smooth step. |
| `LogCopy() / Log()` | Natural log. |
| `Log2Copy() / Log2()` | Log base 2. |
| `PowCopy(float) / Pow(float)` | Power. |
| `ExpCopy() / Exp()` | Exponential. |
| `SinCopy() / Sin()` | Sine. |
| `CosCopy() / Cos()` | Cosine. |
| `TanCopy() / Tan()` | Tangent. |
| `AsinCopy() / Asin()` | Arcsin. |
| `AcosCopy() / Acos()` | Arccos. |
| `AtanCopy() / Atan()` | Arctan. |
| `Atan2Copy(fvec3_t) / Atan2(fvec3_t)` | Atan2 (this as y). |

## Instance Methods - Geometry

3D-specific:

| Method | Returns | Description |
|--------|---------|-------------|
| `Pitch() const` | `radian` | Pitch angle. |
| `Yaw() const` | `radian` | Yaw angle. |
| `PitchYaw() const` | `std::pair<radian, radian>` | Both angles. |
| `RotationTowards(const fvec3_t& dest) const` | `std::pair<fvec3_t, radian>` | Axis-angle to dest. |
| `AngleBetween(const fvec3_t& v) const` | `radian` | Unsigned angle. |
| `SignedAngleBetween(const fvec3_t& v) const` | `radian` | Signed angle. |
| `VectorToLineSegment(const fvec3_t& start, const fvec3_t& end) const` | `fvec3_t` | Vector to segment. |
| `SquareDistToLineSeg(const fvec3_t& start, const fvec3_t& end) const` | `float` | Squared dist to segment. |
| `ClosestPointInLineSegment(const fvec3_t& start, const fvec3_t& end) const` | `fvec3_t` | Closest point. |
| `ClosestPointToRectangle(const fvec3_t& p0, const fvec3_t& e0, const fvec3_t& e1, fvec3_t& outClosestPoint) const` | `float` | Dist to rect, out point. |
| `RotateXCopy(radian rx) const` | `fvec3_t` | Rotate X copy. |
| `RotateYCopy(radian ry) const` | `fvec3_t` | Rotate Y copy. |
| `RotateZCopy(radian rz) const` | `fvec3_t` | Rotate Z copy. |
| `RotateX(radian rx)` | `fvec3_t&` | Rotate X in-place. |
| `RotateY(radian ry)` | `fvec3_t&` | Rotate Y in-place. |
| `RotateZ(radian rz)` | `fvec3_t&` | Rotate Z in-place. |
| `RotateCopy(const radian3& r) const` | `fvec3_t` | Euler rotate copy. |
| `Rotate(const radian3& r)` | `fvec3_t&` | Euler rotate. |
| `RotateInverse(const radian3& r)` | `fvec3_t&` | Inverse Euler. |
| `RotateInverseCopy(const radian3& r) const` | `fvec3_t` | Inverse copy. |
| `Reflection(const fvec3_t& normal) const` | `fvec3_t` | Reflect over normal. |
| `GridSnap(float gridX, float gridY, float gridZ)` | `fvec3_t&` | Snap to grid. |
| `isRightHanded(const fvec3_t& p1, const fvec3_t& p2) const` | `bool` | Right-handed check. |
| `ProjectCopy(const fvec3_t& onto) const` | `fvec3_t` | Project copy. |
| `Project(const fvec3_t& onto)` | `fvec3_t&` | Project in-place. |
| `Perpendicular(const fvec3_t& normal) const` | `fvec3_t` | Perp vector. |
| `ProjectOntoPlane(const fvec3_t& normal) const` | `fvec3_t` | Project to plane. |
| `isNearlyZero(float tolerance) const` | `bool` | Near zero. |
| `isNormalized(float tolerance) const` | `bool` | Unit length check. |
| `MoveTowardsCopy(const fvec3_t& target, float maxDistanceDelta) const` | `fvec3_t` | Move towards copy. |
| `MoveTowards(const fvec3_t& target, float maxDistanceDelta)` | `fvec3_t&` | Move towards. |

## Swizzle Methods

HLSL-style for rearrangement:

- Scalars: `x()`, `y()`, `z()`
- `fvec2`: `xx()`, `xy()`, `xz()`, `yx()`, `yy()`, `yz()`, `zx()`, `zy()`, `zz()`
- `fvec3_t`: `xxx()`, `xxy()`, ..., `zzz()` (all combos)
- `fvec4`: `xxxx()`, `xxxy()`, ..., `zzzz()` (extensive list)

Example:
```cpp
fvec3 v(1, 2, 3);
fvec2 xy = v.xy();  // (1, 2)
fvec3 zyx = v.zyx();  // (3, 2, 1)
fvec4 xyzw = v.xyzx();  // (1, 2, 3, 1)
```

## Operator Overloads

| Operator | Description | Example |
|----------|-------------|---------|
| `+ (const fvec3_t&)` | Add. | `v1 + v2` |
| `- (const fvec3_t&)` | Subtract. | |
| `* (const fvec3_t&)` | Component multiply. | |
| `* (float)` | Scalar multiply. | |
| `/ (float)` | Scalar divide. | |
| `+= (const fvec3_t&)` | Add-assign. | |
| `-= (const fvec3_t&)` | Subtract-assign. | |
| `*= (float)` | Multiply-assign. | |
| `/= (float)` | Divide-assign. | |
| `== (const fvec3_t&)` | Exact equal. | |
| `!= (const fvec3_t&)` | Not equal. | |
| `[] (std::int32_t) const` | Access (0=x,1=y,2=z). | `v[1]` |
| `[] (std::int32_t)` | Mutable access. | |
| `friend * (float, const fvec3_t&)` | Left scalar. | `2 * v` |
| `friend - (const fvec3_t&)` | Negate. | `-v` |

## Tutorials and Examples

### Vector Basics and Physics
Simulate velocity:
```cpp
fvec3 position(0, 0, 0);
fvec3 velocity(1, 0, 1);
velocity.Normalize();  // Unit direction
position += velocity * 0.1f;  // Update
if (velocity.isNormalized()) { /* Check */ }
```

### Geometry and Rotations
Camera direction:
```cpp
fvec3 forward = fvec3::fromForward();
forward.RotateY(xmath::pi_over4_v);  // 45° yaw
auto [pitch, yaw] = forward.PitchYaw();
fvec3 reflected = forward.Reflection(fvec3::fromUp());  // Bounce off ground
```

### Advanced: Line and Plane Interactions
Collision check:
```cpp
fvec3 point(1, 1, 1);
fvec3 start(0, 0, 0), end(2, 0, 0);
fvec3 closest = point.ClosestPointInLineSegment(start, end);  // (1, 0, 0)
float sqDist = point.SquareDistToLineSeg(start, end);
fvec3 projected = point.ProjectOntoPlane(fvec3::fromUp());  // (1, 0, 1)
```

### SIMD vs CPU Usage
```cpp
fvec3d cpu_vec(1,2,3);  // Memory-efficient
fvec3 simd_vec(cpu_vec);  // Convert to SIMD
simd_vec *= 2.0f;  // Faster op
fvec3d back = simd_vec;  // Convert back
```

### Error Handling
- Safe methods prevent NaN (e.g., `NormalizeSafe()` returns zero for zero vectors).
- Use `isFinite()` before ops; `isNearlyZero()` for epsilon checks.

## Related Classes
- `fvec2`, `fvec4`: Dimensional variants.
- `radian3`: For Euler rotations.
- `floatx4`: SIMD internal (SSE).

---