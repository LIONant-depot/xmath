﻿# `fquat_t` Class Documentation

## Overview

The `fquat_t` class is a templated 4D quaternion implementation in the `xmath` namespace, supporting single-precision floating-point (`float`) 
components. It comes in two specializations:
- **`fquat` (likely alias for `fquat`)**: SIMD-optimized using SSE (`floatx4`, assumed `__m128`). Aligned to 16 bytes, pads with `m_W` (real part) 
      and `m_XYZ` (imaginary vector). Faster for rotations but uses more memory and requires aligned allocation.
- **`fquat` CPU variant (likely alias for `fquatd`)**: Standard CPU version with exactly 4 floats, no SIMD or padding. Portable but slower for bulk operations.

Quaternions are ideal for 3D rotations, avoiding gimbal lock, and are used in graphics, physics, and animations. The class assumes right-hand rule for rotations and prioritizes performance with mutable (in-place, chainable) and immutable (copy-suffixed) operations.

### Key Features
- **SIMD vs. CPU**: Use SIMD for speed in loops; CPU for memory efficiency. Conversions between variants supported.
- **Uninitialized Default**: No initialization for speed; set explicitly.
- **Mutable/Immutable**: In-place methods (e.g., `Normalize()`) return `*this`; copies suffixed `_Copy`.
- **Interpolation**: `Lerp`, `Slerp` (spherical), `Squad` for smooth rotations.
- **Construction**: From axis-angle, Euler (`radian3`), from-to vectors, look-at.
- **Validation**: `isFinite`, `isNormalized`, `isNearlyIdentity`, `isNearlyZero`.
- **Operations**: Conjugate, inverse, dot product, angle between, vector rotation.
- **Operators**: Arithmetic, multiplication (composition), vector rotation (`q * v`).
- **No Swizzles**: Quaternions require all components; no HLSL-style swizzles.

Requires `xmath_flinear.h`. Emphasizes performance—use `const` for safety. Assumes right-hand rule.

### Comparison to Similar Classes
- Like GLM's `glm::quat` but with SIMD specialization, safe ops, and variants like `SlerpAccurate` for precision.
- Similar to Eigen's `Eigen::Quaternionf` with better game features (e.g., `LookRotation`, `RotateTowards`); SIMD variant rivals Eigen's but with easier vector rotation (`q * v`).
- Engine analogs: Unreal's `FQuat`, Godot's `Quaternion`, but standalone and templated.

### Performance Considerations
- **SIMD (`fquat`)**: ~4x faster for dot/Slerp on SSE hardware; ensure aligned memory.
- **CPU (`fquatd`)**: Scalar ops; use for non-SSE platforms or memory saving.
- **Interpolation**: `LerpFast` for quick approx; `SlerpAccurate` for precision.
- Tip: Convert CPU to SIMD for hot paths: `fquat simd_q(cpu_q);`.

### Quick Start Example
```cpp
#include "xmath_flinear.h"  // Assumes includes

int main() {
    xmath::fquat q = xmath::fquat::fromAxisAngle(fvec3::fromUnitX(), xmath::pi_over2_v);  // 90° around X
    q.Normalize();  // Ensure unit
    fvec3 rotated = q * fvec3::fromUnitY();  // Rotate (0,1,0) to (0,0,1)
    fquat target = xmath::fquat::fromIdentity();
    fquat interp = q.Slerp(target, 0.5f);  // Halfway

    std::cout << q.ToString() << std::endl;  // e.g., "[0.707, 0, 0, 0.707]"
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `fquat_t()` | Default (uninitialized). | `fquat q;` |
| `fquat_t(float x, float y, float z, float w)` | From components (imaginary x,y,z; real w). | `fquat(0,0,0,1);` // Identity |
| `fquat_t(const fvec3& axis, radian angle)` | From axis-angle. | `fquat(fvec3::fromUnitX(), xmath::pi_over2_v);` |
| `fquat_t(const radian3& euler)` | From Euler angles (ZXY order). | `fquat(radian3(0,0,xmath::pi_over2_v));` |
| `fquat_t(const fvec3& from, const fvec3& to, const fvec3& up)` | From-to rotation. | `fquat(start, end, up);` |
| `fquat_t(const fvec3& forward, const fvec3& up)` | Look rotation. | `fquat(dir, up);` |
| `explicit fquat_t(const floatx4& reg)` (SIMD only) | From SSE register. | Low-level. |
| `fquat_t(const fquat_t<!T_USE_SIMD_V>& other)` | Convert between SIMD/CPU. | `fquat simd(cpu_q);` |
| `fquat_t(const std::array<double, 4>& conversion)` | From doubles. | `fquat(arr);` |

## Assignment and Conversion Operators

| Operator/Method | Description | Example |
|-----------------|-------------|---------|
| `operator std::array<double, 4>() const` | To doubles array. | `std::array<double,4> arr = q;` |
| `operator std::string() const` | To string (e.g., "[x, y, z, w]"). | `std::string s = q;` |
| `std::string ToString() const` | Explicit string. | `q.ToString();` |
| `template <bool V> friend std::ostream& operator<<(std::ostream& os, const fquat_t<V>& quat)` | Stream output. | `std::cout << q;` |

## Static Properties

| Method | Returns | Description |
|--------|---------|-------------|
| `fromIdentity()` | Identity quaternion [0,0,0,1]. | No rotation. |
| `fromZero()` | Zero quaternion [0,0,0,0]. | Invalid for rotations. |
| `fromAxisAngle(const fvec3& axis, radian angle)` | From axis-angle. | Unit axis assumed. |
| `RandomUnitQuaternion()` | Random unit quaternion. | For random rotations. |

## Static Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Dot(const fquat_t& a, const fquat_t& b)` | `float` | Dot product. |
| `Lerp(const fquat_t& a, const fquat_t& b, float t)` | `fquat_t` | Linear interp (normalized). |
| `LerpFast(const fquat_t& Start, const fquat_t& End, float T)` | `fquat_t` | Fast lerp (approx). |
| `LerpUnclamped(const fquat_t& a, const fquat_t& b, float t)` | `fquat_t` | Lerp without t clamp. |
| `Slerp(const fquat_t& a, const fquat_t& b, float t)` | `fquat_t` | Spherical lerp. |
| `SlerpAccurate(const fquat_t& a, const fquat_t& b, float t)` | `fquat_t` | Precise slerp. |
| `SlerpUnclamped(const fquat_t& a, const fquat_t& b, float t)` | `fquat_t` | Slerp without t clamp. |
| `Squad(const fquat_t& a, const fquat_t& a_tangent, const fquat_t& b, const fquat_t& b_tangent, float t)` | `fquat_t` | Squad interpolation. |
| `AngleBetween(const fquat_t& a, const fquat_t& b)` | `radian` | Angle between quaternions. |
| `FromToRotation(const fvec3& from, const fvec3& to, const fvec3& up)` | `fquat_t` | Rotation from-to. |
| `LookRotation(const fvec3& forward, const fvec3& up)` | `fquat_t` | Look-at rotation. |

Instance versions available (e.g., `this->Dot(other)`).

## Instance Methods - Basic Operations

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `Length() const` | `float` | Magnitude. | No |
| `LengthSq() const` | `float` | Squared magnitude. | No |
| `NormalizeCopy() const` | `fquat_t` | Normalized copy. | No |
| `Normalize()` | `fquat_t&` | Normalize in-place. | Yes |
| `NormalizeSafeCopy() const` | `fquat_t` | Safe normalize copy. | No |
| `NormalizeSafe()` | `fquat_t&` | Safe normalize in-place. | Yes |
| `isFinite() const` | `bool` | All finite. | No |
| `isNormalized(float tolerance) const` | `bool` | Length ≈1 within tol. | No |
| `isNearlyIdentity(float tolerance) const` | `bool` | Near [0,0,0,1]. | No |
| `isNearlyZero(float tolerance) const` | `bool` | Near [0,0,0,0]. | No |
| `Equals(const fquat_t& other, float tolerance) const` | `bool` | Approx equal. | No |

## Instance Methods - Quaternion Specifics

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `ConjugateCopy() const` | `fquat_t` | Conjugate copy [-x,-y,-z,w]. | No |
| `Conjugate()` | `fquat_t&` | Conjugate in-place. | Yes |
| `InverseCopy() const` | `fquat_t` | Inverse copy (conjugate / length²). | No |
| `Inverse()` | `fquat_t&` | Inverse in-place. | Yes |
| `Axis() const` | `fvec3` | Rotation axis. | No |
| `Angle() const` | `radian` | Rotation angle. | No |
| `ToEuler() const` | `radian3` | To Euler angles (ZXY). | No |
| `Forward() const` | `fvec3` | Forward direction. | No |
| `Up() const` | `fvec3` | Up direction. | No |
| `Right() const` | `fvec3` | Right direction. | No |
| `ToAxisAngle() const` | `std::pair<fvec3, radian>` | Axis and angle. | No |
| `Delta(const fquat_t& other) const` | `fquat_t` | Delta rotation. | No |
| `LogCopy() const` | `fquat_t` | Logarithm copy. | No |
| `Log()` | `fquat_t&` | Logarithm in-place. | Yes |
| `ExpCopy() const` | `fquat_t` | Exponential copy. | No |
| `Exp()` | `fquat_t&` | Exponential in-place. | Yes |

## Instance Methods - Rotation Operations

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `setupRotationX(radian rx)` | `fquat_t&` | Set to X rotation. | Yes |
| `setupRotationY(radian ry)` | `fquat_t&` | Set to Y rotation. | Yes |
| `setupRotationZ(radian rz)` | `fquat_t&` | Set to Z rotation. | Yes |
| `setupLookRotation(const fvec3& forward, const fvec3& up)` | `fquat_t&` | Set to look-at. | Yes |
| `setupFromToRotation(const fvec3& from, const fvec3& to, const fvec3& up)` | `fquat_t&` | Set to from-to. | Yes |
| `RotateXCopy(radian rx) const` | `fquat_t` | Post-rotate X copy. | No |
| `RotateYCopy(radian ry) const` | `fquat_t` | Post-rotate Y copy. | No |
| `RotateZCopy(radian rz) const` | `fquat_t` | Post-rotate Z copy. | No |
| `RotateX(radian rx)` | `fquat_t&` | Post-rotate X. | Yes |
| `RotateY(radian ry)` | `fquat_t&` | Post-rotate Y. | Yes |
| `RotateZ(radian rz)` | `fquat_t&` | Post-rotate Z. | Yes |
| `PreRotateX(radian rx)` | `fquat_t&` | Pre-rotate X. | Yes |
| `PreRotateY(radian ry)` | `fquat_t&` | Pre-rotate Y. | Yes |
| `PreRotateZ(radian rz)` | `fquat_t&` | Pre-rotate Z. | Yes |
| `PreRotateXCopy(radian rx) const` | `fquat_t` | Pre-rotate X copy. | No |
| `PreRotateYCopy(radian ry) const` | `fquat_t` | Pre-rotate Y copy. | No |
| `PreRotateZCopy(radian rz) const` | `fquat_t` | Pre-rotate Z copy. | No |
| `RotateTowardsCopy(const fquat_t& target, radian maxDelta) const` | `fquat_t` | Rotate towards copy. | No |
| `RotateTowards(const fquat_t& target, radian maxDelta)` | `fquat_t&` | Rotate towards in-place. | Yes |
| `RotateVector(const fvec3& v) const` | `fvec3` | Rotate vector. | No |

## Operator Overloads

| Operator | Description | Example |
|----------|-------------|---------|
| `+ (const fquat_t&)` | Add. | `q1 + q2` |
| `- (const fquat_t&)` | Subtract. | |
| `* (float)` | Scalar mul. | `q * 2.0f` |
| `/ (float)` | Scalar div. | |
| `* (const fquat_t&)` | Composition. | `q1 * q2` |
| `+= (const fquat_t&)` | Add-assign. | |
| `-= (const fquat_t&)` | Sub-assign. | |
| `*=(float)` | Mul-assign. | |
| `/=(float)` | Div-assign. | |
| `*=(const fquat_t&)` | Compose-assign. | |
| `== (const fquat_t&)` | Exact equal. | |
| `!= (const fquat_t&)` | Not equal. | |
| `[] (std::int32_t) const` | Access (0=x,1=y,2=z,3=w). | `q[3]` (w) |
| `[] (std::int32_t)` | Mutable access. | `q[0] = 0.707f;` |
| `friend * (float, const fquat_t&)` | Left scalar. | `2.0f * q` |
| `friend - (const fquat_t&)` | Negate. | `-q` |
| `friend * (const fquat_t&, const fvec3&)` | Rotate vector. | `q * v` |

## Tutorials and Examples

### Basic Quaternion Creation and Rotation
Create and rotate:
```cpp
fquat q = fquat::fromAxisAngle(fvec3::fromUnitY(), xmath::pi_over2_v);  // 90° around Y
fvec3 rotated = q * fvec3::fromUnitX();  // (1,0,0) to (0,0,-1)
q.Inverse();  // Invert rotation
```

### Interpolation for Animation
Smooth rotation:
```cpp
fquat start = fquat::fromIdentity();
fquat end = fquat::fromAxisAngle(fvec3::fromUnitZ(), xmath::pi_v);
fquat slerp = start.Slerp(end, 0.5f);  // Halfway
fquat squad = start.Squad(tangentA, end, tangentB, 0.5f);  // Cubic
```

### From-To and Look Rotation
Camera look-at:
```cpp
fquat look = fquat::LookRotation(target - position, fvec3::fromUp());
fvec3 forward = look.Forward();  // Direction
radian3 euler = look.ToEuler();  // To Euler
```

### Advanced: Delta and Log/Exp
Compute delta rotation:
```cpp
fquat delta = start.Delta(end);
fquat log = delta.LogCopy();  // Log for advanced interp
fquat exp = log.ExpCopy();  // Exp back
```

### Performance Tips
- **SIMD**: Use `fquat` for batched rotations; align arrays to 16 bytes.
- **Normalization**: Call `NormalizeSafe()` once after construction/modification.
- **Interpolation**: `LerpFast` for speed; `SlerpAccurate` for precision.
- **Vector Rotation**: `q * v` is optimized; use for many vectors.

## Related Classes
- `radian`, `radian3`: For angles/Euler.
- `fvec3`: For axes/vectors.
- `fmat3`, `fmat4`: For matrix conversions (assumed external).
- `xmath_trigonometry.h`: For trig functions used in quats.

---
