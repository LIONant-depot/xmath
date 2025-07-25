# `fmat3_t` Class Documentation

## Overview

The `fmat3_t` class is a templated 3x3 matrix implementation in the `xmath` namespace, using single-precision floating-point (`float`) components in column-major order. 
It supports two specializations:
- **SIMD Variant (`fmat3`)**: Optimized with SSE (`floatx4` columns, padded to 3x4), aligned to 16 bytes for faster operations like multiplication.
- **CPU Variant (`fmat3d`)**: Compact 3x3 layout without padding, for memory efficiency or non-SSE platforms.

Designed for linear transformations like rotation and scaling (no translation; use `fmat4_t` for affine). Assumes right-hand rule for rotations. 
Prioritizes performance with mutable (in-place, chainable) and immutable (Copy-suffixed) operations. Uninitialized by default; use `fromIdentity()` or setup methods.

### Key Features
- **Column-Major Storage**: Internal layout, but accessors (`()(row, col)`) feel row-major.
- **SIMD/CPU Toggle**: Template bool for optimization vs. compactness; conversions between variants.
- **Construction**: From identity, zero, quaternion, Euler, scale, rotation.
- **Operations**: Matrix/matrix mul, vector transform, transpose, inverse, determinant.
- **Geometry Helpers**: Extract rotation/scale, direction vectors (forward/up/right), rotate vectors.
- **Mutable Setup**: Chainable methods like `setupRotation`, `Rotate`, `Scale`.
- **Validation**: `isFinite`, `isIdentity`, `SanityCheck`.
- **No Initialization**: Default ctor leaves uninitialized; explicit setup required.
- **Performance**: Inline/constexpr where possible; assumes finite matrices (asserts in debug).

Requires `xmath_flinear.h`. Column-major for compatibility with OpenGL/DirectX conventions.

### Comparison to Similar Classes
- Like GLM's `glm::mat3` but with SIMD specialization, mutable chaining (e.g., `RotateX(angle)`), and safety checks like `SanityCheck`. GLM is more general-purpose; 
      this is optimized for 3D transforms.
- Similar to Eigen's `Eigen::Matrix3f` with better game features (e.g., `ExtractRotation` to quaternion, direction vectors); SIMD variant rivals Eigen's but with easier setup from Euler/quat.
- Engine analogs: Unreal's `FMatrix` (but 4x4), Godot's `Basis` (3x3), but standalone and templated.

### Performance Considerations
- **SIMD (`fmat3`)**: Faster mul/transform (~4x on SSE); align instances to 16 bytes.
- **CPU (`fmat3d`)**: Smaller size (36 vs 48 bytes); use for embedded or non-SSE.
- **Chaining**: Mutable methods return `*this` for fluency: `mat.setupIdentity().RotateX(angle).Scale(s);`.
- **Inverse/Det**: Assume invertible; check det !=0.
- Tip: Convert CPU to SIMD for hot paths: `fmat3_t<true> simd_mat(cpu_mat);`.

### Quick Start Example
```cpp
#include "xmath_flinear.h"  // Assumes includes

int main() {
    xmath::fmat3 mat = xmath::fmat3::fromIdentity();
    mat.RotateX(xmath::pi_over2_v);  // 90° around X
    fvec3 transformed = mat * fvec3::fromUnitY();  // (0,1,0) to (0,0,1)
    fquat rot = mat.ExtractRotation();  // To quat

    std::cout << mat(1,1) << std::endl;  // Access m11
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `fmat3_t()` | Default (uninitialized). | `fmat3 mat;` |
| `fmat3_t(float diagonal)` | Diagonal matrix. | `fmat3(2.0f);` // Scale 2 |
| `fmat3_t(const std::array<float, 9>& arr)` | From array (col-major). | `std::array<float,9> a = {1,0,0,0,1,0,0,0,1}; fmat3(a);` |
| `fmat3_t(std::span<const float, 9> span)` | From span (col-major). | Similar to array. |
| `fmat3_t(const fquat& q)` | From quaternion rotation. | `fmat3(quat);` |
| `fmat3_t(const radian3& Euler)` | From Euler angles (ZXY). | `fmat3(euler);` |
| `fmat3_t(const fquat& rotation, const fvec3& scale)` | From rotation + scale. | `fmat3(q, s);` |
| `explicit fmat3_t(const fmat3_t<!T_USE_SIMD_V>& other)` | Convert SIMD/CPU. | `fmat3(cpu_mat);` |

## Static Constructors

| Method | Returns | Description |
|--------|---------|-------------|
| `fromIdentity()` | Identity matrix. | No transform. |
| `fromZero()` | Zero matrix. | All zeros. |
| `fromRotation(const fquat& q)` | From quaternion. | Rotation matrix. |
| `fromRotation(const radian3& Euler)` | From Euler. | ZXY order. |
| `fromRotation(const fvec3& axis, radian angle)` | From axis-angle. | Rotation matrix. |
| `fromScale(const fvec3& s)` | Scale matrix. | Diagonal scale. |
| `fromRotationX(radian Angle)` | X rotation. | |
| `fromRotationY(radian Angle)` | Y rotation. | |
| `fromRotationZ(radian Angle)` | Z rotation. | |

## Setup Methods (Mutable)

| Method | Returns | Description |
|--------|---------|-------------|
| `setup(const fquat& rotation, const fvec3& scale)` | `fmat3_t&` | Set rotation + scale. |
| `setupIdentity()` | `fmat3_t&` | Set to identity. |
| `setupZero()` | `fmat3_t&` | Set to zero. |
| `setupRotation(const fquat& q)` | `fmat3_t&` | Set from quaternion. |
| `setupRotation(const radian3& euler)` | `fmat3_t&` | Set from Euler. |
| `setupScale(const fvec3& s)` | `fmat3_t&` | Set scale. |
| `setupScale(float s)` | `fmat3_t&` | Uniform scale. |

## Accessors

| Method | Returns | Description |
|--------|---------|-------------|
| `operator[](size_t row) const` | `fvec3` | Row as vector. |
| `operator()(size_t row, size_t col)` | `float&` | Mutable element (row-major feel). |
| `operator()(size_t row, size_t col) const` | `const float&` | Const element. |
| `operator std::span<const float,9>() const` | Span of elements. | Flat access (col-major). |
| `operator fquat() const` | `fquat` | To quaternion. |

## Operators

| Operator | Description | Example |
|----------|-------------|---------|
| `+ (const fmat3_t&)` | Add matrices. | `m1 + m2` |
| `- (const fmat3_t&)` | Subtract. | |
| `* (const fmat3_t&)` | Multiply (composition). | `m1 * m2` |
| `+= (const fmat3_t&)` | Add-assign. | |
| `-= (const fmat3_t&)` | Sub-assign. | |
| `*=(const fmat3_t&)` | Mul-assign. | |
| `* (const fvec3&)` | Transform vector. | `mat * v` |

## Math Functions

| Method | Returns | Description |
|--------|---------|-------------|
| `Transpose()` | `fmat3_t` | Transpose copy. |
| `Inverse()` | `fmat3_t` | Inverse copy (assume invertible). |
| `Determinant()` | `float` | Determinant. |
| `Orthogonalize()` | `fmat3_t&` | Orthogonalize in-place. |
| `Equals(const fmat3_t& other, float tolerance)` | `bool` | Approx equal. |

## Geometry Helpers

| Method | Returns | Description |
|--------|---------|-------------|
| `ExtractRotation()` | `fquat` | Extract quaternion. |
| `ExtractScale()` | `fvec3` | Extract scale. |
| `Forward()` | `fvec3` | Forward direction. |
| `Back()` | `fvec3` | Back direction. |
| `Up()` | `fvec3` | Up direction. |
| `Down()` | `fvec3` | Down direction. |
| `Left()` | `fvec3` | Left direction. |
| `Right()` | `fvec3` | Right direction. |
| `RotateVector(const fvec3& v)` | `fvec3` | Rotate vector. |
| `InvRotateVector(const fvec3& v)` | `fvec3` | Inverse rotate. |
| `TransformDirection(const fvec3& d)` | `fvec3` | Transform direction. |

## Mutable Chaining Methods (Post-Multiply)

| Method | Returns | Description |
|--------|---------|-------------|
| `Rotate(const fquat& q)` | `fmat3_t&` | Post-rotate. |
| `Rotate(const fvec3& axis, radian angle)` | `fmat3_t&` | Post-rotate axis-angle. |
| `RotateX(radian angle)` | `fmat3_t&` | Post-rotate X. |
| `RotateY(radian angle)` | `fmat3_t&` | Post-rotate Y. |
| `RotateZ(radian angle)` | `fmat3_t&` | Post-rotate Z. |
| `Scale(const fvec3& s)` | `fmat3_t&` | Post-scale. |
| `Scale(float s)` | `fmat3_t&` | Uniform post-scale. |

## Pre-Multiply Mutable Chaining

| Method | Returns | Description |
|--------|---------|-------------|
| `PreRotate(const fquat& q)` | `fmat3_t&` | Pre-rotate. |
| `PreRotate(const fvec3& axis, radian angle)` | `fmat3_t&` | Pre-rotate axis-angle. |
| `PreRotateX(radian angle)` | `fmat3_t&` | Pre-rotate X. |
| `PreRotateY(radian angle)` | `fmat3_t&` | Pre-rotate Y. |
| `PreRotateZ(radian angle)` | `fmat3_t&` | Pre-rotate Z. |
| `PreScale(const fvec3& s)` | `fmat3_t&` | Pre-scale. |
| `PreScale(float s)` | `fmat3_t&` | Uniform pre-scale. |

## Clear Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `ClearRotation()` | `fmat3_t&` | Reset rotation to identity. |
| `ClearScale()` | `fmat3_t&` | Reset scale to 1. |

## Immutable Versions (Copy Suffix)

| Method | Returns | Description |
|--------|---------|-------------|
| `RotateCopy(const fquat& q)` | `fmat3_t` | Post-rotate copy. |
| `RotateCopy(const fvec3& axis, radian angle)` | `fmat3_t` | Post-rotate axis-angle copy. |
| `ScaleCopy(const fvec3& s)` | `fmat3_t` | Post-scale copy. |
| `PreRotateCopy(const fquat& q)` | `fmat3_t` | Pre-rotate copy. |
| `PreRotateCopy(const fvec3& axis, radian angle)` | `fmat3_t` | Pre-rotate axis-angle copy. |
| `PreScaleCopy(const fvec3& s)` | `fmat3_t` | Pre-scale copy. |

## Safety and Validation

| Method | Returns | Description |
|--------|---------|-------------|
| `isFinite()` | `bool` | All elements finite. |
| `isIdentity()` | `bool` | Exact identity. |
| `SanityCheck()` | `void` | Assert finite (debug). |

## Tutorials and Examples

### Basic Matrix Creation and Transform
Create rotation matrix:
```cpp
fmat3 mat = fmat3::fromRotationX(xmath::pi_over2_v);  // 90° X
fvec3 transformed = mat * fvec3::fromUnitY();  // (0,1,0) to (0,0,1)
mat.Transpose();  // Transpose in-place
```

### Chaining Operations
Build complex transform:
```cpp
fmat3 mat;
mat.setupIdentity().RotateX(xmath::pi_over4_v).Scale(2.0f);  // Rotate then scale
fquat rot = mat.ExtractRotation();  // To quat
fvec3 scale = mat.ExtractScale();  // (2,2,2)
```

### Pre/Post Multiply
Compose rotations:
```cpp
mat.Rotate(fquat::fromAxisAngle(fvec3::fromUnitY(), xmath::pi_over2_v));  // Post
mat.PreRotate(fquat::fromAxisAngle(fvec3::fromUnitX(), xmath::pi_over4_v));  // Pre
fmat3 copy = mat.RotateCopy(fquat());  // Immutable
```

### Validation and Directions
Check and extract:
```cpp
if (mat.isFinite()) {
    fvec3 fwd = mat.Forward();  // [0,0,1] for identity
    mat.ClearRotation();  // Reset rotation
}
mat.SanityCheck();  // Assert in debug
```

### Advanced: Inverse and Orthogonalize
For physics:
```cpp
fmat3 inv = mat.Inverse();  // Inverse copy
mat.Orthogonalize();  // Make orthogonal
float det = mat.Determinant();  // 1 for rotation
```

## Related Classes
- `fquat`: For rotations used in matrix construction.
- `radian3`: For Euler angles.
- `fvec3`: For vectors transformed by matrix.
- `fmat4_t`: For 4x4 affine transforms.

---
