# `fmat4_t` Class Documentation

## Overview

The `fmat4_t` class is a templated 4x4 matrix implementation in the `xmath` namespace, using single-precision floating-point (`float`) 
components in column-major order. It supports two specializations:
- **SIMD Variant (`fmat4`)**: Optimized with SSE (`floatx4` columns, 4x4 aligned to 16 bytes), for faster transformations in graphics/physics.
- **CPU Variant (`fmat4d`)**: Compact 4x4 layout without padding, for memory efficiency or non-SSE platforms.

Designed for affine 3D transformations (translation, rotation, scale, projection) using homogeneous coordinates. Assumes right-hand rule for rotations. 
Prioritizes performance with mutable (in-place, chainable) and immutable (Copy-suffixed) operations. Uninitialized by default; use `fromIdentity()` 
or setup methods. Internal storage is column-major, but accessors (`()(row, col)`) feel row-major for usability.

### Key Features
- **Column-Major Storage**: Compatible with OpenGL/DirectX; unions for flexible access (columns, cells, elements, named m_rowcol).
- **SIMD/CPU Toggle**: Template bool for optimization vs. compactness; conversions between variants.
- **Construction**: From identity, zero, quaternion, Euler, translation/rotation/scale, perspective/ortho projections, look-at, billboard.
- **Operations**: Matrix/matrix mul, vector transform (fvec3/fvec4), transpose, inverse (general/SRT/RT variants), determinant.
- **Geometry Helpers**: Extract position/rotation/scale, direction vectors (forward/up/right), transform positions/directions.
- **Mutable Setup**: Chainable methods like `setupTranslation`, `Rotate`, `Scale`.
- **Validation**: `isFinite`, `isIdentity`, `SanityCheck` (asserts in debug).
- **Projection Support**: Perspective/ortho for cameras; look-at/billboard for graphics.
- **No Initialization**: Default ctor leaves uninitialized; explicit setup required.
- **Performance**: Inline/constexpr where possible; assumes finite matrices (asserts in debug).

Requires `xmath_flinear.h`. Column-major for graphics APIs.

### Comparison to Similar Classes
- Like GLM's `glm::mat4` but with SIMD specialization, mutable chaining (e.g., `Translate(t).Rotate(q)`), and specialized inverses 
      (`InverseSRT` for scale-rotation-translation). GLM is more general; this is optimized for 3D affine transforms.
- Similar to Eigen's `Eigen::Matrix4f` with better game features (e.g., `ExtractPosition` to fvec3, `fromLookAt`); SIMD variant rivals Eigen's but with easier extraction to quat/Euler.
- Engine analogs: Unreal's `FMatrix`, Godot's `Projection` or `Transform3D`, but standalone and templated.

### Performance Considerations
- **SIMD (`fmat4_t<true>`)**: Faster mul/transform (~4x on SSE); align instances to 16 bytes.
- **CPU (`fmat4_t<false>`)**: Smaller size (64 bytes vs padded); use for tight memory.
- **Chaining**: Mutable methods return `*this`: `mat.setupIdentity().Translate(t).Rotate(q).Scale(s);`.
- **Inverse Variants**: Use `InverseSRT`/`InverseRT` for faster inversion if matrix is SRT/RT.
- **Vector Transform**: `mat * fvec3` assumes w=1 (position); `mat * fvec4` for general.
- Tip: Convert CPU to SIMD for hot paths: `fmat4_t<true> simd_mat(cpu_mat);`.

### Quick Start Example
```cpp
#include "xmath_flinear.h"  // Assumes includes

int main() {
    xmath::fmat4 mat = xmath::fmat4::fromIdentity();
    mat.Translate(fvec3(1,0,0)).RotateX(xmath::pi_over2_v);  // Translate then rotate
    fvec3 transformed = mat * fvec3(0,1,0);  // Position transform
    fquat rot = mat.ExtractRotation();  // To quat

    std::cout << mat(3,0) << std::endl;  // Access m03 (translation x)
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `fmat4_t()` | Default (uninitialized). | `fmat4 mat;` |
| `fmat4_t(float diagonal)` | Diagonal matrix. | `fmat4(2.0f);` // Scale 2, w=1 |
| `fmat4_t(const std::array<float, 16>& arr)` | From array (col-major). | `std::array<float,16> a = {...}; fmat4(a);` |
| `fmat4_t(std::span<const float, 16> span)` | From span (col-major). | Similar to array. |
| `fmat4_t(const fquat& q)` | From quaternion rotation. | `fmat4(quat);` |
| `fmat4_t(const radian3& Euler)` | From Euler angles (ZXY). | `fmat4(euler);` |
| `fmat4_t(const fvec3& translation, const fquat& rotation, const fvec3& scale)` | From TRS. | `fmat4(t, r, s);` |
| `explicit fmat4_t(const fmat4_t<!T_USE_SIMD_V>& other)` | Convert SIMD/CPU. | `fmat4(cpu_mat);` |

## Static Constructors

| Method | Returns | Description |
|--------|---------|-------------|
| `fromIdentity()` | Identity matrix. | No transform. |
| `fromZero()` | Zero matrix. | All zeros. |
| `fromTranslation(const fvec3& t)` | Translation matrix. | |
| `fromRotation(const fquat& q)` | Rotation matrix. | |
| `fromRotation(const fvec3& axis, radian angle)` | Axis-angle rotation. | |
| `fromRotation(const radian3& Euler)` | Euler rotation (ZXY). | |
| `fromScale(const fvec3& s)` | Scale matrix. | |
| `fromPerspective(radian fov, float aspect, float near_plane, float far_plane)` | Perspective projection. | Graphics camera. |
| `fromPerspective(float left, float right, float bottom, float top, float near_plane, float far_plane)` | Custom perspective. | |
| `fromOrtho(float left, float right, float bottom, float top, float near_plane, float far_plane)` | Orthographic projection. | |
| `fromOrtho(float width, float height, float near_plane, float far_plane)` | Ortho by size. | |
| `fromLookAt(const fvec3& eye, const fvec3& target, const fvec3& up)` | Look-at matrix. | Camera view. |
| `fromBillboard(const fvec3& from, const fvec3& to, const fvec3& up)` | Billboard matrix. | Facing camera. |
| `fromRotationX(radian Angle)` | X rotation. | |
| `fromRotationY(radian Angle)` | Y rotation. | |
| `fromRotationZ(radian Angle)` | Z rotation. | |

## Setup Methods (Mutable)

| Method | Returns | Description |
|--------|---------|-------------|
| `setup(const fvec3& translation, const fquat& rotation, const fvec3& scale)` | `fmat4_t&` | Set TRS. |
| `setupIdentity()` | `fmat4_t&` | Set identity. |
| `setupZero()` | `fmat4_t&` | Set zero. |
| `setupTranslation(const fvec3& t)` | `fmat4_t&` | Set translation. |
| `setupRotation(const fquat& q)` | `fmat4_t&` | Set rotation. |
| `setupRotation(const radian3& euler)` | `fmat4_t&` | Set Euler rotation. |
| `setupScale(const fvec3& s)` | `fmat4_t&` | Set scale. |
| `setupScale(float s)` | `fmat4_t&` | Uniform scale. |

## Accessors

| Method | Returns | Description |
|--------|---------|-------------|
| `operator[](size_t row) const` | `fvec4` | Row as vector. |
| `operator()(size_t row, size_t col)` | `float&` | Mutable element. |
| `operator()(size_t row, size_t col) const` | `const float&` | Const element. |
| `operator std::span<const float,16>() const` | Span of elements. | Flat access (col-major). |
| `operator fquat() const` | `fquat` | To quaternion (rotation). |

## Operators

| Operator | Description | Example |
|----------|-------------|---------|
| `+ (const fmat4_t&)` | Add matrices. | `m1 + m2` |
| `- (const fmat4_t&)` | Subtract. | |
| `* (const fmat4_t&)` | Multiply (composition). | `m1 * m2` |
| `+= (const fmat4_t&)` | Add-assign. | |
| `-= (const fmat4_t&)` | Sub-assign. | |
| `*=(const fmat4_t&)` | Mul-assign. | |
| `* (const fvec4&)` | Transform vector4. | `mat * v4` |
| `* (const fvec3&)` | Transform vector3 (w=1). | `mat * v3` |
| `Equals(const fmat4_t& other, float tolerance)` | `bool` | Approx equal. |

## Math Functions

| Method | Returns | Description |
|--------|---------|-------------|
| `Transpose()` | `fmat4_t` | Transpose copy. |
| `Inverse()` | `fmat4_t` | General inverse copy. |
| `InverseSRT()` | `fmat4_t` | SRT inverse copy (faster). |
| `InverseRT()` | `fmat4_t` | RT inverse copy. |
| `Determinant()` | `float` | Determinant. |
| `Orthogonalize()` | `fmat4_t&` | Orthogonalize in-place. |

## Geometry Helpers

| Method | Returns | Description |
|--------|---------|-------------|
| `ExtractPosition()` | `fvec3` | Translation component. |
| `ExtractRotation()` | `fquat` | Rotation quaternion. |
| `ExtractScale()` | `fvec3` | Scale factors. |
| `Forward()` | `fvec3` | Forward direction. |
| `Back()` | `fvec3` | Back direction. |
| `Up()` | `fvec3` | Up direction. |
| `Down()` | `fvec3` | Down direction. |
| `Left()` | `fvec3` | Left direction. |
| `Right()` | `fvec3` | Right direction. |
| `RotateVector(const fvec3& v)` | `fvec3` | Rotate vector. |
| `InvRotateVector(const fvec3& v)` | `fvec3` | Inverse rotate. |
| `TransformPosition(const fvec3& p)` | `fvec3` | Transform position. |
| `TransformDirection(const fvec3& d)` | `fvec3` | Transform direction. |

## Mutable Chaining Methods (Post-Multiply)

| Method | Returns | Description |
|--------|---------|-------------|
| `Translate(const fvec3& t)` | `fmat4_t&` | Post-translate. |
| `Rotate(const fquat& q)` | `fmat4_t&` | Post-rotate. |
| `Rotate(const fvec3& axis, radian angle)` | `fmat4_t&` | Post-rotate axis-angle. |
| `RotateX(radian angle)` | `fmat4_t&` | Post-rotate X. |
| `RotateY(radian angle)` | `fmat4_t&` | Post-rotate Y. |
| `RotateZ(radian angle)` | `fmat4_t&` | Post-rotate Z. |
| `Scale(const fvec3& s)` | `fmat4_t&` | Post-scale. |
| `Scale(float s)` | `fmat4_t&` | Uniform post-scale. |

## Pre-Multiply Mutable Chaining

| Method | Returns | Description |
|--------|---------|-------------|
| `PreTranslate(const fvec3& t)` | `fmat4_t&` | Pre-translate. |
| `PreRotate(const fquat& q)` | `fmat4_t&` | Pre-rotate. |
| `PreRotate(const fvec3& axis, radian angle)` | `fmat4_t&` | Pre-rotate axis-angle. |
| `PreRotateX(radian angle)` | `fmat4_t&` | Pre-rotate X. |
| `PreRotateY(radian angle)` | `fmat4_t&` | Pre-rotate Y. |
| `PreRotateZ(radian angle)` | `fmat4_t&` | Pre-rotate Z. |
| `PreScale(const fvec3& s)` | `fmat4_t&` | Pre-scale. |
| `PreScale(float s)` | `fmat4_t&` | Uniform pre-scale. |

## Clear Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `ClearTranslation()` | `fmat4_t&` | Reset translation to 0. |
| `ClearRotation()` | `fmat4_t&` | Reset rotation to identity. |
| `ClearScale()` | `fmat4_t&` | Reset scale to 1. |

## Immutable Versions (Copy Suffix)

| Method | Returns | Description |
|--------|---------|-------------|
| `TranslateCopy(const fvec3& t)` | `fmat4_t` | Post-translate copy. |
| `RotateCopy(const fquat& q)` | `fmat4_t` | Post-rotate copy. |
| `RotateCopy(const fvec3& axis, radian angle)` | `fmat4_t` | Post-rotate axis-angle copy. |
| `ScaleCopy(const fvec3& s)` | `fmat4_t` | Post-scale copy. |
| `PreTranslateCopy(const fvec3& t)` | `fmat4_t` | Pre-translate copy. |
| `PreRotateCopy(const fquat& q)` | `fmat4_t` | Pre-rotate copy. |
| `PreRotateCopy(const fvec3& axis, radian angle)` | `fmat4_t` | Pre-rotate axis-angle copy. |
| `PreScaleCopy(const fvec3& s)` | `fmat4_t` | Pre-scale copy. |

## Safety and Validation

| Method | Returns | Description |
|--------|---------|-------------|
| `isFinite()` | `bool` | All elements finite. |
| `isIdentity()` | `bool` | Exact identity. |
| `SanityCheck()` | `void` | Assert finite (debug). |

## Tutorials and Examples

### Basic Matrix Creation and Transform
Create view matrix:
```cpp
fmat4_t<true> view = fmat4_t<true>::fromLookAt(eye, target, up);
fvec4 transformed = view * fvec4(pos, 1.0f);  // Homogeneous transform
fvec3 dir = view.TransformDirection(vec);  // Direction
```

### Chaining Operations
Build model matrix:
```cpp
fmat4 model;
model.setupIdentity().Translate(pos).Rotate(q).Scale(s);  // TRS order
fquat rot = model.ExtractRotation();  // To quat
```

### Pre/Post Multiply
Compose matrices:
```cpp
model.Rotate(fquat());  // Post-rotate
model.PreTranslate(offset);  // Pre-translate
fmat4 copy = model.RotateCopy(fquat());  // Immutable
```

### Projection for Graphics
Camera projection:
```cpp
fmat4 proj = fmat4::fromPerspective(fov, aspect, near, far);
fmat4 viewproj = proj * view;  // Compose
model.ClearTranslation();  // Reset position
```

### Advanced: Inverse and Orthogonalize
For inverse transform:
```cpp
fmat4 inv = mat.InverseSRT();  // Faster for SRT
mat.Orthogonalize();  // Make orthogonal
float det = mat.Determinant();  // 1 for rotation/scale
mat.SanityCheck();  // Debug assert
```

## Related Classes
- `fquat`: For rotations in matrix construction.
- `radian3`: For Euler angles.
- `fvec3`, `fvec4`: For vectors transformed by matrix.
- `fmat3_t`: For 3x3 (no translation).

---