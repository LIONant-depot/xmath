# `fvec2` Class Documentation

## Overview

The `fvec2` class represents a 2D vector using single-precision floating-point numbers (`float`). It is part of the `xmath` namespace and is designed for high-performance mathematical operations in applications such as game development, graphics, simulations, and scientific computing. Unlike SIMD-optimized variants, `fvec2` prioritizes simplicity and portability without relying on vector extensions.

### Key Features
- **Uninitialized Default Construction**: For performance, default-constructed instances do not initialize memory. Always initialize explicitly if needed.
- **Mutable vs. Immutable Operations**: Methods that modify the vector in-place (mutable) have shorter names and return `*this` for chaining (e.g., `Normalize()`). Copying variants are suffixed with `_Copy` (e.g., `NormalizeCopy()`).
- **Constexpr Support**: Many operations are `constexpr` for compile-time evaluation, targeting C++20.
- **Validation**: Includes checks like `isFinite()` and assertions for finite components in debug builds.
- **Swizzling**: HLSL-style swizzle methods for component rearrangement, supporting up to `fvec4` returns.
- **Geometry and Math Focus**: Supports vector algebra, trigonometry, clamping, normalization, and geometric utilities like rotation and projection.
- **Interoperability**: Convertible to/from `std::array<double, 2>`, `std::span<float>`, and string representations.

This class is header-only and requires inclusion of `xmath_fvector.h`. It emphasizes performance over safety—use `const` qualifiers for immutability.

### Comparison to Similar Classes
- Similar to GLM's `glm::vec2` but with more component-wise math functions and safe normalization variants.
- Analogous to Eigen's `Eigen::Vector2f` but optimized for games with geometry-specific methods like `Perp()` and `SignedAngleBetween()`.
- Engine integrations: Comparable to Unreal's `FVector2D` or Godot's `Vector2`, but standalone and extensible.

### Quick Start Example
```cpp
#include <xmath_fvector.h>  // Assuming this includes the necessary headers

int main() {
    xmath::fvec2 v1(3.0f, 4.0f);  // Create a vector (3, 4)
    xmath::fvec2 v2 = v1.NormalizeCopy();  // Normalized copy: (0.6, 0.8)
    float len = v1.Length();  // 5.0f
    xmath::fvec2 projected = v1.ProjectCopy(v2);  // Project v1 onto v2

    std::cout << v1.ToString() << std::endl;  // Outputs: "[3, 4]"
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `fvec2()` | Default constructor (uninitialized for performance). | `fvec2 v;` |
| `fvec2(fvec2&&)` | Move constructor (noexcept). | N/A (implicit) |
| `fvec2(const fvec2&)` | Copy constructor (noexcept). | `fvec2 v2(v1);` |
| `fvec2(float x, float y)` | Constructs from individual components. | `fvec2(1.0f, 2.0f);` |
| `fvec2(float value)` | Constructs with all components set to `value`. | `fvec2(5.0f);` // [5, 5] |
| `fvec2(std::span<float> Span)` | Constructs from a span of floats (must have exactly 2 elements). | `std::array<float, 2> arr = {1,2}; fvec2 v(arr);` |
| `fvec2(const std::array<double,2>& Conversion)` | Constructs by converting from doubles. | `std::array<double, 2> d = {1.0, 2.0}; fvec2 v(d);` |

**Notes**: All constructors are `constexpr` and `noexcept` where possible. Ensure spans have exactly 2 elements to avoid runtime errors.

## Assignment and Conversion Operators

| Operator/Method | Description | Example |
|-----------------|-------------|---------|
| `fvec2& operator=(const fvec2&)` | Assignment operator (copy). | `v1 = v2;` |
| `operator std::array<double,2>() const` | Converts to array of doubles. | `std::array<double,2> arr = v;` |
| `operator std::string() const` | Implicit conversion to string (e.g., "[x, y]"). | `std::string s = v;` |
| `std::string ToString() const` | Explicit string conversion. | `v.ToString();` |
| `friend std::ostream& operator<<(std::ostream& os, const fvec2& vec)` | Stream output. | `std::cout << v;` |

## Static Properties

These are convenience factory methods for common vectors:

| Method | Returns | Description |
|--------|---------|-------------|
| `fromZero()` | `fvec2(0, 0)` | Zero vector. |
| `fromOne()` | `fvec2(1, 1)` | One vector. |
| `fromUnitX()` | `fvec2(1, 0)` | Unit X vector. |
| `fromUnitY()` | `fvec2(0, 1)` | Unit Y vector. |
| `fromUp()` | `fvec2(0, 1)` | Up direction (alias for Unit Y). |
| `fromDown()` | `fvec2(0, -1)` | Down direction. |
| `fromLeft()` | `fvec2(-1, 0)` | Left direction. |
| `fromRight()` | `fvec2(1, 0)` | Right direction (alias for Unit X). |
| `fromRandomUnitVector()` | Randomized unit vector. | Uses random generation; not `constexpr`. |

**Example**:
```cpp
auto zero = xmath::fvec2::fromZero();
auto random = xmath::fvec2::fromRandomUnitVector();  // e.g., (0.6, 0.8) with length 1
```

## Static Methods

These operate on two vectors without requiring an instance:

| Method | Returns | Description |
|--------|---------|-------------|
| `Dot(const fvec2& a, const fvec2& b)` | `float` | Dot product: `a.x*b.x + a.y*b.y`. |
| `Min(const fvec2& a, const fvec2& b)` | `fvec2` | Component-wise minimum. |
| `Max(const fvec2& a, const fvec2& b)` | `fvec2` | Component-wise maximum. |
| `Lerp(const fvec2& a, const fvec2& b, float t)` | `fvec2` | Linear interpolation: `a + t*(b - a)`. |
| `Distance(const fvec2& a, const fvec2& b)` | `float` | Euclidean distance. |
| `Cross(const fvec2& a, const fvec2& b)` | `float` | 2D cross product (scalar): `a.x*b.y - a.y*b.x`. |

These have instance overloads (e.g., `this->Dot(a)`) for convenience.

**Example**:
```cpp
fvec2 a(1, 2), b(3, 4);
float dot = fvec2::Dot(a, b);  // 11
fvec2 lerp = fvec2::Lerp(a, b, 0.5f);  // (2, 3)
```

## Instance Methods - Basic Operations

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `Length() const` | `float` | Euclidean length (sqrt(x² + y²)). | No |
| `LengthSq() const` | `float` | Squared length (faster for comparisons). | No |
| `NormalizeCopy() const` | `fvec2` | Normalized copy (divides by length). | No |
| `Normalize()` | `fvec2&` | Normalizes in-place. | Yes |
| `NormalizeSafeCopy() const` | `fvec2` | Safe normalize (handles zero length). | No |
| `NormalizeSafe()` | `fvec2&` | Safe normalize in-place. | Yes |
| `LimitLengthCopy(float MaxLength) const` | `fvec2` | Caps length to `MaxLength`. | No |
| `isFinite() const` | `bool` | Checks if all components are finite. | No |
| `isInRange(float min, float max) const` | `bool` | Checks if all components in [min, max]. | No |
| `Equals(const fvec2& other, float tolerance) const` | `bool` | Approximate equality with tolerance. | No |

**Normalization Example**:
```cpp
fvec2 v(3, 4);
float len = v.Length();  // 5
v.Normalize();  // Now (0.6, 0.8)
if (!v.NormalizeSafe()) { /* Handle zero vector */ }
```

## Instance Methods - Component-Wise Math

These apply functions component-wise. Copy versions return new vectors; mutable ones modify in-place.

| Method (Copy/Mutable) | Description | Example |
|-----------------------|-------------|---------|
| `AbsCopy() / Abs()` | Absolute value. | `fvec2(-1, 2).AbsCopy() == (1, 2)` |
| `OneOverCopy() / OneOver()` | Reciprocal (1/x, 1/y). | Handles division by zero? Use with care. |
| `SqrtCopy() / Sqrt()` | Square root. | `fvec2(4, 9).SqrtCopy() == (2, 3)` |
| `InvSqrtCopy() / InvSqrt()` | Inverse square root (1/sqrt). | Useful for normalization. |
| `SignCopy() / Sign()` | Sign (-1, 0, 1). | `fvec2(-3, 0).SignCopy() == (-1, 0)` |
| `FloorCopy() / Floor()` | Floor. | `fvec2(1.9, -1.1).FloorCopy() == (1, -2)` |
| `CeilCopy() / Ceil()` | Ceiling. | |
| `FractCopy() / Fract()` | Fractional part. | `fvec2(1.9, -1.1).FractCopy() == (0.9, 0.9)` (note: handles negatives). |
| `RoundCopy() / Round()` | Round to nearest. | |
| `TruncCopy() / Trunc()` | Truncate towards zero. | |
| `ModCopy(float divisor) / Mod(float divisor)` | Modulo. | `fvec2(5, -5).ModCopy(3) == (2, 1)` (C++ % behavior for negatives). |
| `ClampCopy(float min_val, float max_val) / Clamp(float min_val, float max_val)` | Clamp to scalar range. | |
| `ClampCopy(const fvec2& min, const fvec2& max) / Clamp(const fvec2& min, const fvec2& max)` | Component-wise clamp. | |
| `Step(float edge) const` | Step function (0 if < edge, 1 otherwise). | Copy only. |
| `SmoothStep(float edge0, float edge1) const` | Smooth Hermite interpolation. | Copy only. |
| `LogCopy() / Log()` | Natural log. | |
| `Log2Copy() / Log2()` | Log base 2. | |
| `PowCopy(float exp) / Pow(float exp)` | Raise to power. | |
| `SinCopy() / Sin()` | Sine (radians). | |
| `CosCopy() / Cos()` | Cosine. | |
| `TanCopy() / Tan()` | Tangent. | |
| `AsinCopy() / Asin()` | Arc sine. | |
| `AcosCopy() / Acos()` | Arc cosine. | |
| `AtanCopy() / Atan()` | Arc tangent. | |
| `Atan2Copy(const fvec2& x) / Atan2(const fvec2& x)` | Component-wise atan2(y, x). | Note: `this` is y. |
| `SignedAngleBetween(const fvec2& v) const` | `radian` | Signed angle between vectors. |
| `ExpCopy() / Exp()` | Exponential (e^x). | |

**Trigonometry Example**:
```cpp
fvec2 angles(0.0f, XMATH_PI / 2.0f);
fvec2 sins = angles.SinCopy();  // (0, 1)
```

## Instance Methods - Geometry

| Method | Returns | Description |
|--------|---------|-------------|
| `Reflection(const fvec2& normal) const` | `fvec2` | Reflects over normal (assumes unit normal). |
| `DistanceSquare(const fvec2& v) const` | `float` | Squared distance to `v`. |
| `AngleBetween(const fvec2& v) const` | `radian` | Unsigned angle between vectors. |
| `GridSnap(float gridX, float gridY)` | `fvec2&` | Snaps to grid (mutable). |
| `Perp() const` | `fvec2` | Perpendicular vector (-y, x). |
| `WhichSideOfLine(const fvec2& V0, const fvec2& V1) const` | `float` | Side of line V0-V1 (>0 left, <0 right, 0 on line). |
| `ClosestPointInLine(const fvec2& V0, const fvec2& V1) const` | `fvec2` | Closest point on infinite line. |
| `ClosestPointInLineSegment(const fvec2& V0, const fvec2& V1) const` | `fvec2` | Closest point on segment. |
| `RotateCopy(radian angle) const` | `fvec2` | Rotated copy. |
| `Rotate(radian angle)` | `fvec2&` | Rotate in-place. |
| `ProjectCopy(const fvec2& onto) const` | `fvec2` | Projection onto vector. |
| `Project(const fvec2& onto)` | `fvec2&` | Project in-place. |

**Geometry Example**:
```cpp
fvec2 dir(1, 0);
fvec2 rotated = dir.RotateCopy(XMATH_PI / 2.0f);  // (0, 1)
fvec2 normal(0, 1);
fvec2 reflected = dir.Reflection(normal);  // (1, 0) unchanged
```

## Swizzle Methods

Swizzles return copies with rearranged components. Supported for scalar, `fvec2`, `fvec3`, `fvec4`.

- Scalars: `x()`, `y()`
- `fvec2`: `xx()`, `xy()`, `yx()`, `yy()`
- `fvec3`: `xxx()`, `xxy()`, `xyx()`, `xyy()`, `yxx()`, `yxy()`, `yyx()`, `yyy()`
- `fvec4`: `xxxx()`, `xxxy()`, ..., `yyyy()` (all combinations)

**Example**:
```cpp
fvec2 v(1, 2);
fvec2 swapped = v.yx();  // (2, 1)
fvec3 repeated = v.xxx();  // (1, 1, 1)
fvec4 extended = v.xxyy();  // (1, 1, 2, 2)
```

**Notes**: All `constexpr` where possible. Useful for shader-like code or constructing higher-dimensional vectors.

## Operator Overloads

| Operator | Description | Example |
|----------|-------------|---------|
| `+ (const fvec2& other) const` | Addition. | `v1 + v2` |
| `- (const fvec2& other) const` | Subtraction. | |
| `* (float scalar) const` | Scalar multiply. | `v * 2.0f` |
| `/ (float scalar) const` | Scalar divide. | |
| `+= (const fvec2& other)` | Add-assign. | |
| `-= (const fvec2& other)` | Subtract-assign. | |
| `*= (float scalar)` | Multiply-assign. | |
| `/= (float scalar)` | Divide-assign. | |
| `== (const fvec2& other) const` | Exact equality. | Use `Equals` for tolerance. |
| `!= (const fvec2& other) const` | Inequality. | |
| `[] (std::int32_t index) const` | Const access (0=x, 1=y). | `v[0]` |
| `[] (std::int32_t index)` | Mutable access. | `v[0] = 3.0f;` |
| `friend * (float scalar, const fvec2& v)` | Left scalar multiply. | `2.0f * v` |
| `friend - (const fvec2& v)` | Unary negation. | `-v` |

**Operator Example**:
```cpp
fvec2 v(1, 2);
v *= 3.0f;  // (3, 6)
if (v == fvec2(3, 6)) { /* True */ }
```

## Tutorials and Examples

### Basic Vector Algebra
Vectors can represent positions, directions, or colors. Here's a simple physics simulation snippet:
```cpp
fvec2 position(0, 0);
fvec2 velocity(1, 0.5f);
float deltaTime = 0.1f;

position += velocity * deltaTime;  // Update position
velocity.Normalize();  // Keep direction unit
```

### Geometric Operations
For 2D collision detection:
```cpp
fvec2 point(5, 5);
fvec2 lineStart(0, 0), lineEnd(10, 0);
fvec2 closest = point.ClosestPointInLineSegment(lineStart, lineEnd);  // (5, 0)
float side = point.WhichSideOfLine(lineStart, lineEnd);  // Positive if above
```

### Advanced: Rotation and Projection
Simulate orbiting:
```cpp
fvec2 center(0, 0);
fvec2 satellite(10, 0);
radian angle = XMATH_PI / 180.0f;  // 1 degree
satellite.Rotate(angle);  // Rotate around origin
fvec2 projected = satellite.ProjectCopy(fvec2(1, 0));  // Project onto x-axis
```

### Performance Tips
- Prefer mutable methods for chaining: `v.Abs().Normalize().Clamp(0, 1);`
- Use `LengthSq()` over `Length()` for comparisons to avoid sqrt.
- For random vectors, `fromRandomUnitVector()` uses std random; seed if needed.

### Error Handling
- Methods like `Normalize()` assume non-zero length; use safe variants for robustness.
- Trigonometric functions follow std::math behavior (e.g., domain errors for asin >1).

## Related Classes
- `fvec3`, `fvec4`: Higher-dimensional analogs.
- `radian`: Angle type used in rotations.
- See `xmath` namespace for matrices and quaternions.

---
