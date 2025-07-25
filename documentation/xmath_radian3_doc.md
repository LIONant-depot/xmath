# `radian3` Class Documentation

## Overview

The `radian3` struct in the `xmath` namespace represents a 3D rotation using Euler angles in the `radian` type, a specialized wrapper for `float` values representing angles in radians. The angles are:
- **Pitch**: Rotation around the X-axis.
- **Yaw**: Rotation around the Y-axis.
- **Roll**: Rotation around the Z-axis.

The rotation order is ZXY (Roll, Pitch, Yaw), suitable for applications like aerospace, robotics, and game development (e.g., heading-attitude-bank). This class is lightweight, using three `radian` values (internally `float`) without SIMD or special alignment, prioritizing minimal memory footprint over computational speed.

Euler angles are intuitive but prone to **gimbal lock** (loss of a degree of freedom when pitch is near ±π/2). `radian3` provides utilities like `hasGimbalLock()` and `MinAngleDiff()` to mitigate such issues, making it ideal for scenarios where memory efficiency or human-readable angles are preferred over quaternion or matrix representations.

### Key Features
- **Compact Storage**: Three `radian` values (~12 bytes), no padding.
- **Angle Management**: Wrap angles to [-π, π] (`ModAngle()`) or [0, 2π] (`ModAngle2()`), compute minimal differences for smooth transitions.
- **Gimbal Lock Detection**: `hasGimbalLock()` flags problematic configurations.
- **Validation**: `isValid()` checks for finite angles; `isInrange()` enforces bounds.
- **Interpolation**: `Lerp()` for linear blending (use `fquat` for spherical interpolation).
- **Operators**: Component-wise arithmetic for combining rotations.
- **No Hardware Acceleration**: Scalar operations, no alignment beyond atomic types.

Requires `xmath.h`. For performance-critical paths, convert to `fquat` (quaternions) or `fmat4`/`fmat3` (matrices). **See Also**: `radian`, `fquat`, `fmat3`, `fmat4`.

### Comparison to Similar Classes
- Unlike GLM's `glm::vec3` for Euler angles, `radian3` uses the `radian` type with dedicated angle operations (e.g., `ModAngle()`, `MinAngleDiff()`), tailored for rotations.
- Compared to Eigen's `Eigen::AngleAxis`, `radian3` is lighter (no matrix dependencies) and game-focused with features like `hasGimbalLock()`. Eigen integrates better with linear algebra, but `radian3` excels in angle manipulation and memory efficiency.
- Engine analogs: Similar to Unity's `Transform.eulerAngles` (but using `radian`) or Unreal's `FRotator`, but standalone with advanced wrapping and difference calculations.

### Usage Notes
- **Radian Type**: `radian` is a `float` wrapper for angles. Construct with `radian(float)`, e.g., `radian(0.4f)`.
- **Pi Constant**: Use `xmath::pi_v` (type `radian`) for π, e.g., `radian3(xmath::pi_v, radian(0), radian(0))`.
- **Gimbal Lock**: Occurs when pitch is near ±π/2. Use `hasGimbalLock()` to detect and switch to `fquat` for robustness.
- **Performance**: Scalar operations are slower than SIMD-based classes like `fvec4`. Convert to `fmat4` for fast transformations in loops.

### Quick Start Example
```cpp
#include "xmath.h"  // Includes radian3, radian, and xmath::pi_v

int main() {
    xmath::radian3 rot(xmath::pi_v / radian(2.0f), radian(0.0f), radian(0.0f));  // 90° pitch
    if (rot.hasGimbalLock()) {
        // Convert to fquat (assumed external toFquat())
    }
    rot.ModAngle();  // Wrap to [-π, π]
    xmath::radian3 target(radian(0.0f), xmath::pi_v, radian(0.0f));  // 180° yaw
    xmath::radian3 interp = rot.Lerp(0.5f, target);  // Interpolate halfway

    std::cout << "Pitch: " << float(interp.m_Pitch) << " radians" << std::endl;
    return 0;
}
```

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `radian3()` | Default (uninitialized). | `radian3 r;` |
| `radian3(const radian Angles)` | All angles set to `Angles`. | `radian3(xmath::pi_v);` // [π, π, π] |
| `radian3(radian Pitch, radian Yaw, radian Roll)` | From individual angles. | `radian3(xmath::pi_v / radian(2), radian(0), radian(0));` // 90° pitch |

**Notes**:
- All `constexpr noexcept`. Initialize explicitly to avoid undefined values.
- Use `radian(float)` for float literals, e.g., `radian(0.4f)`, or `xmath::pi_v` for π.

## Static Properties

| Method | Returns | Description |
|--------|---------|-------------|
| `fromZero()` | `radian3(radian(0), radian(0), radian(0))` | Zero rotation. |

**Example**:
```cpp
auto zero = xmath::radian3::fromZero();  // [0, 0, 0]
```

## Instance Methods

| Method | Returns | Description | Mutable? |
|--------|---------|-------------|----------|
| `Zero()` | `void` | Set all angles to `radian(0)`. | Yes |
| `ModAngle()` | `radian3&` | Wrap angles to [-π, π] in-place. | Yes |
| `ModAngleCopy() const` | `radian3` | Wrapped copy [-π, π]. | No |
| `ModAngle2()` | `radian3&` | Wrap to [0, 2π] in-place. | Yes |
| `ModAngle2Copy() const` | `radian3` | Wrapped copy [0, 2π]. | No |
| `MinAngleDiff(const radian3& R) const` | `radian3` | Minimal angular differences, accounting for 2π periodicity. | No |
| `Difference(const radian3& R) const` | `radian` | Sum of absolute angle differences (no wrapping). | No |
| `isInrange(radian Min, radian Max) const` | `bool` | All angles in [Min, Max]. | No |
| `isValid() const` | `bool` | All angles finite (not NaN/Inf). | No |
| `hasGimbalLock() const` | `bool` | True if pitch near ±π/2 (gimbal lock risk). | No |
| `Normalize()` | `radian3&` | Normalize angles (calls `ModAngle()`). | Yes |
| `NormalizeCopy() const` | `radian3` | Normalized copy. | No |
| `Clamp(radian Min, radian Max)` | `radian3&` | Clamp each angle in-place. | Yes |
| `ClampCopy(radian Min, radian Max) const` | `radian3` | Clamped copy. | No |
| `Lerp(float t, const radian3& Other) const` | `radian3` | Linear interpolation (t=[0,1]). | No |

**Notes**:
- `ModAngle()` wraps to [-π, π]; `ModAngle2()` to [0, 2π]. Choose based on your application's angle convention.
- `MinAngleDiff()` ensures shortest rotation path, critical for smooth animations or AI.
- `Lerp()` interpolates linearly; for spherical interpolation, convert to `fquat`.
- `hasGimbalLock()` likely checks if pitch is within a small epsilon of ±π/2 (e.g., ±0.01 radians).

## Operator Overloads

| Operator | Description | Example |
|----------|-------------|---------|
| `== (const radian3& R) const` | Exact equality of angles. | `r1 == r2` |
| `+= (const radian3& R)` | Add angles component-wise. | `r += other;` |
| `-= (const radian3& R)` | Subtract angles. | |
| `*= (float Scalar)` | Scale angles (float, not radian). | `r *= 2.0f;` |
| `/= (float Scalar)` | Divide angles. | |
| `+ (const radian3& R) const` | Add copy. | `r1 + r2` |
| `- (const radian3& R) const` | Subtract copy. | |
| `- () const` | Negate copy. | `-r` |
| `/ (float S) const` | Divide copy. | `r / 2.0f` |
| `* (float S) const` | Scale copy. | `r * 2.0f` |
| `friend * (float S, const radian3& R)` | Left scale. | `2.0f * r` |

**Notes**:
- Operators do not automatically wrap angles; use `ModAngle()` or `ModAngle2()` after arithmetic.
- Scaling uses raw `float` (not `radian`) for simplicity; wrap results if needed.

## Member Variables
- `radian m_Pitch`: Rotation around X-axis.
- `radian m_Yaw`: Rotation around Y-axis.
- `radian m_Roll`: Rotation around Z-axis.

Direct access is allowed, but prefer methods for operations to maintain type safety and handle wrapping.

## Tutorials and Examples

### Basic Rotation Composition
Combine rotations for a camera orientation:
```cpp
xmath::radian3 cam(xmath::pi_v / radian(4.0f), radian(0), radian(0));  // 45° pitch
xmath::radian3 turn(radian(0), xmath::pi_v / radian(2.0f), radian(0));  // 90° yaw
cam += turn;  // Combine
cam.ModAngle();  // Wrap to [-π, π]
if (!cam.isValid()) {
    cam.Zero();  // Reset on invalid
}
```

### Handling Gimbal Lock
Mitigate gimbal lock in a flight simulator:
```cpp
xmath::radian3 plane(xmath::pi_v / radian(2.0f) + radian(0.01f), radian(0), radian(0));  // Near 90° pitch
if (plane.hasGimbalLock()) {
    // Convert to fquat (assumed external method, e.g., toFquat())
    // Use quaternion for rotation to avoid lock
}
xmath::radian3 safe = plane.ClampCopy(
    -xmath::pi_v / radian(2.0f) + radian(0.01f),
    xmath::pi_v / radian(2.0f) - radian(0.01f)
);  // Avoid lock
```

### Interpolation for Animation
Smoothly transition an object's rotation:
```cpp
xmath::radian3 start(radian(0), radian(0), radian(0));
xmath::radian3 end(xmath::pi_v, xmath::pi_v / radian(2.0f), radian(0));
for (float t = 0; t <= 1; t += 0.1f) {
    xmath::radian3 frame = start.Lerp(t, end);
    frame.Normalize();  // Wrap to [-π, π]
    // Convert to fmat4 for rendering (assumed toFmat4())
}
```

### Advanced: Minimal Angle Difference for AI
Smooth character turning to face a target:
```cpp
xmath::radian3 current(radian(0), radian(3.0f) * xmath::pi_v, radian(0));  // Overwrapped yaw (3π)
current.ModAngle2();  // Yaw = π
xmath::radian3 target(radian(0), xmath::pi_v / radian(2.0f), radian(0));  // 90° yaw
xmath::radian3 diff = current.MinAngleDiff(target);  // Yaw diff = -π/2
current += diff * 0.1f;  // Turn 10% per frame
current.ModAngle();  // Keep in [-π, π]
```

### Error Handling
- **Invalid Angles**: Use `isValid()` after user input or calculations to catch NaN/Inf:
  ```cpp
  if (!rot.isValid()) {
      rot = xmath::radian3::fromZero();  // Reset
  }
  ```
- **Gimbal Lock**: Check `hasGimbalLock()` in applications like flight sims or 3D modeling. Convert to `fquat` for robustness:
  ```cpp
  if (rot.hasGimbalLock()) {
      // Convert to fquat (needs fquat definition)
  }
  ```
- **Angle Wrapping**: Use `MinAngleDiff()` for smooth rotations to avoid 2π jumps:
  ```cpp
  xmath::radian3 diff = start.MinAngleDiff(end);
  start += diff * 0.05f;  // Gradual turn
  ```
- **Difference**: Use `Difference()` for raw error metrics (no periodicity):
  ```cpp
  radian error = rot.Difference(target);  // Sum of abs diffs
  ```

## Performance Tips
- **Memory**: Minimal footprint (~12 bytes), ideal for storage or serialization.
- **Compute**: Scalar operations are slower than SIMD classes (`fvec4`, `fquat`). For rendering loops, convert to `fmat4` once and reuse.
- **Caching**: Cache `ModAngle()` or `Normalize()` results in loops to avoid redundant wrapping.
- **Interpolation**: `Lerp()` is linear; for spherical interpolation, convert to `fquat` (assumes `fquat` has slerp).
- **Conversion**: If applying rotations to many vectors, convert to `fmat4` or `fmat3` for SIMD efficiency.

## Related Classes
- **`radian`**: Angle type (wraps `float`). Construct with `radian(float)`, e.g., `radian(0.4f)`.
- **`fquat`**: Quaternion for gimbal-free rotations. Use for slerp or robust transforms.
- **`fmat3`**: 3x3 matrix, often for rotation-only transforms.
- **`fmat4`**: 4x4 matrix for full transforms (rotation + translation).
- **`fvec3`**: For Euler angles as vectors, but lacks angle-specific ops like `ModAngle()`.

## Example with Visual Context
Consider a character facing yaw = 3π (overwrapped) needing to turn to yaw = π/2 for AI pathing. Without `MinAngleDiff()`, interpolation might take the long path (3π → π/2, crossing 2π). Using it:
```cpp
xmath::radian3 current(radian(0), radian(3.0f) * xmath::pi_v, radian(0));  // Yaw = 3π
current.ModAngle2();  // Yaw = π
xmath::radian3 target(radian(0), xmath::pi_v / radian(2.0f), radian(0));  // Yaw = π/2
xmath::radian3 diff = current.MinAngleDiff(target);  // Yaw diff = -π/2
current += diff * 0.1f;  // Turn -π/20 per frame
current.ModAngle();  // Stay in [-π, π]
```
---
