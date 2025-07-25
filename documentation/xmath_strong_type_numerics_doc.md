# `xmath_strong_type_numberics.h` Documentation

## Overview

The `xmath_strong_type_numberics.h` header in the `xmath` namespace defines a CRTP (Curiously Recurring Template Pattern) base class `strong_typing_numerics_t` 
for creating strongly-typed numeric wrappers. This template enables type-safe arithmetic on custom numeric types (e.g., `radian`, `degree`) while inheriting 
standard operator behavior from the underlying numeric type `T` (e.g., `float`, `double`, `int`).

The class promotes safety by preventing implicit mixing of unrelated numeric types (e.g., adding radians to degrees), while supporting common arithmetic, 
compound assignment, unary negation (for signed types), and comparisons via the spaceship operator. It is designed for performance-critical applications 
like math libraries in games or simulations, where type errors can be costly.

### Key Features
- **CRTP Design**: Derived classes (e.g., `radian`) inherit operators without virtual overhead.
- **Type Safety**: Wrappers like `radian` ensure operations are only between compatible types.
- **Arithmetic Support**: Binary operators (+, -, *, /) with scalars and same-type instances; compound assignments (+=, -=, *=, /=).
- **Unary Negation**: For signed `T` (e.g., `float`, not `unsigned int`).
- **Comparisons**: Full ordering via `<=>` (spaceship operator), enabling ==, !=, <, >, <=, >=.
- **Zero Factory**: Static `fromZero()` for convenience.
- **Performance**: All methods/operators are `constexpr noexcept` where possible, header-only.

Requires no external dependencies. Used in other `xmath` headers (e.g., for `radian`, `degree` in trigonometry). **See Also**: `xmath_trigonometry.h` for examples of derived types.

### Comparison to Similar Utilities
- Unlike Boost's `strong_typedef`, this is lighter and C++20-focused with spaceship operator and constexpr support.
- Compared to Eigen's numeric types, it's more generic (CRTP for any numeric `T`) and avoids matrix dependencies, but lacks Eigen's vectorized operations.
- Engine analogs: Similar to Unreal's typed numerics (e.g., `TAngle`), but standalone and templated for broader use.

### Usage Notes
- **Derivation**: Derive from `strong_typing_numerics_t<Derived, T>` and add custom methods (e.g., conversions).
- **Signed/Unsigned**: Unary `-` requires signed `T` (compile-time check).
- **Performance**: Inline and constexpr; suitable for compile-time evaluation.
- **Type Safety**: Operators only work with same `T_DERIVED`; scalars are raw `T`.

### Quick Start Example
```cpp
#include "xmath_strong_type_numberics.h"

// Example derived type (similar to radian)
template <typename T>
struct custom_angle : xmath::strong_typing_numerics_t<custom_angle<T>, T> {
    using base = xmath::strong_typing_numerics_t<custom_angle<T>, T>;
    using base::base;  // Inherit constructors
};

int main() {
    custom_angle<float> a(1.0f), b(2.0f);
    auto sum = a + b;  // 3.0f
    a *= 3.0f;  // a = 3.0f
    if (a > b) { /* True */ }
    auto zero = custom_angle<float>::fromZero();  // 0.0f

    return 0;
}
```

## Template Parameters
- `T_DERIVED`: The derived class type (CRTP pattern).
- `T`: Underlying numeric type (e.g., `float`, `double`, must support arithmetic).

## Constructors

| Constructor | Description | Example |
|-------------|-------------|---------|
| `strong_typing_numerics_t()` | Default (uninitialized `m_Value`). | `Derived d;` |
| `explicit strong_typing_numerics_t(T value)` | Construct from raw value. | `Derived d(5.0f);` |

**Notes**: Both `constexpr noexcept`. Use explicit to avoid implicit conversions.

## Static Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `fromZero()` | `T_DERIVED` | Zero value instance. |

**Example**:
```cpp
auto zero = Derived::fromZero();  // Derived{ T{0} }
```

## Member Variables
- `T m_Value`: Underlying numeric value (public access).

## Arithmetic Operators

| Operator | Returns | Description | Example |
|----------|---------|-------------|---------|
| `T_DERIVED operator+(const T_DERIVED rhs) const` | `T_DERIVED` | Addition. | `a + b` |
| `T_DERIVED operator-(const T_DERIVED rhs) const` | `T_DERIVED` | Subtraction. | `a - b` |
| `T_DERIVED operator*(T scalar) const` | `T_DERIVED` | Scalar multiply (right). | `a * 2.0f` |
| `friend T_DERIVED operator*(T scalar, T_DERIVED rhs)` | `T_DERIVED` | Scalar multiply (left). | `2.0f * a` |
| `T_DERIVED operator/(T scalar) const` | `T_DERIVED` | Scalar divide. | `a / 2.0f` |

**Notes**: All `constexpr noexcept`. Only same-type or scalar operations.

## Compound Assignment Operators

| Operator | Returns | Description | Example |
|----------|---------|-------------|---------|
| `T_DERIVED& operator+=(const T_DERIVED rhs)` | `T_DERIVED&` | Add-assign. | `a += b` |
| `T_DERIVED& operator-=(const T_DERIVED rhs)` | `T_DERIVED&` | Subtract-assign. | `a -= b` |
| `T_DERIVED& operator*=(T scalar)` | `T_DERIVED&` | Multiply-assign. | `a *= 2.0f` |
| `T_DERIVED& operator/=(T scalar)` | `T_DERIVED&` | Divide-assign. | `a /= 2.0f` |

**Notes**: All `constexpr noexcept`, chainable.

## Unary Operators

| Operator | Returns | Description | Example |
|----------|---------|-------------|---------|
| `T_DERIVED operator-() const` | `T_DERIVED` | Negation (signed `T` only). | `-a` |

**Notes**: `constexpr noexcept`, requires `std::is_unsigned_v<T> == false`.

## Comparison Operators

| Operator | Description | Example |
|----------|-------------|---------|
| `auto operator<=>(const strong_typing_numerics_t& rhs) const` | Three-way comparison. | `a <=> b` |

**Notes**: `constexpr noexcept`, enables ==, !=, <, >, <=, >= via spaceship.

## Tutorials and Examples

### Creating a Derived Type
Define a custom angle type:
```cpp
template <typename T>
struct angle : xmath::strong_typing_numerics_t<angle<T>, T> {
    using base = xmath::strong_typing_numerics_t<angle<T>, T>;
    using base::base;
    // Add custom methods if needed
};
angle<float> a(45.0f), b(90.0f);
angle<float> sum = a + b;  // 135.0f
sum -= angle<float>(30.0f);  // 105.0f
```

### Arithmetic Chaining
Chain operations for physics calculations:
```cpp
Derived v(10.0f);
v += Derived(5.0f) * 2.0f;  // v = 10 + (5 * 2) = 20
v /= 4.0f;  // v = 5.0f
if (v == Derived(5.0f)) { /* True */ }
```

### Signed vs Unsigned Types
For signed `float`:
```cpp
Derived neg = -Derived(5.0f);  // -5.0f
```
For unsigned (e.g., `unsigned int`):
```cpp
// -Derived(u) compiles only if T signed
```

### Comparisons in Algorithms
Sort or check ranges:
```cpp
if (a < b) { /* If a.m_Value < b.m_Value */ }
bool eq = (a <=> b) == 0;  // Equivalent to ==
```

### Error Handling and Safety
- **Type Mismatch**: Won't compile if mixing unrelated derived types.
- **Division by Zero**: Relies on underlying `T` behavior (e.g., float NaN).
- **Overflow**: Handled by `T` (e.g., integer wraparound).
- **Comparisons**: Spaceship handles NaN/Inf per `T`.

## Performance Tips
- **Constexpr**: Use in compile-time contexts:
  ```cpp
  constexpr Derived c = Derived(1.0f) + Derived(2.0f);  // Computed at compile-time
  ```
- **No Overhead**: CRTP avoids virtual dispatch; operators inline.
- **Zero Initialization**: Use `fromZero()` for constants.
- **Derived Extensions**: Add custom ops in derived classes (e.g., wrapping for angles).

## Related Classes (Brief Context)
- **`radian_t<T>`**, **`degree_t<T>`**: Derived from this for angles in `xmath_trigonometry.h`.
- **`radian3`**: Uses `radian` for Euler angles.
- **`fquat`**, **`fmat3`**, **`fmat4`**: May use typed numerics for safety.

---
