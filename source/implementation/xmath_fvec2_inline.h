#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
    #error "You must include xmath_vector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------
    //
    // Constructs a vector with specified x and y components.
    //
    // Parameters:
    //  x - The x component.
    //  y - The y component.
    //
    constexpr fvec2::fvec2(float x, float y) noexcept
        : m_X(x), m_Y(y)
    {
    }

    //------------------------------------------------------------------------------
    // Constructs a vector with all components set to the same value.
    //
    // Parameters:
    //  value - The value for all components.
    //
    constexpr fvec2::fvec2(float value) noexcept
        : m_X(value), m_Y(value)
    {
    }

    //------------------------------------------------------------------------------
    // Constructs a vector from a span of floats (at least 2 elements).
    //
    // Parameters:
    //  Span - The span containing at least x and y values.
    //
    constexpr fvec2::fvec2(std::span<float> Span) noexcept
    {
        if (Span.size() >= 2)
        {
            m_X = Span[0];
            m_Y = Span[1];
        }
    }

    //------------------------------------------------------------------------------
    // Constructs a vector from an array of doubles, casting to float.
    //
    // Parameters:
    //  Conversion - The array with x and y values.
    //
    constexpr fvec2::fvec2(const std::array<double, 2>& Conversion ) noexcept
    {
        m_X = static_cast<float>(Conversion[0]);
        m_Y = static_cast<float>(Conversion[1]);
    }

    //------------------------------------------------------------------------------
    // Conversion to array of doubles.
    //
    // Returns:
    //  An array with x and y as doubles.
    //
    constexpr fvec2::operator std::array<double,2> (void) const noexcept
    {
        return { static_cast<double>(m_X), static_cast<double>(m_Y) };
    }

    //------------------------------------------------------------------------------
    // Conversion to string representation.
    //
    // Returns:
    //  A string in "(x, y)" format.
    //
    inline fvec2::operator std::string(void) const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // Returns a string representation of the vector.
    //
    // Returns:
    //  A string in "(x, y)" format.
    //
    std::string fvec2::ToString(void) const noexcept
    {
        return std::format("({}, {})", m_X, m_Y);
    }

    //------------------------------------------------------------------------------
    // operator<<
    //
    // Overloads the stream output operator to print the vector in "(x, y)" format.
    //
    // Parameters:
    //    os   - The output stream.
    //    vec  - The vector to print.
    //
    // Returns:
    //    Reference to the output stream.
    //
    inline std::ostream& operator<< (std::ostream& os, const fvec2& vec) noexcept
    {
        return os << '(' << vec.m_X << ", " << vec.m_Y << ')';
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Returns a vector with all components set to 0.
    //
    constexpr fvec2 fvec2::fromZero(void) noexcept
    {
        return fvec2(0.f, 0.f);
    }

    //------------------------------------------------------------------------------
    // fromOne
    //------------------------------------------------------------------------------
    //
    // Returns a vector with all components set to 1.
    //
    constexpr fvec2 fvec2::fromOne(void) noexcept
    {
        return fvec2(1.f, 1.f);
    }

    //------------------------------------------------------------------------------
    // fromUnitX
    //------------------------------------------------------------------------------
    //
    // Returns the unit vector along the X-axis (1, 0).
    //
    constexpr fvec2 fvec2::fromUnitX(void) noexcept
    {
        return fvec2(1.f, 0.f);
    }

    //------------------------------------------------------------------------------
    // fromUnitY
    //------------------------------------------------------------------------------
    //
    // Returns the unit vector along the Y-axis (0, 1).
    //
    constexpr fvec2 fvec2::fromUnitY(void) noexcept
    {
        return fvec2(0.f, 1.f);
    }

    //------------------------------------------------------------------------------
    // fromUp
    //------------------------------------------------------------------------------
    //
    // Returns the up direction vector (0, 1).
    //
    constexpr fvec2 fvec2::fromUp(void) noexcept
    {
        return fvec2(0.f, 1.f);
    }

    //------------------------------------------------------------------------------
    // fromDown
    //------------------------------------------------------------------------------
    //
    // Returns the down direction vector (0, -1).
    //
    constexpr fvec2 fvec2::fromDown(void) noexcept
    {
        return fvec2(0.f, -1.f);
    }

    //------------------------------------------------------------------------------
    // fromLeft
    //------------------------------------------------------------------------------
    //
    // Returns the left direction vector (-1, 0).
    //
    constexpr fvec2 fvec2::fromLeft(void) noexcept
    {
        return fvec2(-1.f, 0.f);
    }

    //------------------------------------------------------------------------------
    // fromRight
    //------------------------------------------------------------------------------
    //
    // Returns the right direction vector (1, 0).
    //
    constexpr fvec2 fvec2::fromRight(void) noexcept
    {
        return fvec2(1.f, 0.f);
    }

    //------------------------------------------------------------------------------
    // fromRandomUnitVector
    //------------------------------------------------------------------------------
    //
    // Returns a random unit vector.
    //
    // Returns:
    //  A normalized random vector.
    //
    // Notes:
    //  Uses rand() for randomness; seed if needed.
    //
    inline fvec2 fvec2::fromRandomUnitVector(void) noexcept
    {
        float angle = static_cast<float>(rand()) / RAND_MAX * 2.f * 3.1415926535f;
        return fvec2(std::cos(angle), std::sin(angle));
    }

    //------------------------------------------------------------------------------
    // Static methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes the dot product of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  The dot product.
    //
    constexpr float fvec2::Dot(const fvec2& a, const fvec2& b) noexcept
    {
        return a.m_X * b.m_X + a.m_Y * b.m_Y;
    }

    //------------------------------------------------------------------------------
    // Min
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise minimum of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  A vector with the minimum components.
    //
    constexpr fvec2 fvec2::Min(const fvec2& a, const fvec2& b) noexcept
    {
        return fvec2(std::min(a.m_X, b.m_X), std::min(a.m_Y, b.m_Y));
    }

    //------------------------------------------------------------------------------
    // Max
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise maximum of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  A vector with the maximum components.
    //
    constexpr fvec2 fvec2::Max(const fvec2& a, const fvec2& b) noexcept
    {
        return fvec2(std::max(a.m_X, b.m_X), std::max(a.m_Y, b.m_Y));
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Performs linear interpolation between two vectors.
    //
    // Parameters:
    //  a - The starting vector.
    //  b - The ending vector.
    //  t - The interpolation factor (0 to 1).
    //
    // Returns:
    //  The interpolated vector.
    //
    constexpr fvec2 fvec2::Lerp(const fvec2& a, const fvec2& b, float t) noexcept
    {
        return a + (b - a) * t;
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean distance between two points.
    //
    // Parameters:
    //  a - The first point.
    //  b - The second point.
    //
    // Returns:
    //  The distance.
    //
    inline float fvec2::Distance(const fvec2& a, const fvec2& b) noexcept
    {
        return (a - b).Length();
    }

    //------------------------------------------------------------------------------
    // Cross
    //------------------------------------------------------------------------------
    //
    // Computes the 2D cross product (determinant) of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  The cross product scalar.
    //
    // Notes:
    //  Equivalent to a.x * b.y - a.y * b.x.
    //
    constexpr float fvec2::Cross(const fvec2& a, const fvec2& b) noexcept
    {
        return a.m_X * b.m_Y - a.m_Y * b.m_X;
    }

    //------------------------------------------------------------------------------
    // Static methods as members
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes the dot product with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  The dot product.
    //
    constexpr float fvec2::Dot(const fvec2& a) const noexcept
    {
        return Dot(*this, a);
    }

    //------------------------------------------------------------------------------
    // Min
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise minimum with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  A vector with the minimum components.
    //
    constexpr fvec2 fvec2::Min(const fvec2& a) const noexcept
    {
        return Min(*this, a);
    }

    //------------------------------------------------------------------------------
    // Max
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise maximum with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  A vector with the maximum components.
    //
    constexpr fvec2 fvec2::Max(const fvec2& a) const noexcept
    {
        return Max(*this, a);
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Performs linear interpolation to another vector.
    //
    // Parameters:
    //  a - The ending vector.
    //  t - The interpolation factor (0 to 1).
    //
    // Returns:
    //  The interpolated vector.
    //
    constexpr fvec2 fvec2::Lerp(const fvec2& a, float t) const noexcept
    {
        return Lerp(*this, a, t);
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean distance to another point.
    //
    // Parameters:
    //  a - The other point.
    //
    // Returns:
    //  The distance.
    //
    inline float fvec2::Distance(const fvec2& a) const noexcept
    {
        return Distance(*this, a);
    }

    //------------------------------------------------------------------------------
    // Cross
    //------------------------------------------------------------------------------
    //
    // Computes the 2D cross product with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  The cross product scalar.
    //
    constexpr float fvec2::Cross(const fvec2& a) const noexcept
    {
        return Cross(*this, a);
    }

    //------------------------------------------------------------------------------
    // Instance methods - Basic operations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Length
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean length (magnitude) of the vector.
    //
    // Returns:
    //  The length.
    //
    inline float fvec2::Length(void) const noexcept
    {
        return xmath::Sqrt(LengthSq());
    }

    //------------------------------------------------------------------------------
    // LengthSq
    //------------------------------------------------------------------------------
    //
    // Computes the squared length of the vector.
    //
    // Returns:
    //  The squared length.
    //
    // Notes:
    //  Faster than Length() for comparisons.
    //
    constexpr float fvec2::LengthSq(void) const noexcept
    {
        return m_X * m_X + m_Y * m_Y;
    }

    //------------------------------------------------------------------------------
    // NormalizeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy of the vector.
    //
    // Returns:
    //  The normalized vector.
    //
    // Notes:
    //  Returns zero vector if length is zero.
    //
    inline fvec2 fvec2::NormalizeCopy(void) const noexcept
    {
        float len = Length();
        return len > 0.f ? *this / len : fvec2(0.f);
    }

    //------------------------------------------------------------------------------
    // Normalize
    //------------------------------------------------------------------------------
    //
    // Normalizes the vector in-place.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Assumes non-zero length; use NormalizeSafe for zero-check.
    //
    inline fvec2& fvec2::Normalize(void) noexcept
    {
        float len = Length();
        assert(len > 0.f);
        *this /= len;
        return *this;
    }

    //------------------------------------------------------------------------------
    // NormalizeSafeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a safe normalized copy, handling infinite/NaN and zero length.
    //
    // Returns:
    //  The normalized vector or zero.
    //
    inline fvec2 fvec2::NormalizeSafeCopy(void) const noexcept
    {
        if (!isFinite()) return fvec2(0.f);
        float len = Length();
        return len > 0.f ? *this / len : fvec2(0.f);
    }

    //------------------------------------------------------------------------------
    // NormalizeSafe
    //------------------------------------------------------------------------------
    //
    // Safely normalizes the vector in-place, handling infinite/NaN and zero.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    inline fvec2& fvec2::NormalizeSafe(void) noexcept
    {
        if (!isFinite())
        {
            *this = fvec2(0.f);
            return *this;
        }
        float len = Length();
        if (len > 0.f) *this /= len;
        else *this = fvec2(0.f);
        return *this;
    }

    //------------------------------------------------------------------------------
    // LimitLengthCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with length limited to MaxLength.
    //
    // Parameters:
    //  MaxLength - The maximum length.
    //
    // Returns:
    //  The limited vector.
    //
    inline fvec2 fvec2::LimitLengthCopy(float MaxLength) const noexcept
    {
        float l = LengthSq();
        if (l <= (MaxLength * MaxLength)) return *this;
        return (*this) * (MaxLength * xmath::InvSqrt(l));
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all components are finite (not NaN or infinite).
    //
    // Returns:
    //  True if finite, false otherwise.
    //
    inline bool fvec2::isFinite(void) const noexcept
    {
        return std::isfinite(m_X) && std::isfinite(m_Y);
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if all components are within a given range.
    //
    // Parameters:
    //  min - The minimum value.
    //  max - The maximum value.
    //
    // Returns:
    //  True if all components are in [min, max].
    //
    constexpr bool fvec2::isInRange(float min, float max) const noexcept
    {
        return m_X >= min && m_X <= max && m_Y >= min && m_Y <= max;
    }

    //------------------------------------------------------------------------------
    // Equals
    //------------------------------------------------------------------------------
    //
    // Checks if this vector is approximately equal to another within a tolerance.
    //
    // Parameters:
    //  other - The vector to compare.
    //  tolerance - The epsilon tolerance.
    //
    // Returns:
    //  True if equal within tolerance, false otherwise.
    //
    inline bool fvec2::Equals(const fvec2& other, float tolerance) const noexcept
    {
        return xmath::Abs(m_X - other.m_X) <= tolerance && xmath::Abs(m_Y - other.m_Y) <= tolerance;
    }

    //------------------------------------------------------------------------------
    // Instance methods - Component-wise math
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // AbsCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with absolute value components.
    //
    // Returns:
    //  The absolute vector.
    //
    constexpr fvec2 fvec2::AbsCopy(void) const noexcept
    {
        return fvec2(xmath::Abs(m_X), xmath::Abs(m_Y));
    }

    //------------------------------------------------------------------------------
    // Abs
    //------------------------------------------------------------------------------
    //
    // Computes absolute values in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    constexpr fvec2& fvec2::Abs(void) noexcept
    {
        m_X = xmath::Abs(m_X);
        m_Y = xmath::Abs(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // OneOverCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with reciprocal components.
    //
    // Returns:
    //  The reciprocal vector.
    //
    // Notes:
    //  Assumes non-zero components; may produce inf/NaN.
    //
    constexpr fvec2 fvec2::OneOverCopy(void) const noexcept
    {
        return fvec2(1.f / m_X, 1.f / m_Y);
    }

    //------------------------------------------------------------------------------
    // OneOver
    //------------------------------------------------------------------------------
    //
    // Computes reciprocal in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Assumes non-zero components; may produce inf/NaN.
    //
    constexpr fvec2& fvec2::OneOver(void) noexcept
    {
        m_X = 1.f / m_X;
        m_Y = 1.f / m_Y;
        return *this;
    }

    //------------------------------------------------------------------------------
    // SqrtCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise square root of the vector.
    //
    // Returns:
    //  A vector with the square root of each component.
    //
    inline fvec2 fvec2::SqrtCopy(void) const noexcept
    {
        return fvec2(xmath::Sqrt(m_X), xmath::Sqrt(m_Y));
    }

    //------------------------------------------------------------------------------
    // Sqrt
    //------------------------------------------------------------------------------
    //
    // Computes the square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Sqrt(void) noexcept
    {
        m_X = xmath::Sqrt(m_X);
        m_Y = xmath::Sqrt(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // InvSqrtCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise inverse square root (1/sqrt(x)) of the vector.
    //
    // Returns:
    //  A vector with the inverse square root of each component.
    //
    inline fvec2 fvec2::InvSqrtCopy(void) const noexcept
    {
        return fvec2(xmath::InvSqrt(m_X), xmath::InvSqrt(m_Y));
    }

    //------------------------------------------------------------------------------
    // InvSqrt
    //------------------------------------------------------------------------------
    //
    // Computes the inverse square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::InvSqrt(void) noexcept
    {
        m_X = xmath::InvSqrt(m_X);
        m_Y = xmath::InvSqrt(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // SignCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise sign of the vector.
    //
    // Returns:
    //  A vector with -1 (negative), 0 (zero), or 1 (positive) per component.
    //
    inline fvec2 fvec2::SignCopy(void) const noexcept
    {
        return fvec2((m_X > 0.f) ? 1.f : (m_X < 0.f) ? -1.f : 0.f,
                     (m_Y > 0.f) ? 1.f : (m_Y < 0.f) ? -1.f : 0.f);
    }

    //------------------------------------------------------------------------------
    // Sign
    //------------------------------------------------------------------------------
    //
    // Computes the sign of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Sign(void) noexcept
    {
        m_X = (m_X > 0.f) ? 1.f : (m_X < 0.f) ? -1.f : 0.f;
        m_Y = (m_Y > 0.f) ? 1.f : (m_Y < 0.f) ? -1.f : 0.f;
        return *this;
    }

    //------------------------------------------------------------------------------
    // FloorCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise floor of the vector.
    //
    // Returns:
    //  A vector with the floor of each component.
    //
    inline fvec2 fvec2::FloorCopy(void) const noexcept
    {
        return fvec2(xmath::Floor(m_X), xmath::Floor(m_Y));
    }

    //------------------------------------------------------------------------------
    // Floor
    //------------------------------------------------------------------------------
    //
    // Computes the floor of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Floor(void) noexcept
    {
        m_X = xmath::Floor(m_X);
        m_Y = xmath::Floor(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // CeilCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise ceiling of the vector.
    //
    // Returns:
    //  A vector with the ceiling of each component.
    //
    inline fvec2 fvec2::CeilCopy(void) const noexcept
    {
        return fvec2(xmath::Ceil(m_X), xmath::Ceil(m_Y));
    }

    //------------------------------------------------------------------------------
    // Ceil
    //------------------------------------------------------------------------------
    //
    // Computes the ceiling of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Ceil(void) noexcept
    {
        m_X = xmath::Ceil(m_X);
        m_Y = xmath::Ceil(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // FractCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise fractional part of the vector.
    //
    // Returns:
    //  A vector with the fractional part of each component.
    //
    inline fvec2 fvec2::FractCopy(void) const noexcept
    {
        return fvec2(m_X - xmath::Floor(m_X), m_Y - xmath::Floor(m_Y));
    }

    //------------------------------------------------------------------------------
    // Fract
    //------------------------------------------------------------------------------
    //
    // Computes the fractional part of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Fract(void) noexcept
    {
        m_X -= xmath::Floor(m_X);
        m_Y -= xmath::Floor(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RoundCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise rounding of the vector to the nearest integer.
    //
    // Returns:
    //  A vector with the rounded value of each component.
    //
    inline fvec2 fvec2::RoundCopy(void) const noexcept
    {
        return fvec2(xmath::Round(m_X), xmath::Round(m_Y));
    }

    //------------------------------------------------------------------------------
    // Round
    //------------------------------------------------------------------------------
    //
    // Computes the rounding of the vector in-place to the nearest integer.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Round(void) noexcept
    {
        m_X = xmath::Round(m_X);
        m_Y = xmath::Round(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // TruncCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise truncation of the vector (towards zero).
    //
    // Returns:
    //  A vector with the truncated value of each component.
    //
    inline fvec2 fvec2::TruncCopy(void) const noexcept
    {
        return fvec2(xmath::Trunc(m_X), xmath::Trunc(m_Y));
    }

    //------------------------------------------------------------------------------
    // Trunc
    //------------------------------------------------------------------------------
    //
    // Computes the truncation of the vector in-place (towards zero).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Trunc(void) noexcept
    {
        m_X = xmath::Trunc(m_X);
        m_Y = xmath::Trunc(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ModCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise modulo of the vector by a scalar divisor.
    //
    // Parameters:
    //  divisor - The scalar divisor.
    //
    // Returns:
    //  A vector with the modulo of each component.
    //
    inline fvec2 fvec2::ModCopy(float divisor) const noexcept
    {
        return fvec2(xmath::FMod(m_X, divisor), xmath::FMod(m_Y, divisor));
    }

    //------------------------------------------------------------------------------
    // Mod
    //------------------------------------------------------------------------------
    //
    // Computes the modulo of the vector in-place by a scalar divisor.
    //
    // Parameters:
    //  divisor - The scalar divisor.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Mod(float divisor) noexcept
    {
        m_X = xmath::FMod(m_X, divisor);
        m_Y = xmath::FMod(m_Y, divisor);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClampCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise clamping of the vector between min and max.
    //
    // Parameters:
    //  min_val - The minimum value.
    //  max_val - The maximum value.
    //
    // Returns:
    //  A vector with each component clamped.
    //
    inline fvec2 fvec2::ClampCopy(float min_val, float max_val) const noexcept
    {
        return fvec2(xmath::Clamp(m_X, min_val, max_val), xmath::Clamp(m_Y, min_val, max_val));
    }

    //------------------------------------------------------------------------------
    // Clamp
    //------------------------------------------------------------------------------
    //
    // Computes the clamping of the vector in-place between min and max.
    //
    // Parameters:
    //  min_val - The minimum value.
    //  max_val - The maximum value.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Clamp(float min_val, float max_val) noexcept
    {
        m_X = xmath::Clamp(m_X, min_val, max_val);
        m_Y = xmath::Clamp(m_Y, min_val, max_val);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Step
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise step function.
    //
    // Parameters:
    //  edge - The edge value.
    //
    // Returns:
    //  A vector with 0 (if component < edge) or 1 (otherwise).
    //
    inline fvec2 fvec2::Step(float edge) const noexcept
    {
        return fvec2((m_X < edge) ? 0.f : 1.f, (m_Y < edge) ? 0.f : 1.f);
    }

    //------------------------------------------------------------------------------
    // SmoothStep
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise smoothstep interpolation.
    //
    // Parameters:
    //  edge0 - The lower edge.
    //  edge1 - The upper edge.
    //
    // Returns:
    //  A vector with smoothstep interpolation per component.
    //
    inline fvec2 fvec2::SmoothStep(float edge0, float edge1) const noexcept
    {
        float tx = xmath::Clamp((m_X - edge0) / (edge1 - edge0), 0.f, 1.f);
        float ty = xmath::Clamp((m_Y - edge0) / (edge1 - edge0), 0.f, 1.f);
        return fvec2(tx * tx * (3.f - 2.f * tx), ty * ty * (3.f - 2.f * ty));
    }

    //------------------------------------------------------------------------------
    // LogCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise natural logarithm of the vector.
    //
    // Returns:
    //  A vector with the natural log of each component.
    //
    inline fvec2 fvec2::LogCopy(void) const noexcept
    {
        return fvec2(xmath::Log(m_X), xmath::Log(m_Y));
    }

    //------------------------------------------------------------------------------
    // Log
    //------------------------------------------------------------------------------
    //
    // Computes the natural logarithm of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Log(void) noexcept
    {
        m_X = xmath::Log(m_X);
        m_Y = xmath::Log(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Log2Copy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise base-2 logarithm of the vector.
    //
    // Returns:
    //  A vector with the base-2 log of each component.
    //
    inline fvec2 fvec2::Log2Copy(void) const noexcept
    {
        return fvec2(xmath::Log2(m_X), xmath::Log2(m_Y));
    }

    //------------------------------------------------------------------------------
    // Log2
    //------------------------------------------------------------------------------
    //
    // Computes the base-2 logarithm of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Log2(void) noexcept
    {
        m_X = xmath::Log2(m_X);
        m_Y = xmath::Log2(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PowCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise power (base^exp) of the vector.
    //
    // Parameters:
    //  exp - The exponent (scalar).
    //
    // Returns:
    //  A vector with each component raised to the exponent.
    //
    inline fvec2 fvec2::PowCopy(float exp) const noexcept
    {
        return fvec2(xmath::Pow(m_X, exp), xmath::Pow(m_Y, exp));
    }

    //------------------------------------------------------------------------------
    // Pow
    //------------------------------------------------------------------------------
    //
    // Computes the power of the vector in-place.
    //
    // Parameters:
    //  exp - The exponent (scalar).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Pow(float exp) noexcept
    {
        m_X = xmath::Pow(m_X, exp);
        m_Y = xmath::Pow(m_Y, exp);
        return *this;
    }

    //------------------------------------------------------------------------------
    // SinCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise sine of the vector (in radians).
    //
    // Returns:
    //  A vector with the sine of each component.
    //
    inline fvec2 fvec2::SinCopy(void) const noexcept
    {
        return fvec2(xmath::Sin(radian{m_X}), xmath::Sin(radian{m_Y}));
    }

    //------------------------------------------------------------------------------
    // Sin
    //------------------------------------------------------------------------------
    //
    // Computes the sine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Sin(void) noexcept
    {
        m_X = xmath::Sin(radian{m_X});
        m_Y = xmath::Sin(radian{m_Y});
        return *this;
    }

    //------------------------------------------------------------------------------
    // CosCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise cosine of the vector (in radians).
    //
    // Returns:
    //  A vector with the cosine of each component.
    //
    inline fvec2 fvec2::CosCopy(void) const noexcept
    {
        return fvec2(xmath::Cos(radian{m_X}), xmath::Cos(radian{m_Y}));
    }

    //------------------------------------------------------------------------------
    // Cos
    //------------------------------------------------------------------------------
    //
    // Computes the cosine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Cos(void) noexcept
    {
        m_X = xmath::Cos(radian{m_X});
        m_Y = xmath::Cos(radian{m_Y});
        return *this;
    }

    //------------------------------------------------------------------------------
    // TanCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise tangent of the vector (in radians).
    //
    // Returns:
    //  A vector with the tangent of each component.
    //
    inline fvec2 fvec2::TanCopy(void) const noexcept
    {
        return fvec2(xmath::Tan(radian{m_X}), xmath::Tan(radian{m_Y}));
    }

    //------------------------------------------------------------------------------
    // Tan
    //------------------------------------------------------------------------------
    //
    // Computes the tangent of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Tan(void) noexcept
    {
        m_X = xmath::Tan(radian{m_X});
        m_Y = xmath::Tan(radian{m_Y});
        return *this;
    }

    //------------------------------------------------------------------------------
    // AsinCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arcsine of the vector.
    //
    // Returns:
    //  A vector with the arcsine of each component (in radians).
    //
    inline fvec2 fvec2::AsinCopy(void) const noexcept
    {
        return fvec2(xmath::Asin(m_X).m_Value, xmath::Asin(m_Y).m_Value);
    }

    //------------------------------------------------------------------------------
    // Asin
    //------------------------------------------------------------------------------
    //
    // Computes the arcsine of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Asin(void) noexcept
    {
        m_X = xmath::Asin(m_X).m_Value;
        m_Y = xmath::Asin(m_Y).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // AcosCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arccosine of the vector.
    //
    // Returns:
    //  A vector with the arccosine of each component (in radians).
    //
    inline fvec2 fvec2::AcosCopy(void) const noexcept
    {
        return fvec2(xmath::Acos(m_X).m_Value, xmath::Acos(m_Y).m_Value);
    }

    //------------------------------------------------------------------------------
    // Acos
    //------------------------------------------------------------------------------
    //
    // Computes the arccosine of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Acos(void) noexcept
    {
        m_X = xmath::Acos(m_X).m_Value;
        m_Y = xmath::Acos(m_Y).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // AtanCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arctangent of the vector.
    //
    // Returns:
    //  A vector with the arctangent of each component (in radians).
    //
    inline fvec2 fvec2::AtanCopy(void) const noexcept
    {
        return fvec2(xmath::Atan(m_X).m_Value, xmath::Atan(m_Y).m_Value);
    }

    //------------------------------------------------------------------------------
    // Atan
    //------------------------------------------------------------------------------
    //
    // Computes the arctangent of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Atan(void) noexcept
    {
        m_X = xmath::Atan(m_X).m_Value;
        m_Y = xmath::Atan(m_Y).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Atan2Copy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arctangent of y/x, considering the quadrant.
    //
    // Parameters:
    //  x - The vector containing x components.
    //
    // Returns:
    //  A vector with the arctangent of each y/x pair (in radians).
    //
    inline fvec2 fvec2::Atan2Copy(const fvec2& x) const noexcept
    {
        return fvec2(xmath::Atan2(m_Y, x.m_X).m_Value, xmath::Atan2(m_Y, x.m_Y).m_Value);
    }

    //------------------------------------------------------------------------------
    // Atan2
    //------------------------------------------------------------------------------
    //
    // Computes the arctangent of y/x in-place, considering the quadrant.
    //
    // Parameters:
    //  x - The vector containing x components.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Atan2(const fvec2& x) noexcept
    {
        return *this = Atan2Copy(x);
    }

    //------------------------------------------------------------------------------
    // SignedAngleBetween
    //------------------------------------------------------------------------------
    //
    // Computes the signed angle between this vector and another (in radians).
    //
    // Parameters:
    //  v - The other vector.
    //
    // Returns:
    //  The signed angle.
    //
    inline radian fvec2::SignedAngleBetween(const fvec2& v) const noexcept
    {
        return xmath::Atan2(Cross(v), Dot(v));
    }

    //------------------------------------------------------------------------------
    // ExpCopy
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise exponential (e^x) of the vector.
    //
    // Returns:
    //  A vector with the exponential of each component.
    //
    inline fvec2 fvec2::ExpCopy(void) const noexcept
    {
        return fvec2(xmath::Exp(m_X), xmath::Exp(m_Y));
    }

    //------------------------------------------------------------------------------
    // Exp
    //------------------------------------------------------------------------------
    //
    // Computes the exponential of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Exp(void) noexcept
    {
        m_X = xmath::Exp(m_X);
        m_Y = xmath::Exp(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Instance methods - Geometry
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Reflection
    //------------------------------------------------------------------------------
    //
    // Computes the reflection of the vector over a normal.
    //
    // Parameters:
    //  normal - The reflection normal (should be normalized).
    //
    // Returns:
    //  The reflected vector.
    //
    constexpr fvec2 fvec2::Reflection(const fvec2& normal) const noexcept
    {
        return *this - 2.f * Dot(normal) * normal;
    }

    //------------------------------------------------------------------------------
    // DistanceSquare
    //------------------------------------------------------------------------------
    //
    // Computes the squared distance to another point.
    //
    // Parameters:
    //  v - The other point.
    //
    // Returns:
    //  The squared distance.
    //
    constexpr float fvec2::DistanceSquare(const fvec2& v) const noexcept
    {
        return (*this - v).LengthSq();
    }

    //------------------------------------------------------------------------------
    // AngleBetween
    //------------------------------------------------------------------------------
    //
    // Computes the angle between this vector and another.
    //
    // Parameters:
    //  v - The other vector.
    //
    // Returns:
    //  The angle in radians.
    //
    // Notes:
    //  Vectors should be non-zero; returns 0 if either is zero.
    //
    inline radian fvec2::AngleBetween(const fvec2& v) const noexcept
    {
        float dot = Dot(v);
        float len = Length() * v.Length();
        if (len == 0.f) return radian{0.f};
        float cos = dot / len;
        cos = xmath::Clamp(cos, -1.f, 1.f);
        return xmath::Acos(cos);
    }

    //------------------------------------------------------------------------------
    // GridSnap
    //------------------------------------------------------------------------------
    //
    // Snaps the vector components to a grid.
    //
    // Parameters:
    //  gridX - Grid size for X.
    //  gridY - Grid size for Y.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    constexpr fvec2& fvec2::GridSnap(float gridX, float gridY) noexcept
    {
        m_X = xmath::Round(m_X / gridX) * gridX;
        m_Y = xmath::Round(m_Y / gridY) * gridY;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Perp
    //------------------------------------------------------------------------------
    //
    // Returns a perpendicular vector (rotated 90 degrees counterclockwise).
    //
    // Returns:
    //  The perpendicular vector (-y, x).
    //
    constexpr fvec2 fvec2::Perp(void) const noexcept
    {
        return fvec2(-m_Y, m_X);
    }

    //------------------------------------------------------------------------------
    // WhichSideOfLine
    //------------------------------------------------------------------------------
    //
    // Determines which side of a line the point is on.
    //
    // Parameters:
    //  V0 - First endpoint of the line.
    //  V1 - Second endpoint of the line.
    //
    // Returns:
    //  >0 if left, <0 if right, =0 if on the line.
    //
    // Notes:
    //  The return value divided by distance(V0,V1) is the min distance to the line.
    //
    constexpr float fvec2::WhichSideOfLine(const fvec2& V0, const fvec2& V1) const noexcept
    {
        return ((m_Y - V0.m_Y) * (V1.m_X - V0.m_X) - (m_X - V0.m_X) * (V1.m_Y - V0.m_Y));
    }

    //------------------------------------------------------------------------------
    // ClosestPointInLine
    //------------------------------------------------------------------------------
    //
    // Returns the closest point on the infinite line to this point.
    //
    // Parameters:
    //  V0 - First point on the line.
    //  V1 - Second point on the line.
    //
    // Returns:
    //  The closest point on the line.
    //
    // Notes:
    //  Not limited to the segment between V0 and V1.
    //
    inline fvec2 fvec2::ClosestPointInLine(const fvec2& V0, const fvec2& V1) const noexcept
    {
        // safety checks
        assert((V0.m_X != V1.m_X) && (V0.m_Y != V1.m_Y));

        float u = (m_X - V1.m_X) * (V1.m_X - V0.m_X) + (m_Y - V0.m_Y) * (V1.m_Y - V0.m_Y);
        u /= (V0 - V1).LengthSq();

        return V0.Lerp(V1, u);
    }

    //------------------------------------------------------------------------------
    // ClosestPointInLineSegment
    //------------------------------------------------------------------------------
    //
    // Returns the closest point on the line segment to this point.
    //
    // Parameters:
    //  V0 - First endpoint of the segment.
    //  V1 - Second endpoint of the segment.
    //
    // Returns:
    //  The closest point on the segment.
    //
    inline fvec2 fvec2::ClosestPointInLineSegment(const fvec2& V0, const fvec2& V1) const noexcept
    {
        // degenerate case
        if (V0 == V1)
        {
            return V0;
        }

        float u = (m_X - V1.m_X) * (V1.m_X - V0.m_X) + (m_Y - V0.m_Y) * (V1.m_Y - V0.m_Y);
        u /= (V0 - V1).LengthSq();

        // cap u to the range [0..1]
        u = xmath::Clamp(u, 0.0f, 1.0f);

        return V0.Lerp(V1, u);
    }

    //------------------------------------------------------------------------------
    // RotateCopy
    //------------------------------------------------------------------------------
    //
    // Computes the rotation of the vector by an angle (in radians).
    //
    // Parameters:
    //  angle - The rotation angle.
    //
    // Returns:
    //  A rotated vector.
    //
    inline fvec2 fvec2::RotateCopy(radian angle) const noexcept
    {
        float c = xmath::Cos(angle);
        float s = xmath::Sin(angle);
        return fvec2(m_X * c - m_Y * s, m_X * s + m_Y * c);
    }

    //------------------------------------------------------------------------------
    // Rotate
    //------------------------------------------------------------------------------
    //
    // Rotates the vector in-place by an angle (in radians).
    //
    // Parameters:
    //  angle - The rotation angle.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Rotate(radian angle) noexcept
    {
        float c = xmath::Cos(angle);
        float s = xmath::Sin(angle);
        float new_x = m_X * c - m_Y * s;
        m_Y = m_X * s + m_Y * c;
        m_X = new_x;
        return *this;
    }

    //------------------------------------------------------------------------------
    // ProjectCopy
    //------------------------------------------------------------------------------
    //
    // Computes the projection of this vector onto another.
    //
    // Parameters:
    //  onto - The vector to project onto.
    //
    // Returns:
    //  The projected vector.
    //
    inline fvec2 fvec2::ProjectCopy(const fvec2& onto) const noexcept
    {
        return onto * (Dot(onto) / onto.LengthSq());
    }

    //------------------------------------------------------------------------------
    // Project
    //------------------------------------------------------------------------------
    //
    // Projects the vector in-place onto another.
    //
    // Parameters:
    //  onto - The vector to project onto.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Project(const fvec2& onto) noexcept
    {
        *this = ProjectCopy(onto);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Swizzle methods for float (HLSL-style, return copy with swizzled components)
    //------------------------------------------------------------------------------

    constexpr float fvec2::x(void) const noexcept
    {
        return m_X;
    }

    constexpr float fvec2::y(void) const noexcept
    {
        return m_Y;
    }

    //------------------------------------------------------------------------------
    // Swizzle methods for fvec2 (HLSL-style, return copy with swizzled components)
    //------------------------------------------------------------------------------

    constexpr fvec2 fvec2::xx(void) const noexcept
    {
        return fvec2(m_X, m_X);
    }

    constexpr fvec2 fvec2::xy(void) const noexcept
    {
        return fvec2(m_X, m_Y);
    }

    constexpr fvec2 fvec2::yx(void) const noexcept
    {
        return fvec2(m_Y, m_X);
    }

    constexpr fvec2 fvec2::yy(void) const noexcept
    {
        return fvec2(m_Y, m_Y);
    }

    //------------------------------------------------------------------------------
    // Swizzle methods for fvec3 (HLSL-style, return copy with swizzled components)
    //------------------------------------------------------------------------------

    constexpr fvec3 fvec2::xxx(void) const noexcept
    {
        return fvec3(m_X, m_X, m_X);
    }

    constexpr fvec3 fvec2::xxy(void) const noexcept
    {
        return fvec3(m_X, m_X, m_Y);
    }

    constexpr fvec3 fvec2::xyx(void) const noexcept
    {
        return fvec3(m_X, m_Y, m_X);
    }

    constexpr fvec3 fvec2::xyy(void) const noexcept
    {
        return fvec3(m_X, m_Y, m_Y);
    }

    constexpr fvec3 fvec2::yxx(void) const noexcept
    {
        return fvec3(m_Y, m_X, m_X);
    }

    constexpr fvec3 fvec2::yxy(void) const noexcept
    {
        return fvec3(m_Y, m_X, m_Y);
    }

    constexpr fvec3 fvec2::yyx(void) const noexcept
    {
        return fvec3(m_Y, m_Y, m_X);
    }

    constexpr fvec3 fvec2::yyy(void) const noexcept
    {
        return fvec3(m_Y, m_Y, m_Y);
    }

    //------------------------------------------------------------------------------
    // Swizzle methods for fvec4 (HLSL-style, return copy with swizzled components)
    //------------------------------------------------------------------------------

    inline fvec4 fvec2::xxxx(void) const noexcept
    {
        return fvec4(m_X, m_X, m_X, m_X);
    }

    inline fvec4 fvec2::xxxy(void) const noexcept
    {
        return fvec4(m_X, m_X, m_X, m_Y);
    }

    inline fvec4 fvec2::xxyx(void) const noexcept
    {
        return fvec4(m_X, m_X, m_Y, m_X);
    }

    inline fvec4 fvec2::xxyy(void) const noexcept
    {
        return fvec4(m_X, m_X, m_Y, m_Y);
    }

    inline fvec4 fvec2::xyxx(void) const noexcept
    {
        return fvec4(m_X, m_Y, m_X, m_X);
    }

    inline fvec4 fvec2::xyxy(void) const noexcept
    {
        return fvec4(m_X, m_Y, m_X, m_Y);
    }

    inline fvec4 fvec2::xyyx(void) const noexcept
    {
        return fvec4(m_X, m_Y, m_Y, m_X);
    }

    inline fvec4 fvec2::xyyy(void) const noexcept
    {
        return fvec4(m_X, m_Y, m_Y, m_Y);
    }

    inline fvec4 fvec2::yxxx(void) const noexcept
    {
        return fvec4(m_Y, m_X, m_X, m_X);
    }

    inline fvec4 fvec2::yxxy(void) const noexcept
    {
        return fvec4(m_Y, m_X, m_X, m_Y);
    }

    inline fvec4 fvec2::yxyx(void) const noexcept
    {
        return fvec4(m_Y, m_X, m_Y, m_X);
    }

    inline fvec4 fvec2::yxyy(void) const noexcept
    {
        return fvec4(m_Y, m_X, m_Y, m_Y);
    }

    inline fvec4 fvec2::yyxx(void) const noexcept
    {
        return fvec4(m_Y, m_Y, m_X, m_X);
    }

    inline fvec4 fvec2::yyxy(void) const noexcept
    {
        return fvec4(m_Y, m_Y, m_X, m_Y);
    }

    inline fvec4 fvec2::yyyx(void) const noexcept
    {
        return fvec4(m_Y, m_Y, m_Y, m_X);
    }

    inline fvec4 fvec2::yyyy(void) const noexcept
    {
        return fvec4(m_Y, m_Y, m_Y, m_Y);
    }

    //------------------------------------------------------------------------------
    // Operator overloads
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator+
    //------------------------------------------------------------------------------
    //
    // Adds two vectors.
    //
    // Parameters:
    //  other - The vector to add.
    //
    // Returns:
    //  The sum vector.
    //
    constexpr fvec2 fvec2::operator+(const fvec2& other) const noexcept
    {
        return fvec2(m_X + other.m_X, m_Y + other.m_Y);
    }

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Subtracts two vectors.
    //
    // Parameters:
    //  other - The vector to subtract.
    //
    // Returns:
    //  The difference vector.
    //
    constexpr fvec2 fvec2::operator-(const fvec2& other) const noexcept
    {
        return fvec2(m_X - other.m_X, m_Y - other.m_Y);
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Multiplies by a scalar.
    //
    // Parameters:
    //  scalar - The scalar value.
    //
    // Returns:
    //  The scaled vector.
    //
    constexpr fvec2 fvec2::operator*(float scalar) const noexcept
    {
        return fvec2(m_X * scalar, m_Y * scalar);
    }

    //------------------------------------------------------------------------------
    // operator/
    //------------------------------------------------------------------------------
    //
    // Divides by a scalar.
    //
    // Parameters:
    //  scalar - The scalar value.
    //
    // Returns:
    //  The scaled vector.
    //
    // Notes:
    //  Assumes scalar != 0.
    //
    constexpr fvec2 fvec2::operator/(float scalar) const noexcept
    {
        float inv = 1.f / scalar;
        return *this * inv;
    }

    //------------------------------------------------------------------------------
    // operator+=
    //------------------------------------------------------------------------------
    //
    // Adds another vector in-place.
    //
    // Parameters:
    //  other - The vector to add.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    constexpr fvec2& fvec2::operator+=(const fvec2& other) noexcept
    {
        m_X += other.m_X;
        m_Y += other.m_Y;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator-=
    //------------------------------------------------------------------------------
    //
    // Subtracts another vector in-place.
    //
    // Parameters:
    //  other - The vector to subtract.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    constexpr fvec2& fvec2::operator-=(const fvec2& other) noexcept
    {
        m_X -= other.m_X;
        m_Y -= other.m_Y;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*=
    //------------------------------------------------------------------------------
    //
    // Multiplies by a scalar in-place.
    //
    // Parameters:
    //  scalar - The scalar value.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    constexpr fvec2& fvec2::operator*=(float scalar) noexcept
    {
        m_X *= scalar;
        m_Y *= scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator/=
    //------------------------------------------------------------------------------
    //
    // Divides by a scalar in-place.
    //
    // Parameters:
    //  scalar - The scalar value.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Assumes scalar != 0.
    //
    constexpr fvec2& fvec2::operator/=(float scalar) noexcept
    {
        float inv = 1.f / scalar;
        return *this *= inv;
    }

    //------------------------------------------------------------------------------
    // operator==
    //------------------------------------------------------------------------------
    //
    // Checks for equality with another vector.
    //
    // Parameters:
    //  other - The other vector.
    //
    // Returns:
    //  True if equal, false otherwise.
    //
    constexpr bool fvec2::operator==(const fvec2& other) const noexcept
    {
        return m_X == other.m_X && m_Y == other.m_Y;
    }

    //------------------------------------------------------------------------------
    // operator!=
    //------------------------------------------------------------------------------
    //
    // Checks for inequality with another vector.
    //
    // Parameters:
    //  other - The other vector.
    //
    // Returns:
    //  True if not equal, false otherwise.
    //
    constexpr bool fvec2::operator!=(const fvec2& other) const noexcept
    {
        return !(*this == other);
    }

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses a component by index (const).
    //
    // Parameters:
    //  index - The index (0 for x, 1 for y).
    //
    // Returns:
    //  The component value.
    //
    // Notes:
    //  Asserts index in [0,1].
    //
    constexpr float fvec2::operator[](std::int32_t index) const noexcept
    {
        assert(index >= 0 && index < 2);
        return m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses a component by index (mutable).
    //
    // Parameters:
    //  index - The index (0 for x, 1 for y).
    //
    // Returns:
    //  Reference to the component.
    //
    // Notes:
    //  Asserts index in [0,1].
    //
    constexpr float& fvec2::operator[](std::int32_t index) noexcept
    {
        assert(index >= 0 && index < 2);
        return m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // Friend operators
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Multiplies a scalar by a vector.
    //
    // Parameters:
    //  scalar - The scalar.
    //  v - The vector.
    //
    // Returns:
    //  The scaled vector.
    //
    constexpr fvec2 operator*(float scalar, const fvec2& v) noexcept
    {
        return v * scalar;
    }

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Negates a vector.
    //
    // Parameters:
    //  v - The vector.
    //
    // Returns:
    //  The negated vector.
    //
    constexpr fvec2 operator- (const fvec2& v) noexcept
    {
        return fvec2(-v.m_X, -v.m_Y);
    }
}

