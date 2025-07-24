#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
    #error "You must include xmath_vector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------
    constexpr fvec2::fvec2(float x, float y) noexcept
        : m_X(x), m_Y(y)
    {
    }

    //------------------------------------------------------------------------------

    constexpr fvec2::fvec2(float value) noexcept
        : m_X(value), m_Y(value)
    {
    }

    //------------------------------------------------------------------------------

    constexpr fvec2::fvec2(std::span<float> Span) noexcept
    {
        if (Span.size() >= 2)
        {
            m_X = Span[0];
            m_Y = Span[1];
        }
    }

    //------------------------------------------------------------------------------

    constexpr fvec2::fvec2(const std::array<double, 2>& Conversion) noexcept
    {
        m_X = static_cast<float>(Conversion[0]);
        m_Y = static_cast<float>(Conversion[1]);
    }

    //------------------------------------------------------------------------------

    constexpr fvec2::operator std::array<double, 2>() const noexcept
    {
        return { m_X, m_Y };
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
    // UnitX
    //------------------------------------------------------------------------------
    //
    // Returns the unit vector along the X-axis (1, 0).
    //
    constexpr fvec2 fvec2::fromUnitX(void) noexcept
    {
        return fvec2(1.f, 0.f);
    }

    //------------------------------------------------------------------------------
    // UnitY
    //------------------------------------------------------------------------------
    //
    // Returns the unit vector along the Y-axis (0, 1).
    //
    constexpr fvec2 fvec2::fromUnitY(void) noexcept
    {
        return fvec2(0.f, 1.f);
    }

    //------------------------------------------------------------------------------
    // Up
    //------------------------------------------------------------------------------
    //
    // Returns the up direction vector (0, 1).
    //
    constexpr fvec2 fvec2::fromUp(void) noexcept
    {
        return fvec2(0.f, 1.f);
    }

    //------------------------------------------------------------------------------
    // Down
    //------------------------------------------------------------------------------
    //
    // Returns the down direction vector (0, -1).
    //
    constexpr fvec2 fvec2::fromDown(void) noexcept
    {
        return fvec2(0.f, -1.f);
    }

    //------------------------------------------------------------------------------
    // Left
    //------------------------------------------------------------------------------
    //
    // Returns the left direction vector (-1, 0).
    //
    constexpr fvec2 fvec2::fromLeft(void) noexcept
    {
        return fvec2(-1.f, 0.f);
    }

    //------------------------------------------------------------------------------
    // Right
    //------------------------------------------------------------------------------
    //
    // Returns the right direction vector (1, 0).
    //
    constexpr fvec2 fvec2::fromRight(void) noexcept
    {
        return fvec2(1.f, 0.f);
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
    // Instance methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Sin
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise sine of the vector (in radians).
    //
    // Returns:
    //  A vector with the sine of each component.
    //
    inline fvec2 fvec2::SinCopy() const noexcept
    {
        return fvec2(xmath::Sin(xmath::radian{m_X}), xmath::Sin(xmath::radian{m_Y}));
    }

    //------------------------------------------------------------------------------
    // Cos
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise cosine of the vector (in radians).
    //
    // Returns:
    //  A vector with the cosine of each component.
    //
    inline fvec2 fvec2::CosCopy() const noexcept
    {
        return fvec2(xmath::Cos(xmath::radian{ m_X }), xmath::Cos(xmath::radian{ m_Y }));
    }

    //------------------------------------------------------------------------------
    // Tan
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise tangent of the vector (in radians).
    //
    // Returns:
    //  A vector with the tangent of each component.
    //
    inline fvec2 fvec2::TanCopy() const noexcept
    {
        return fvec2(xmath::Tan(xmath::radian{ m_X }), xmath::Tan(xmath::radian{ m_Y }));
    }

    //------------------------------------------------------------------------------
    // Asin
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arcsine of the vector.
    //
    // Returns:
    //  A vector with the arcsine of each component (in radians).
    //
    inline fvec2 fvec2::AsinCopy() const noexcept
    {
        return fvec2(xmath::Asin(m_X).m_Value, xmath::Asin(m_Y).m_Value);
    }

    //------------------------------------------------------------------------------
    // Acos
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arccosine of the vector.
    //
    // Returns:
    //  A vector with the arccosine of each component (in radians).
    //
    inline fvec2 fvec2::AcosCopy() const noexcept
    {
        return fvec2(xmath::Acos(m_X).m_Value, xmath::Acos(m_Y).m_Value);
    }

    //------------------------------------------------------------------------------
    // Atan
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise arctangent of the vector.
    //
    // Returns:
    //  A vector with the arctangent of each component (in radians).
    //
    inline fvec2 fvec2::AtanCopy() const noexcept
    {
        return fvec2(xmath::Atan(m_X).m_Value, xmath::Atan(m_Y).m_Value);
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
        return fvec2(xmath::Atan2(m_X, x.m_X).m_Value, xmath::Atan2(m_Y, x.m_Y).m_Value);
    }

    //------------------------------------------------------------------------------

    inline fvec2& fvec2::Atan2(const fvec2& x) noexcept
    {
        return *this = Atan2(x);
    }

    //------------------------------------------------------------------------------
    // Exp
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise exponential (e^x) of the vector.
    //
    // Returns:
    //  A vector with the exponential of each component.
    //
    inline fvec2 fvec2::ExpCopy() const noexcept
    {
        return fvec2(xmath::Exp(m_X), xmath::Exp(m_Y));
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
    inline fvec2 fvec2::LogCopy() const noexcept
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
    inline fvec2& fvec2::Log() noexcept
    {
        return *this = LogCopy();
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
    inline fvec2 fvec2::Log2Copy() const noexcept
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
    inline fvec2& fvec2::Log2() noexcept
    {
        return *this = Log2Copy();
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
        return *this = PowCopy(exp);
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
    inline fvec2 fvec2::SqrtCopy() const noexcept
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
    inline fvec2& fvec2::Sqrt() noexcept
    {
        return *this = SqrtCopy();
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
    inline fvec2 fvec2::InvSqrtCopy() const noexcept
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
    inline fvec2& fvec2::InvSqrt() noexcept
    {
        return *this = InvSqrtCopy();
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
    inline fvec2 fvec2::SignCopy() const noexcept
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
    inline fvec2& fvec2::Sign() noexcept
    {
        return *this = SignCopy();
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
    inline fvec2 fvec2::FloorCopy() const noexcept
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
    inline fvec2& fvec2::Floor() noexcept
    {
        return *this = FloorCopy();
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
    inline fvec2 fvec2::CeilCopy() const noexcept
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
    inline fvec2& fvec2::Ceil() noexcept
    {
        return *this = CeilCopy();
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
    inline fvec2 fvec2::FractCopy() const noexcept
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
    inline fvec2& fvec2::Fract() noexcept
    {
        return *this = FractCopy();
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
        return fvec2(xmath::FMod( m_X, divisor), xmath::FMod(m_Y, divisor));
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
        return *this = ModCopy(divisor);
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
        return *this = ClampCopy(min_val, max_val);
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
    // Sin
    //------------------------------------------------------------------------------
    //
    // Computes the sine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec2& fvec2::Sin() noexcept
    {
        m_X = xmath::Sin(xmath::radian{m_X});
        m_Y = xmath::Sin(xmath::radian{m_Y});
        return *this;
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
    inline fvec2& fvec2::Cos() noexcept
    {
        m_X = xmath::Cos(xmath::radian{ m_X});
        m_Y = xmath::Cos(xmath::radian{ m_Y});
        return *this;
    }

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
    // length
    //------------------------------------------------------------------------------
    //
    // Computes magnitude (Euclidean length) of vector.
    //
    // Returns:
    //  sqrt(m_X^2 + m_Y^2 + z^2) >=0.
    //
    // Notes:
    //  SIMD uses dot + sqrt; scalar fallback.
    //  Assumes finite; asserts validity.
    //
    inline fvec2 fvec2::LimitLengthCopy(float MaxLength) const noexcept
    {
        auto l = LengthSq();
        if (l <= (MaxLength * MaxLength)) return *this;
        return (*this) * (MaxLength * xmath::InvSqrt(l));
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
    //  Reference to self (chainable).
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
        return fvec2(std::abs(m_X), std::abs(m_Y));
    }

    //------------------------------------------------------------------------------
    // Abs
    //------------------------------------------------------------------------------
    //
    // Computes absolute values in-place.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    constexpr fvec2& fvec2::Abs(void) noexcept
    {
        m_X = xmath::Abs(m_X);
        m_Y = xmath::Abs(m_Y);
        return *this;
    }

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
    //  Reference to self (chainable).
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
    // Determines which side of a line a is point on.
    //
    // Note that the value returned divided by the distance from line_v1 to line_v2
    // is the minimum distance from the point to the line.  That may be useful 
    // for determining the distance relationship between points, without actually
    // having to calculate the distance.
    //
    // pt:		Point
    // line_v1: Endpoint of the line
    // line_v2: Endpoint of the line
    //
    // returns: > 0.0f if pt is to the left of the line (line_v1->line_v2)
    //				< 0.0f if pt is to the right of the line (line_v1->line_v2)
    //				= 0.0f if pt is on the line
    //
    //------------------------------------------------------------------------------
    constexpr float fvec2::WhichSideOfLine(const fvec2& V0, const fvec2& V1) const noexcept
    {
        return ((m_Y - V0.m_Y) * (V1.m_X - V0.m_X) - (m_X - V0.m_X) * (V1.m_Y - V0.m_Y));
    }

    //------------------------------------------------------------------------------
    // Returns the closest point on the 2D LINE defined by line_v1 and line_v2
    // to the point pt.  Note that output is NOT necissarily between line_v1
    // and line_v2.
    //
    // Reference: http://astronomy.swin.edu.au/~pbourke/geometry/pointline/
    //
    // out:		(output) Closest point on the line.
    // pt:		Point
    // line_v1:	First point on the line.
    // line_v2: Second point on the line.
    //------------------------------------------------------------------------------
    inline fvec2 fvec2::ClosestPointInLine(const fvec2& V0, const fvec2& V1) const noexcept
    {
        // safety checks
        assert((V0.m_X != V1.m_X) && (V0.m_Y != V1.m_Y));

        float u = (m_X - V1.m_X) * (V1.m_X - V0.m_X) + (m_Y - V0.m_Y) * (V1.m_Y - V0.m_Y);
        u /= (V0 - V1).LengthSq();

        return V0.Lerp( V1, u );
    }

    //------------------------------------------------------------------------------
    // Returns the closest point on the 2D LINESEGMENT defined by line_v1 and line_v2
    // to the point pt.  Note that output WILL BE between line_v1 and line_v2 (or 
    // equal to one of them).
    //
    // out:		(output) Closest point on the line segment.
    // pt:		Point
    // line_v1: Endpoint of the line segment.
    // line_v2: Endpoint of the line segment.
    //------------------------------------------------------------------------------
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
        u = xmath::Range(u, 0.0f, 1.0f);

        return V0.Lerp(V1, u );
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
        return *this = RotateCopy(angle);
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
        return *this = ProjectCopy(onto);
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
    inline fvec2 fvec2::RoundCopy() const noexcept
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
    inline fvec2& fvec2::Round() noexcept
    {
        m_X = xmath::Round(m_X); m_Y = xmath::Round(m_Y); return *this;
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
    inline fvec2 fvec2::TruncCopy() const noexcept
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
    inline fvec2& fvec2::Trunc() noexcept
    {
        return *this = TruncCopy();
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
    // Swizzle methods for float
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // x
    //------------------------------------------------------------------------------
    //
    // Returns the x component.
    //
    constexpr float fvec2::x(void) const noexcept { return m_X; }
    constexpr float fvec2::y(void) const noexcept { return m_Y; }

    //------------------------------------------------------------------------------
    // Swizzle methods for fvec2
    //------------------------------------------------------------------------------

    constexpr fvec2 fvec2::xx(void) const noexcept { return fvec2(m_X, m_X); }
    constexpr fvec2 fvec2::xy(void) const noexcept { return fvec2(m_X, m_Y); }
    constexpr fvec2 fvec2::yx(void) const noexcept { return fvec2(m_Y, m_X); }
    constexpr fvec2 fvec2::yy(void) const noexcept { return fvec2(m_Y, m_Y); }

    //------------------------------------------------------------------------------
    // Swizzle methods for fvec3
    //------------------------------------------------------------------------------

    constexpr fvec3 fvec2::xxx(void) const noexcept { return fvec3(m_X, m_X, m_X); }
    constexpr fvec3 fvec2::xxy(void) const noexcept { return fvec3(m_X, m_X, m_Y); }
    constexpr fvec3 fvec2::xyx(void) const noexcept { return fvec3(m_X, m_Y, m_X); }
    constexpr fvec3 fvec2::xyy(void) const noexcept { return fvec3(m_X, m_Y, m_Y); }
    constexpr fvec3 fvec2::yxx(void) const noexcept { return fvec3(m_Y, m_X, m_X); }
    constexpr fvec3 fvec2::yxy(void) const noexcept { return fvec3(m_Y, m_X, m_Y); }
    constexpr fvec3 fvec2::yyx(void) const noexcept { return fvec3(m_Y, m_Y, m_X); }
    constexpr fvec3 fvec2::yyy(void) const noexcept { return fvec3(m_Y, m_Y, m_Y); }

    //------------------------------------------------------------------------------
    // Swizzle methods for fvec4
    //------------------------------------------------------------------------------

    inline fvec4 fvec2::xxxx(void) const noexcept { return fvec4(m_X, m_X, m_X, m_X); }
    inline fvec4 fvec2::xxxy(void) const noexcept { return fvec4(m_X, m_X, m_X, m_Y); }
    inline fvec4 fvec2::xxyx(void) const noexcept { return fvec4(m_X, m_X, m_Y, m_X); }
    inline fvec4 fvec2::xxyy(void) const noexcept { return fvec4(m_X, m_X, m_Y, m_Y); }
    inline fvec4 fvec2::xyxx(void) const noexcept { return fvec4(m_X, m_Y, m_X, m_X); }
    inline fvec4 fvec2::xyxy(void) const noexcept { return fvec4(m_X, m_Y, m_X, m_Y); }
    inline fvec4 fvec2::xyyx(void) const noexcept { return fvec4(m_X, m_Y, m_Y, m_X); }
    inline fvec4 fvec2::xyyy(void) const noexcept { return fvec4(m_X, m_Y, m_Y, m_Y); }
    inline fvec4 fvec2::yxxx(void) const noexcept { return fvec4(m_Y, m_X, m_X, m_X); }
    inline fvec4 fvec2::yxxy(void) const noexcept { return fvec4(m_Y, m_X, m_X, m_Y); }
    inline fvec4 fvec2::yxyx(void) const noexcept { return fvec4(m_Y, m_X, m_Y, m_X); }
    inline fvec4 fvec2::yxyy(void) const noexcept { return fvec4(m_Y, m_X, m_Y, m_Y); }
    inline fvec4 fvec2::yyxx(void) const noexcept { return fvec4(m_Y, m_Y, m_X, m_X); }
    inline fvec4 fvec2::yyxy(void) const noexcept { return fvec4(m_Y, m_Y, m_X, m_Y); }
    inline fvec4 fvec2::yyyx(void) const noexcept { return fvec4(m_Y, m_Y, m_Y, m_X); }
    inline fvec4 fvec2::yyyy(void) const noexcept { return fvec4(m_Y, m_Y, m_Y, m_Y); }

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
    constexpr fvec2 fvec2::operator+ (const fvec2& other) const noexcept
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
    constexpr fvec2 fvec2::operator- (const fvec2& other) const noexcept
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
    constexpr fvec2 fvec2::operator* (float scalar) const noexcept
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
    constexpr fvec2 fvec2::operator/ (float scalar) const noexcept
    {
        return fvec2(m_X / scalar, m_Y / scalar);
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
    //  Reference to self (chainable).
    //
    constexpr fvec2& fvec2::operator+= (const fvec2& other) noexcept
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
    //  Reference to self (chainable).
    //
    constexpr fvec2& fvec2::operator-= (const fvec2& other) noexcept
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
    //  Reference to self (chainable).
    //
    constexpr fvec2& fvec2::operator*= (float scalar) noexcept
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
    //  Reference to self (chainable).
    //
    // Notes:
    //  Assumes scalar != 0.
    //
    constexpr fvec2& fvec2::operator/= (float scalar) noexcept
    {
        m_X /= scalar;
        m_Y /= scalar;
        return *this;
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
    constexpr bool fvec2::operator== (const fvec2& other) const noexcept
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
    constexpr bool fvec2::operator!= (const fvec2& other) const noexcept
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
    constexpr float fvec2::operator[] (std::int32_t index) const noexcept
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
    constexpr float& fvec2::operator[] (std::int32_t index) noexcept
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
    constexpr fvec2 operator* (float scalar, const fvec2& v) noexcept
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
