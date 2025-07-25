#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
#error "You must include xmath_vector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Constructs a vector with specified x, y, z, w components.
    //
    // Parameters:
    //  x - The x component.
    //  y - The y component.
    //  z - The z component.
    //  w - The w component.
    //
    inline fvec4::fvec4(float x, float y, float z, float w) noexcept
        : m_XYZW(_mm_set_ps(w, z, y, x))
    {
    }

    //------------------------------------------------------------------------------
    // Constructs a vector with all components set to the same value.
    //
    // Parameters:
    //  value - The value for all components.
    //
    inline fvec4::fvec4(float value) noexcept
        : m_XYZW(_mm_set_ps1(value))
    {
    }

    //------------------------------------------------------------------------------
    // Constructs a vector from a fvec3 and w component.
    //
    // Parameters:
    //  other - The fvec3 for x, y, z.
    //  w - The w component.
    //
    inline fvec4::fvec4(const fvec3& other, float w) noexcept
        : m_XYZW(_mm_insert_ps(other.m_XYZW, _mm_set_ss(w), 0x30)) // SSE4.1, assume available or alternative
    {
    }

    //------------------------------------------------------------------------------
    // Constructs a vector from x and a fvec3 for y, z, w.
    //
    // Parameters:
    //  x - The x component.
    //  other - The fvec3 for y, z, w.
    //
    // Notes:
    //  Bug: Original code sets ps with other.m_Z as w, but should set w as other.m_Z if intended, but according to class, it's x, y=other.m_X, z=other.m_Y, w=other.m_Z? Wait, the constructor is float x, const fvec3& other, so likely y=other.x, z=other.y, w=other.z. Fixed to member assignment.
    //
    inline fvec4::fvec4(float x, const fvec3& other) noexcept
        : m_X(x), m_Y(other.m_X), m_Z(other.m_Y), m_W(other.m_Z)
    {
    }

    //------------------------------------------------------------------------------
    // Constructs from a SIMD register.
    //
    // Parameters:
    //  reg - The __m128 register.
    //
    constexpr fvec4::fvec4(const floatx4& reg) noexcept
        : m_XYZW(reg)
    {
    }

    //------------------------------------------------------------------------------
    // Constructs from two fvec2 for xy and zw.
    //
    // Parameters:
    //  xy - The fvec2 for x, y.
    //  zw - The fvec2 for z, w.
    //
    inline fvec4::fvec4(const fvec2& xy, const fvec2& zw) noexcept
        : m_X(xy.m_X), m_Y(xy.m_Y), m_Z(zw.m_X), m_W(zw.m_Y)
    {
    }

    //------------------------------------------------------------------------------
    // Constructs from a span of floats (at least 4 elements).
    //
    // Parameters:
    //  Span - The span containing x, y, z, w.
    //
    inline fvec4::fvec4(std::span<float> Span) noexcept
    {
        assert(Span.size() >= 4);
        m_XYZW = _mm_loadu_ps(Span.data());
    }

    //------------------------------------------------------------------------------
    // Conversion to array of doubles.
    //
    // Returns:
    //  An array with x, y, z, w as doubles.
    //
    constexpr fvec4::operator std::array<double, 4>(void) const noexcept
    {
        return { static_cast<double>(m_X), static_cast<double>(m_Y), static_cast<double>(m_Z), static_cast<double>(m_W) };
    }

    //------------------------------------------------------------------------------
    // Conversion to string representation.
    //
    // Returns:
    //  A string in "(x, y, z, w)" format.
    //
    inline fvec4::operator std::string(void) const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // Returns a string representation of the vector.
    //
    // Returns:
    //  A string in "(x, y, z, w)" format.
    //
    std::string fvec4::ToString(void) const noexcept
    {
        return std::format("({}, {}, {}, {})", m_X, m_Y, m_Z, m_W);
    }

    //------------------------------------------------------------------------------
    // operator<<
    //
    // Overloads the stream output operator to print the vector in "(x, y, z, w)" format.
    //
    // Parameters:
    //    os   - The output stream.
    //    vec  - The vector to print.
    //
    // Returns:
    //    Reference to the output stream.
    //
    inline std::ostream& operator<< (std::ostream& os, const fvec4& vec) noexcept
    {
        return os << '(' << vec.m_X << ", " << vec.m_Y << ", " << vec.m_Z << ", " << vec.m_W << ')';
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Returns a vector with all components set to 0.
    //
    constexpr fvec4 fvec4::fromZero(void) noexcept
    {
        return fvec4(_mm_setzero_ps());
    }

    //------------------------------------------------------------------------------
    // Returns a vector with all components set to 1.
    //
    constexpr fvec4 fvec4::fromOne(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 1,1,1,1 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the unit vector along the X-axis (1, 0, 0, 0).
    //
    constexpr fvec4 fvec4::fromUnitX(void) noexcept
    {
        return fromRight();
    }

    //------------------------------------------------------------------------------
    // Returns the unit vector along the Y-axis (0, 1, 0, 0).
    //
    constexpr fvec4 fvec4::fromUnitY(void) noexcept
    {
        return fromUp();
    }

    //------------------------------------------------------------------------------
    // Returns the unit vector along the Z-axis (0, 0, 1, 0).
    //
    constexpr fvec4 fvec4::fromUnitZ(void) noexcept
    {
        return fromForward();
    }

    //------------------------------------------------------------------------------
    // Returns the unit vector along the W-axis (0, 0, 0, 1).
    //
    constexpr fvec4 fvec4::fromUnitW(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,0,0,1 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the up direction vector (0, 1, 0, 0).
    //
    constexpr fvec4 fvec4::fromUp(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,1,0,0 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the down direction vector (0, -1, 0, 0).
    //
    constexpr fvec4 fvec4::fromDown(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,-1,0,0 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the left direction vector (-1, 0, 0, 0).
    //
    constexpr fvec4 fvec4::fromLeft(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ -1, 0, 0, 0 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the right direction vector (1, 0, 0, 0).
    //
    constexpr fvec4 fvec4::fromRight(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 1,0,0,0 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the forward direction vector (0, 0, 1, 0).
    //
    constexpr fvec4 fvec4::fromForward(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,0,1,0 }} };
    }

    //------------------------------------------------------------------------------
    // Returns the back direction vector (0, 0, -1, 0).
    //
    constexpr fvec4 fvec4::fromBack(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,0,-1,0 }} };
    }

    //------------------------------------------------------------------------------
    // Returns a random unit vector.
    //
    // Returns:
    //  A normalized random vector.
    //
    // Notes:
    //  Uses random angles; scalar operation.
    //
    inline fvec4 fvec4::fromRandomUnitVector(void) noexcept
    {
        const radian theta{ static_cast<float>(rand()) / RAND_MAX * 2.f * xmath::pi_v.m_Value };
        const radian phi  { static_cast<float>(rand()) / RAND_MAX * xmath::pi_v.m_Value };
        const float sin_phi   = xmath::Sin(phi);
        const float cos_theta = xmath::Cos(theta);
        const float sin_theta = xmath::Sin(theta);
        return fvec4(sin_phi * cos_theta, sin_phi * sin_theta, xmath::Cos(phi), 1.f).Normalize();
    }

    //------------------------------------------------------------------------------
    // Homogenize
    //------------------------------------------------------------------------------
    //
    // Makes the instance class homogeneous
    //
    // Params:
    //  none
    //
    // Returns:
    //  instance of the class 
    //
    // Notes:
    //  
    //
    inline
    fvec4& fvec4::Homogenize(void) noexcept
    {
        return *this /= m_W;
    }

    //------------------------------------------------------------------------------
    // HomogeneousCopy
    //------------------------------------------------------------------------------
    //
    // Makes a copy of the vector homogeneous
    //
    // Params:
    //  none
    //
    // Returns:
    //  instance a homogeneous vector
    //
    // Notes:
    //  
    //
    inline
    fvec4 fvec4::HomogeneousCopy(void) const noexcept
    {
        return *this / m_W;
    }

    //------------------------------------------------------------------------------
    // Static methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Computes the dot product of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  The dot product.
    //
    inline float fvec4::Dot(const fvec4& a, const fvec4& b) noexcept
    {
        __m128 mul = _mm_mul_ps(a.m_XYZW, b.m_XYZW);
        mul = _mm_hadd_ps(mul, mul);
        mul = _mm_hadd_ps(mul, mul);
        return _mm_cvtss_f32(mul);
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise minimum of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  A vector with the minimum components.
    //
    inline fvec4 fvec4::Min(const fvec4& a, const fvec4& b) noexcept
    {
        return fvec4(_mm_min_ps(a.m_XYZW, b.m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise maximum of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  A vector with the maximum components.
    //
    inline fvec4 fvec4::Max(const fvec4& a, const fvec4& b) noexcept
    {
        return fvec4(_mm_max_ps(a.m_XYZW, b.m_XYZW));
    }

    //------------------------------------------------------------------------------
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
    inline fvec4 fvec4::Lerp(const fvec4& a, const fvec4& b, float t) noexcept
    {
        __m128 tt = _mm_set1_ps(t);
        return fvec4(_mm_add_ps(a.m_XYZW, _mm_mul_ps(tt, _mm_sub_ps(b.m_XYZW, a.m_XYZW))));
    }

    //------------------------------------------------------------------------------
    // Computes the Euclidean distance between two points.
    //
    // Parameters:
    //  a - The first point.
    //  b - The second point.
    //
    // Returns:
    //  The distance.
    //
    inline float fvec4::Distance(const fvec4& a, const fvec4& b) noexcept
    {
        return (a - b).Length();
    }

    //------------------------------------------------------------------------------
    // Static methods as members
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Computes the dot product with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  The dot product.
    //
    inline float fvec4::Dot(const fvec4& a) const noexcept
    {
        return Dot(*this, a);
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise minimum with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  A vector with the minimum components.
    //
    inline fvec4 fvec4::Min(const fvec4& a) const noexcept
    {
        return Min(*this, a);
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise maximum with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  A vector with the maximum components.
    //
    inline fvec4 fvec4::Max(const fvec4& a) const noexcept
    {
        return Max(*this, a);
    }

    //------------------------------------------------------------------------------
    // Performs linear interpolation to another vector.
    //
    // Parameters:
    //  a - The ending vector.
    //  t - The interpolation factor (0 to 1).
    //
    // Returns:
    //  The interpolated vector.
    //
    inline fvec4 fvec4::Lerp(const fvec4& a, float t) const noexcept
    {
        return Lerp(*this, a, t);
    }

    //------------------------------------------------------------------------------
    // Computes the Euclidean distance to another point.
    //
    // Parameters:
    //  a - The other point.
    //
    // Returns:
    //  The distance.
    //
    inline float fvec4::Distance(const fvec4& a) const noexcept
    {
        return Distance(*this, a);
    }

    //------------------------------------------------------------------------------
    // Instance methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Computes the Euclidean length of the vector.
    //
    // Returns:
    //  The length.
    //
    inline float fvec4::Length(void) const noexcept
    {
        return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(m_XYZW, m_XYZW, 0xF1)));
    }

    //------------------------------------------------------------------------------
    // Computes the squared length of the vector.
    //
    // Returns:
    //  The squared length.
    //
    // Notes:
    //  Faster than Length() for comparisons.
    //
    inline float fvec4::LengthSq(void) const noexcept
    {
        return _mm_cvtss_f32(_mm_dp_ps(m_XYZW, m_XYZW, 0xF1));
    }

    //------------------------------------------------------------------------------
    // Returns a copy with length limited to MaxLength.
    //
    // Parameters:
    //  MaxLength - The maximum length.
    //
    // Returns:
    //  The limited vector.
    //
    inline fvec4 fvec4::LimitLengthCopy(float MaxLength) const noexcept
    {
        float l2 = LengthSq();
        if (l2 > MaxLength * MaxLength)
        {
            float invLen = xmath::InvSqrt(l2);
            return *this * (MaxLength * invLen);
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Returns a normalized copy of the vector.
    //
    // Returns:
    //  The normalized vector.
    //
    // Notes:
    //  Returns zero vector if length is zero.
    //
    inline fvec4 fvec4::NormalizeCopy(void) const noexcept
    {
        float len = Length();
        if (len > 0.f) return *this / len;
        return fromZero();
    }

    //------------------------------------------------------------------------------
    // Normalizes the vector in-place.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Assumes non-zero length; use NormalizeSafe for zero-check.
    //
    inline fvec4& fvec4::Normalize(void) noexcept
    {
        float len = Length();
        assert(len > 0.f);
        *this /= len;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Returns a safe normalized copy, handling infinite/NaN and zero length.
    //
    // Returns:
    //  The normalized vector or zero.
    //
    inline fvec4 fvec4::NormalizeSafeCopy(void) const noexcept
    {
        if (!isFinite()) return fromZero();
        float len = Length();
        if (len > 0.f) return *this / len;
        return fromZero();
    }

    //------------------------------------------------------------------------------
    // Safely normalizes the vector in-place, handling infinite/NaN and zero.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    inline fvec4& fvec4::NormalizeSafe(void) noexcept
    {
        if (!isFinite())
        {
            *this = fromZero();
            return *this;
        }
        float len = Length();
        if (len > 0.f) *this /= len;
        else *this = fromZero();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Checks if all components are finite (not NaN or infinite).
    //
    // Returns:
    //  True if finite, false otherwise.
    //
    inline bool fvec4::isFinite(void) const noexcept
    {
        return std::isfinite(m_X) && std::isfinite(m_Y) && std::isfinite(m_Z) && std::isfinite(m_W);
    }

    //------------------------------------------------------------------------------
    // Checks if all components are within a given range.
    //
    // Parameters:
    //  min - The minimum value.
    //  max - The maximum value.
    //
    // Returns:
    //  True if all components are in [min, max].
    //
    inline bool fvec4::isInRange(float min, float max) const noexcept
    {
        return m_X >= min && m_X <= max && m_Y >= min && m_Y <= max && m_Z >= min && m_Z <= max && m_W >= min && m_W <= max;
    }

    //------------------------------------------------------------------------------
    // Returns a copy with reciprocal components.
    //
    // Returns:
    //  The reciprocal vector.
    //
    // Notes:
    //  Assumes non-zero components; may produce inf/NaN.
    //
    inline fvec4 fvec4::OneOverCopy(void) const noexcept
    {
        return fvec4(_mm_div_ps(_mm_set1_ps(1.f), m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes reciprocal in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Assumes non-zero components; may produce inf/NaN.
    //
    inline fvec4& fvec4::OneOver(void) noexcept
    {
        m_XYZW = _mm_div_ps(_mm_set1_ps(1.f), m_XYZW);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Returns a copy with absolute value components.
    //
    // Returns:
    //  The absolute vector.
    //
    inline fvec4 fvec4::AbsCopy(void) const noexcept
    {
        return fvec4(_mm_andnot_ps(_mm_set1_ps(-0.f), m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes absolute values in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Abs(void) noexcept
    {
        m_XYZW = _mm_andnot_ps(_mm_set1_ps(-0.f), m_XYZW);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the reflection of the vector over a normal.
    //
    // Parameters:
    //  normal - The reflection normal (should be normalized).
    //
    // Returns:
    //  The reflected vector.
    //
    inline fvec4 fvec4::Reflection(const fvec4& normal) const noexcept
    {
        return *this - 2.f * Dot(normal) * normal;
    }

    //------------------------------------------------------------------------------
    // Computes the squared distance to another point.
    //
    // Parameters:
    //  v - The other point.
    //
    // Returns:
    //  The squared distance.
    //
    inline float fvec4::DistanceSquare(const fvec4& v) const noexcept
    {
        return (*this - v).LengthSq();
    }

    //------------------------------------------------------------------------------
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
    inline radian fvec4::AngleBetween(const fvec4& v) const noexcept
    {
        float len = Length() * v.Length();
        if (len == 0.f) return radian(0.f);
        float cos = Dot(v) / len;
        cos = xmath::Clamp(cos, -1.f, 1.f);
        return xmath::Acos(cos);
    }

    //------------------------------------------------------------------------------
    // Snaps the vector components to a grid.
    //
    // Parameters:
    //  gridX - Grid size for X.
    //  gridY - Grid size for Y.
    //  gridZ - Grid size for Z.
    //  gridW - Grid size for W.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::GridSnap(float gridX, float gridY, float gridZ, float gridW) noexcept
    {
        m_X = xmath::Round(m_X / gridX) * gridX;
        m_Y = xmath::Round(m_Y / gridY) * gridY;
        m_Z = xmath::Round(m_Z / gridZ) * gridZ;
        m_W = xmath::Round(m_W / gridW) * gridW;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise square root of the vector.
    //
    // Returns:
    //  A vector with the square root of each component.
    //
    inline fvec4 fvec4::SqrtCopy(void) const noexcept
    {
        return fvec4(_mm_sqrt_ps(m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes the square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Sqrt(void) noexcept
    {
        m_XYZW = _mm_sqrt_ps(m_XYZW);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise inverse square root (1/sqrt(x)) of the vector.
    // High-precision version using Newton refinement.
    //
    // Returns:
    //  A vector with the inverse square root of each component.
    //
    // Notes:
    //  Suitable for high-accuracy requirements; slightly slower than InvSqrtCopyFast.
    //
    inline fvec4 fvec4::InvSqrtCopy(void) const noexcept
    {
        __m128 approx = _mm_rsqrt_ps(m_XYZW);
        __m128 half = _mm_set1_ps(0.5f);
        __m128 three_half = _mm_set1_ps(1.5f);
        __m128 refined = _mm_mul_ps(approx, _mm_sub_ps(three_half, _mm_mul_ps(half, _mm_mul_ps(m_XYZW, _mm_mul_ps(approx, approx)))));
        return fvec4(refined);
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise inverse square root (1/sqrt(x)) of the vector.
    // Fast version using raw SSE approximation.
    //
    // Returns:
    //  A vector with the approximate inverse square root of each component.
    //
    // Notes:
    //  Faster but less precise (~12 bits) than InvSqrtCopy; suitable for performance-critical code.
    //
    inline fvec4 fvec4::InvSqrtFastCopy(void) const noexcept
    {
        return fvec4(_mm_rsqrt_ps(m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes the inverse square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::InvSqrt(void) noexcept
    {
        return *this = InvSqrtCopy();
    }

    //------------------------------------------------------------------------------
    // Computes the inverse square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
//
    // Notes:
    //  Faster but less precise (~12 bits) than InvSqrtCopy; suitable for performance-critical code.
    //
    inline fvec4& fvec4::InvSqrtFast(void) noexcept
    {
        return *this = InvSqrtFastCopy();
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise sign of the vector.
    //
    // Returns:
    //  A vector with -1 (negative), 0 (zero), or 1 (positive) per component.
    //
    inline fvec4 fvec4::SignCopy(void) const noexcept
    {
        __m128 zero = _mm_setzero_ps();
        __m128 one = _mm_set1_ps(1.f);
        __m128 neg_one = _mm_set1_ps(-1.f);
        __m128 gt_zero = _mm_cmpgt_ps(m_XYZW, zero);
        __m128 lt_zero = _mm_cmplt_ps(m_XYZW, zero);
        __m128 sign_pos = _mm_and_ps(gt_zero, one);
        __m128 sign_neg = _mm_and_ps(lt_zero, neg_one);
        return fvec4(_mm_or_ps(sign_pos, sign_neg));
    }

    //------------------------------------------------------------------------------
    // Computes the sign of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Sign(void) noexcept
    {
        *this = SignCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise floor of the vector.
    //
    // Returns:
    //  A vector with the floor of each component.
    //
    inline fvec4 fvec4::FloorCopy(void) const noexcept
    {
        return fvec4(_mm_floor_ps(m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes the floor of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Floor(void) noexcept
    {
        m_XYZW = _mm_floor_ps(m_XYZW);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise ceiling of the vector.
    //
    // Returns:
    //  A vector with the ceiling of each component.
    //
    inline fvec4 fvec4::CeilCopy(void) const noexcept
    {
        return fvec4(_mm_ceil_ps(m_XYZW));
    }

    //------------------------------------------------------------------------------
    // Computes the ceiling of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Ceil(void) noexcept
    {
        m_XYZW = _mm_ceil_ps(m_XYZW);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise fractional part of the vector.
    //
    // Returns:
    //  A vector with the fractional part of each component.
    //
    inline fvec4 fvec4::FractCopy(void) const noexcept
    {
        return *this - TruncCopy();
    }

    //------------------------------------------------------------------------------
    // Computes the fractional part of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Fract(void) noexcept
    {
        *this -= TruncCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise rounding of the vector to the nearest integer.
    //
    // Returns:
    //  A vector with the rounded value of each component.
    //
    inline fvec4 fvec4::RoundCopy(void) const noexcept
    {
        return fvec4(_mm_round_ps(m_XYZW, _MM_FROUND_TO_NEAREST_INT));
    }

    //------------------------------------------------------------------------------
    // Computes the rounding of the vector in-place to the nearest integer.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Round(void) noexcept
    {
        m_XYZW = _mm_round_ps(m_XYZW, _MM_FROUND_TO_NEAREST_INT);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise truncation of the vector (towards zero).
    //
    // Returns:
    //  A vector with the truncated value of each component.
    //
    inline fvec4 fvec4::TruncCopy(void) const noexcept
    {
        return fvec4(_mm_round_ps(m_XYZW, _MM_FROUND_TO_ZERO));
    }

    //------------------------------------------------------------------------------
    // Computes the truncation of the vector in-place (towards zero).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Trunc(void) noexcept
    {
        m_XYZW = _mm_round_ps(m_XYZW, _MM_FROUND_TO_ZERO);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise modulo of the vector by a scalar divisor.
    //
    // Parameters:
    //  divisor - The scalar divisor.
    //
    // Returns:
    //  A vector with the modulo of each component.
    //
    inline fvec4 fvec4::ModCopy(float divisor) const noexcept
    {
        return fvec4(xmath::FMod(m_X, divisor), xmath::FMod(m_Y, divisor), xmath::FMod(m_Z, divisor), xmath::FMod(m_W, divisor));
    }

    //------------------------------------------------------------------------------
    // Computes the modulo of the vector in-place by a scalar divisor.
    //
    // Parameters:
    //  divisor - The scalar divisor.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Mod(float divisor) noexcept
    {
        m_X = xmath::FMod(m_X, divisor);
        m_Y = xmath::FMod(m_Y, divisor);
        m_Z = xmath::FMod(m_Z, divisor);
        m_W = xmath::FMod(m_W, divisor);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise clamping of the vector between min and max.
    //
    // Parameters:
    //  min_val - The minimum value.
    //  max_val - The maximum value.
    //
    // Returns:
    //  A vector with each component clamped.
    //
    inline fvec4 fvec4::ClampCopy(float min_val, float max_val) const noexcept
    {
        __m128 minv = _mm_set1_ps(min_val);
        __m128 maxv = _mm_set1_ps(max_val);
        return fvec4(_mm_max_ps(minv, _mm_min_ps(m_XYZW, maxv)));
    }

    //------------------------------------------------------------------------------
    // Computes the clamping of the vector in-place between min and max.
    //
    // Parameters:
    //  min_val - The minimum value.
    //  max_val - The maximum value.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Clamp(float min_val, float max_val) noexcept
    {
        *this = ClampCopy(min_val, max_val);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise clamping of the vector between min and max vectors.
    //
    // Parameters:
    //  min - The min vector.
    //  max - The max vector.
    //
    // Returns:
    //  A vector with each component clamped.
    //
    inline fvec4 fvec4::ClampCopy(const fvec4& min, const fvec4& max) const noexcept
    {
        return fvec4(_mm_max_ps(min.m_XYZW, _mm_min_ps(m_XYZW, max.m_XYZW)));
    }

    //------------------------------------------------------------------------------
    // Computes the clamping of the vector in-place between min and max vectors.
    //
    // Parameters:
    //  min - The min vector.
    //  max - The max vector.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Clamp(const fvec4& min, const fvec4& max) noexcept
    {
        m_XYZW = _mm_max_ps(min.m_XYZW, _mm_min_ps(m_XYZW, max.m_XYZW));
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise step function.
    //
    // Parameters:
    //  edge - The edge value.
    //
    // Returns:
    //  A vector with 0 (if component < edge) or 1 (otherwise).
    //
    inline fvec4 fvec4::Step(float edge) const noexcept
    {
        __m128 edgev = _mm_set1_ps(edge);
        __m128 cmp = _mm_cmpge_ps(m_XYZW, edgev);
        return fvec4(_mm_and_ps(cmp, _mm_set1_ps(1.f)));
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise smoothstep interpolation.
    //
    // Parameters:
    //  edge0 - The lower edge (scalar).
    //  edge1 - The upper edge (scalar).
    //
    // Returns:
    //  A vector with smoothstep interpolation per component.
    //
    // Notes:
    //  Clamps t to [0,1] before applying 3t^2 - 2t^3.
    //
    inline fvec4 fvec4::SmoothStep(float edge0, float edge1) const noexcept
    {
        __m128 edge0_vec = _mm_set1_ps(edge0);
        __m128 edge1_vec = _mm_set1_ps(edge1);
        __m128 t = _mm_div_ps(_mm_sub_ps(m_XYZW, edge0_vec), _mm_sub_ps(edge1_vec, edge0_vec));
        t = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(t, _mm_set1_ps(1.f))); // Clamp to [0,1]
        __m128 t2 = _mm_mul_ps(t, t);
        __m128 t3 = _mm_mul_ps(t2, t);
        __m128 result = _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(3.f), t2), _mm_mul_ps(_mm_set1_ps(2.f), t3));
        return fvec4(result);
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise natural logarithm of the vector.
    //
    // Returns:
    //  A vector with the natural log of each component.
    //
    inline fvec4 fvec4::LogCopy(void) const noexcept
    {
        return fvec4(xmath::Log(m_X), xmath::Log(m_Y), xmath::Log(m_Z), xmath::Log(m_W));
    }

    //------------------------------------------------------------------------------
    // Computes the natural logarithm of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Log(void) noexcept
    {
        m_X = xmath::Log(m_X);
        m_Y = xmath::Log(m_Y);
        m_Z = xmath::Log(m_Z);
        m_W = xmath::Log(m_W);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise base-2 logarithm of the vector.
    //
    // Returns:
    //  A vector with the base-2 log of each component.
    //
    inline fvec4 fvec4::Log2Copy(void) const noexcept
    {
        return fvec4(xmath::Log2(m_X), xmath::Log2(m_Y), xmath::Log2(m_Z), xmath::Log2(m_W));
    }

    //------------------------------------------------------------------------------
    // Computes the base-2 logarithm of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Log2(void) noexcept
    {
        m_X = xmath::Log2(m_X);
        m_Y = xmath::Log2(m_Y);
        m_Z = xmath::Log2(m_Z);
        m_W = xmath::Log2(m_W);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise power (base^exp) of the vector.
    //
    // Parameters:
    //  exp - The exponent (scalar).
    //
    // Returns:
    //  A vector with each component raised to the exponent.
    //
    inline fvec4 fvec4::PowCopy(float exp) const noexcept
    {
        return fvec4(xmath::Pow(m_X, exp), xmath::Pow(m_Y, exp), xmath::Pow(m_Z, exp), xmath::Pow(m_W, exp));
    }

    //------------------------------------------------------------------------------
    // Computes the power of the vector in-place.
    //
    // Parameters:
    //  exp - The exponent (scalar).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Pow(float exp) noexcept
    {
        m_X = xmath::Pow(m_X, exp);
        m_Y = xmath::Pow(m_Y, exp);
        m_Z = xmath::Pow(m_Z, exp);
        m_W = xmath::Pow(m_W, exp);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise sine of the vector (in radians).
    //
    // Returns:
    //  A vector with the sine of each component.
    //
    inline fvec4 fvec4::SinCopy(void) const noexcept
    {
        return fvec4(xmath::Sin(radian{ m_X }), xmath::Sin(radian{ m_Y }), xmath::Sin(radian{ m_Z }), xmath::Sin(radian{ m_W }));
    }

    //------------------------------------------------------------------------------
    // Computes the sine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Sin(void) noexcept
    {
        m_X = xmath::Sin(radian{ m_X });
        m_Y = xmath::Sin(radian{ m_Y });
        m_Z = xmath::Sin(radian{ m_Z });
        m_W = xmath::Sin(radian{ m_W });
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise cosine of the vector (in radians).
    //
    // Returns:
    //  A vector with the cosine of each component.
    //
    inline fvec4 fvec4::CosCopy(void) const noexcept
    {
        return fvec4(xmath::Cos(radian{ m_X }), xmath::Cos(radian{ m_Y }), xmath::Cos(radian{ m_Z }), xmath::Cos(radian{ m_W }));
    }

    //------------------------------------------------------------------------------
    // Computes the cosine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Cos(void) noexcept
    {
        m_X = xmath::Cos(radian{ m_X });
        m_Y = xmath::Cos(radian{ m_Y });
        m_Z = xmath::Cos(radian{ m_Z });
        m_W = xmath::Cos(radian{ m_W });
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise tangent of the vector (in radians).
    //
    // Returns:
    //  A vector with the tangent of each component.
    //
    inline fvec4 fvec4::TanCopy(void) const noexcept
    {
        return fvec4(xmath::Tan(radian{ m_X }), xmath::Tan(radian{ m_Y }), xmath::Tan(radian{ m_Z }), xmath::Tan(radian{ m_W }));
    }

    //------------------------------------------------------------------------------
    // Computes the tangent of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Tan(void) noexcept
    {
        m_X = xmath::Tan(radian{ m_X });
        m_Y = xmath::Tan(radian{ m_Y });
        m_Z = xmath::Tan(radian{ m_Z });
        m_W = xmath::Tan(radian{ m_W });
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arcsine of the vector.
    //
    // Returns:
    //  A vector with the arcsine of each component (in radians).
    //
    inline fvec4 fvec4::AsinCopy(void) const noexcept
    {
        return fvec4(xmath::Asin(m_X).m_Value, xmath::Asin(m_Y).m_Value, xmath::Asin(m_Z).m_Value, xmath::Asin(m_W).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arcsine of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Asin(void) noexcept
    {
        m_X = xmath::Asin(m_X).m_Value;
        m_Y = xmath::Asin(m_Y).m_Value;
        m_Z = xmath::Asin(m_Z).m_Value;
        m_W = xmath::Asin(m_W).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arccosine of the vector.
    //
    // Returns:
    //  A vector with the arccosine of each component (in radians).
    //
    inline fvec4 fvec4::AcosCopy(void) const noexcept
    {
        return fvec4(xmath::Acos(m_X).m_Value, xmath::Acos(m_Y).m_Value, xmath::Acos(m_Z).m_Value, xmath::Acos(m_W).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arccosine of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Acos(void) noexcept
    {
        m_X = xmath::Acos(m_X).m_Value;
        m_Y = xmath::Acos(m_Y).m_Value;
        m_Z = xmath::Acos(m_Z).m_Value;
        m_W = xmath::Acos(m_W).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arctangent of the vector.
    //
    // Returns:
    //  A vector with the arctangent of each component (in radians).
    //
    inline fvec4 fvec4::AtanCopy(void) const noexcept
    {
        return fvec4(xmath::Atan(m_X).m_Value, xmath::Atan(m_Y).m_Value, xmath::Atan(m_Z).m_Value, xmath::Atan(m_W).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arctangent of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Atan(void) noexcept
    {
        m_X = xmath::Atan(m_X).m_Value;
        m_Y = xmath::Atan(m_Y).m_Value;
        m_Z = xmath::Atan(m_Z).m_Value;
        m_W = xmath::Atan(m_W).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arctangent of y/x, considering the quadrant.
    //
    // Parameters:
    //  x - The vector containing x components.
    //
    // Returns:
    //  A vector with the arctangent of each y/x pair (in radians).
    //
    inline fvec4 fvec4::Atan2Copy(const fvec4& x) const noexcept
    {
        return fvec4(xmath::Atan2(m_X, x.m_X).m_Value, xmath::Atan2(m_Y, x.m_Y).m_Value, xmath::Atan2(m_Z, x.m_Z).m_Value, xmath::Atan2(m_W, x.m_W).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arctangent of y/x in-place, considering the quadrant.
    //
    // Parameters:
    //  x - The vector containing x components.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    inline fvec4& fvec4::Atan2(const fvec4& x) noexcept
    {
        *this = Atan2Copy(x);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Swizzle methods for float (HLSL-style, return copy with swizzled components)
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Swizzle macros and instantiations
    //------------------------------------------------------------------------------

    #define SWIZZLE1(func, a)                           \
        inline float fvec4::func() const noexcept {    \
            return this->m_##a;                         \
        }

    SWIZZLE1(x, X)
    SWIZZLE1(y, Y)
    SWIZZLE1(z, Z)
    SWIZZLE1(w, W)

    #undef SWIZZLE1

    #define SWIZZLE2(func, a, b)                        \
        inline fvec2 fvec4::func() const noexcept {   \
            return { this->m_##a, this->m_##b };        \
        }

        SWIZZLE2(xx, X, X)
        SWIZZLE2(xy, X, Y)
        SWIZZLE2(xz, X, Z)
        SWIZZLE2(xw, X, W)
        SWIZZLE2(yx, Y, X)
        SWIZZLE2(yy, Y, Y)
        SWIZZLE2(yz, Y, Z)
        SWIZZLE2(yw, Y, W)
        SWIZZLE2(zx, Z, X)
        SWIZZLE2(zy, Z, Y)
        SWIZZLE2(zz, Z, Z)
        SWIZZLE2(zw, Z, W)
        SWIZZLE2(wx, W, X)
        SWIZZLE2(wy, W, Y)
        SWIZZLE2(wz, W, Z)
        SWIZZLE2(ww, W, W)

    #undef SWIZZLE2

    #define SWIZZLE3(func, a, b, c) \
        inline fvec3 fvec4::func() const noexcept {                                  \
            enum : int { X = 0, Y = 1, Z = 2, W = 3 };                                 \
            return fvec3{ _mm_shuffle_ps( m_XYZW, m_XYZW, _MM_SHUFFLE(W, c, b, a)) }; \
        }

        SWIZZLE3(xxx, X, X, X)
        SWIZZLE3(xxy, X, X, Y)
        SWIZZLE3(xxz, X, X, Z)
        SWIZZLE3(xxw, X, X, W)
        SWIZZLE3(xyx, X, Y, X)
        SWIZZLE3(xyy, X, Y, Y)
        SWIZZLE3(xyz, X, Y, Z)
        SWIZZLE3(xyw, X, Y, W)
        SWIZZLE3(xzx, X, Z, X)
        SWIZZLE3(xzy, X, Z, Y)
        SWIZZLE3(xzz, X, Z, Z)
        SWIZZLE3(xzw, X, Z, W)
        SWIZZLE3(xwx, X, W, X)
        SWIZZLE3(xwy, X, W, Y)
        SWIZZLE3(xwz, X, W, Z)
        SWIZZLE3(xww, X, W, W)
        SWIZZLE3(yxx, Y, X, X)
        SWIZZLE3(yxy, Y, X, Y)
        SWIZZLE3(yxz, Y, X, Z)
        SWIZZLE3(yxw, Y, X, W)
        SWIZZLE3(yyx, Y, Y, X)
        SWIZZLE3(yyy, Y, Y, Y)
        SWIZZLE3(yyz, Y, Y, Z)
        SWIZZLE3(yyw, Y, Y, W)
        SWIZZLE3(yzx, Y, Z, X)
        SWIZZLE3(yzy, Y, Z, Y)
        SWIZZLE3(yzz, Y, Z, Z)
        SWIZZLE3(yzw, Y, Z, W)
        SWIZZLE3(ywx, Y, W, X)
        SWIZZLE3(ywy, Y, W, Y)
        SWIZZLE3(ywz, Y, W, Z)
        SWIZZLE3(yww, Y, W, W)
        SWIZZLE3(zxx, Z, X, X)
        SWIZZLE3(zxy, Z, X, Y)
        SWIZZLE3(zxz, Z, X, Z)
        SWIZZLE3(zxw, Z, X, W)
        SWIZZLE3(zyx, Z, Y, X)
        SWIZZLE3(zyy, Z, Y, Y)
        SWIZZLE3(zyz, Z, Y, Z)
        SWIZZLE3(zyw, Z, Y, W)
        SWIZZLE3(zzx, Z, Z, X)
        SWIZZLE3(zzy, Z, Z, Y)
        SWIZZLE3(zzz, Z, Z, Z)
        SWIZZLE3(zzw, Z, Z, W)
        SWIZZLE3(zwx, Z, W, X)
        SWIZZLE3(zwy, Z, W, Y)
        SWIZZLE3(zwz, Z, W, Z)
        SWIZZLE3(zww, Z, W, W)
        SWIZZLE3(wxx, W, X, X)
        SWIZZLE3(wxy, W, X, Y)
        SWIZZLE3(wxz, W, X, Z)
        SWIZZLE3(wxw, W, X, W)
        SWIZZLE3(wyx, W, Y, X)
        SWIZZLE3(wyy, W, Y, Y)
        SWIZZLE3(wyz, W, Y, Z)
        SWIZZLE3(wyw, W, Y, W)
        SWIZZLE3(wzx, W, Z, X)
        SWIZZLE3(wzy, W, Z, Y)
        SWIZZLE3(wzz, W, Z, Z)
        SWIZZLE3(wzw, W, Z, W)
        SWIZZLE3(wwx, W, W, X)
        SWIZZLE3(wwy, W, W, Y)
        SWIZZLE3(wwz, W, W, Z)
        SWIZZLE3(www, W, W, W)

    #undef SWIZZLE3

#define SWIZZLE4(func, a, b, c, d)                                                \
    inline fvec4 fvec4::func() const noexcept {                                 \
        enum : int { X = 0, Y = 1, Z = 2, W = 3 };                                \
        return fvec4{ _mm_shuffle_ps(m_XYZW, m_XYZW, _MM_SHUFFLE(d, c, b, a)) }; \
    }

        SWIZZLE4(xxxx, X, X, X, X)
        SWIZZLE4(xxxy, X, X, X, Y)
        SWIZZLE4(xxxz, X, X, X, Z)
        SWIZZLE4(xxxw, X, X, X, W)
        SWIZZLE4(xxyx, X, X, Y, X)
        SWIZZLE4(xxyy, X, X, Y, Y)
        SWIZZLE4(xxyz, X, X, Y, Z)
        SWIZZLE4(xxyw, X, X, Y, W)
        SWIZZLE4(xxzx, X, X, Z, X)
        SWIZZLE4(xxzy, X, X, Z, Y)
        SWIZZLE4(xxzz, X, X, Z, Z)
        SWIZZLE4(xxzw, X, X, Z, W)
        SWIZZLE4(xxwx, X, X, W, X)
        SWIZZLE4(xxwy, X, X, W, Y)
        SWIZZLE4(xxwz, X, X, W, Z)
        SWIZZLE4(xxww, X, X, W, W)
        SWIZZLE4(xyxx, X, Y, X, X)
        SWIZZLE4(xyxy, X, Y, X, Y)
        SWIZZLE4(xyxz, X, Y, X, Z)
        SWIZZLE4(xyxw, X, Y, X, W)
        SWIZZLE4(xyyx, X, Y, Y, X)
        SWIZZLE4(xyyy, X, Y, Y, Y)
        SWIZZLE4(xyyz, X, Y, Y, Z)
        SWIZZLE4(xyyw, X, Y, Y, W)
        SWIZZLE4(xyzx, X, Y, Z, X)
        SWIZZLE4(xyzy, X, Y, Z, Y)
        SWIZZLE4(xyzz, X, Y, Z, Z)
        SWIZZLE4(xyzw, X, Y, Z, W)
        SWIZZLE4(xywx, X, Y, W, X)
        SWIZZLE4(xywy, X, Y, W, Y)
        SWIZZLE4(xywz, X, Y, W, Z)
        SWIZZLE4(xyww, X, Y, W, W)
        SWIZZLE4(xzxx, X, Z, X, X)
        SWIZZLE4(xzxy, X, Z, X, Y)
        SWIZZLE4(xzxz, X, Z, X, Z)
        SWIZZLE4(xzxw, X, Z, X, W)
        SWIZZLE4(xzyx, X, Z, Y, X)
        SWIZZLE4(xzyy, X, Z, Y, Y)
        SWIZZLE4(xzyz, X, Z, Y, Z)
        SWIZZLE4(xzyw, X, Z, Y, W)
        SWIZZLE4(xzzx, X, Z, Z, X)
        SWIZZLE4(xzzy, X, Z, Z, Y)
        SWIZZLE4(xzzz, X, Z, Z, Z)
        SWIZZLE4(xzzw, X, Z, Z, W)
        SWIZZLE4(xzwx, X, Z, W, X)
        SWIZZLE4(xzwy, X, Z, W, Y)
        SWIZZLE4(xzwz, X, Z, W, Z)
        SWIZZLE4(xzww, X, Z, W, W)
        SWIZZLE4(xwxx, X, W, X, X)
        SWIZZLE4(xwxy, X, W, X, Y)
        SWIZZLE4(xwxz, X, W, X, Z)
        SWIZZLE4(xwxw, X, W, X, W)
        SWIZZLE4(xwyx, X, W, Y, X)
        SWIZZLE4(xwyy, X, W, Y, Y)
        SWIZZLE4(xwyz, X, W, Y, Z)
        SWIZZLE4(xwyw, X, W, Y, W)
        SWIZZLE4(xwzx, X, W, Z, X)
        SWIZZLE4(xwzy, X, W, Z, Y)
        SWIZZLE4(xwzz, X, W, Z, Z)
        SWIZZLE4(xwzw, X, W, Z, W)
        SWIZZLE4(xwwx, X, W, W, X)
        SWIZZLE4(xwwy, X, W, W, Y)
        SWIZZLE4(xwwz, X, W, W, Z)
        SWIZZLE4(xwww, X, W, W, W)
        SWIZZLE4(yxxx, Y, X, X, X)
        SWIZZLE4(yxxy, Y, X, X, Y)
        SWIZZLE4(yxxz, Y, X, X, Z)
        SWIZZLE4(yxxw, Y, X, X, W)
        SWIZZLE4(yxyx, Y, X, Y, X)
        SWIZZLE4(yxyy, Y, X, Y, Y)
        SWIZZLE4(yxyz, Y, X, Y, Z)
        SWIZZLE4(yxyw, Y, X, Y, W)
        SWIZZLE4(yxzx, Y, X, Z, X)
        SWIZZLE4(yxzy, Y, X, Z, Y)
        SWIZZLE4(yxzz, Y, X, Z, Z)
        SWIZZLE4(yxzw, Y, X, Z, W)
        SWIZZLE4(yxwx, Y, X, W, X)
        SWIZZLE4(yxwy, Y, X, W, Y)
        SWIZZLE4(yxwz, Y, X, W, Z)
        SWIZZLE4(yxww, Y, X, W, W)
        SWIZZLE4(yyxx, Y, Y, X, X)
        SWIZZLE4(yyxy, Y, Y, X, Y)
        SWIZZLE4(yyxz, Y, Y, X, Z)
        SWIZZLE4(yyxw, Y, Y, X, W)
        SWIZZLE4(yyyx, Y, Y, Y, X)
        SWIZZLE4(yyyy, Y, Y, Y, Y)
        SWIZZLE4(yyyz, Y, Y, Y, Z)
        SWIZZLE4(yyyw, Y, Y, Y, W)
        SWIZZLE4(yyzx, Y, Y, Z, X)
        SWIZZLE4(yyzy, Y, Y, Z, Y)
        SWIZZLE4(yyzz, Y, Y, Z, Z)
        SWIZZLE4(yyzw, Y, Y, Z, W)
        SWIZZLE4(yywx, Y, Y, W, X)
        SWIZZLE4(yywy, Y, Y, W, Y)
        SWIZZLE4(yywz, Y, Y, W, Z)
        SWIZZLE4(yyww, Y, Y, W, W)
        SWIZZLE4(yzxx, Y, Z, X, X)
        SWIZZLE4(yzxy, Y, Z, X, Y)
        SWIZZLE4(yzxz, Y, Z, X, Z)
        SWIZZLE4(yzxw, Y, Z, X, W)
        SWIZZLE4(yzyx, Y, Z, Y, X)
        SWIZZLE4(yzyy, Y, Z, Y, Y)
        SWIZZLE4(yzyz, Y, Z, Y, Z)
        SWIZZLE4(yzyw, Y, Z, Y, W)
        SWIZZLE4(yzzx, Y, Z, Z, X)
        SWIZZLE4(yzzy, Y, Z, Z, Y)
        SWIZZLE4(yzzz, Y, Z, Z, Z)
        SWIZZLE4(yzzw, Y, Z, Z, W)
        SWIZZLE4(yzwx, Y, Z, W, X)
        SWIZZLE4(yzwy, Y, Z, W, Y)
        SWIZZLE4(yzwz, Y, Z, W, Z)
        SWIZZLE4(yzww, Y, Z, W, W)
        SWIZZLE4(ywxx, Y, W, X, X)
        SWIZZLE4(ywxy, Y, W, X, Y)
        SWIZZLE4(ywxz, Y, W, X, Z)
        SWIZZLE4(ywxw, Y, W, X, W)
        SWIZZLE4(ywyx, Y, W, Y, X)
        SWIZZLE4(ywyy, Y, W, Y, Y)
        SWIZZLE4(ywyz, Y, W, Y, Z)
        SWIZZLE4(ywyw, Y, W, Y, W)
        SWIZZLE4(ywzx, Y, W, Z, X)
        SWIZZLE4(ywzy, Y, W, Z, Y)
        SWIZZLE4(ywzz, Y, W, Z, Z)
        SWIZZLE4(ywzw, Y, W, Z, W)
        SWIZZLE4(ywwx, Y, W, W, X)
        SWIZZLE4(ywwy, Y, W, W, Y)
        SWIZZLE4(ywwz, Y, W, W, Z)
        SWIZZLE4(ywww, Y, W, W, W)
        SWIZZLE4(zxxx, Z, X, X, X)
        SWIZZLE4(zxxy, Z, X, X, Y)
        SWIZZLE4(zxxz, Z, X, X, Z)
        SWIZZLE4(zxxw, Z, X, X, W)
        SWIZZLE4(zxyx, Z, X, Y, X)
        SWIZZLE4(zxyy, Z, X, Y, Y)
        SWIZZLE4(zxyz, Z, X, Y, Z)
        SWIZZLE4(zxyw, Z, X, Y, W)
        SWIZZLE4(zxzx, Z, X, Z, X)
        SWIZZLE4(zxzy, Z, X, Z, Y)
        SWIZZLE4(zxzz, Z, X, Z, Z)
        SWIZZLE4(zxzw, Z, X, Z, W)
        SWIZZLE4(zxwx, Z, X, W, X)
        SWIZZLE4(zxwy, Z, X, W, Y)
        SWIZZLE4(zxwz, Z, X, W, Z)
        SWIZZLE4(zxww, Z, X, W, W)
        SWIZZLE4(zyxx, Z, Y, X, X)
        SWIZZLE4(zyxy, Z, Y, X, Y)
        SWIZZLE4(zyxz, Z, Y, X, Z)
        SWIZZLE4(zyxw, Z, Y, X, W)
        SWIZZLE4(zyyx, Z, Y, Y, X)
        SWIZZLE4(zyyy, Z, Y, Y, Y)
        SWIZZLE4(zyyz, Z, Y, Y, Z)
        SWIZZLE4(zyyw, Z, Y, Y, W)
        SWIZZLE4(zyzx, Z, Y, Z, X)
        SWIZZLE4(zyzy, Z, Y, Z, Y)
        SWIZZLE4(zyzz, Z, Y, Z, Z)
        SWIZZLE4(zyzw, Z, Y, Z, W)
        SWIZZLE4(zywx, Z, Y, W, X)
        SWIZZLE4(zywy, Z, Y, W, Y)
        SWIZZLE4(zywz, Z, Y, W, Z)
        SWIZZLE4(zyww, Z, Y, W, W)
        SWIZZLE4(zzxx, Z, Z, X, X)
        SWIZZLE4(zzxy, Z, Z, X, Y)
        SWIZZLE4(zzxz, Z, Z, X, Z)
        SWIZZLE4(zzxw, Z, Z, X, W)
        SWIZZLE4(zzyx, Z, Z, Y, X)
        SWIZZLE4(zzyy, Z, Z, Y, Y)
        SWIZZLE4(zzyz, Z, Z, Y, Z)
        SWIZZLE4(zzyw, Z, Z, Y, W)
        SWIZZLE4(zzzx, Z, Z, Z, X)
        SWIZZLE4(zzzy, Z, Z, Z, Y)
        SWIZZLE4(zzzz, Z, Z, Z, Z)
        SWIZZLE4(zzzw, Z, Z, Z, W)
        SWIZZLE4(zzwx, Z, Z, W, X)
        SWIZZLE4(zzwy, Z, Z, W, Y)
        SWIZZLE4(zzwz, Z, Z, W, Z)
        SWIZZLE4(zzww, Z, Z, W, W)
        SWIZZLE4(zwxx, Z, W, X, X)
        SWIZZLE4(zwxy, Z, W, X, Y)
        SWIZZLE4(zwxz, Z, W, X, Z)
        SWIZZLE4(zwxw, Z, W, X, W)
        SWIZZLE4(zwyx, Z, W, Y, X)
        SWIZZLE4(zwyy, Z, W, Y, Y)
        SWIZZLE4(zwyz, Z, W, Y, Z)
        SWIZZLE4(zwyw, Z, W, Y, W)
        SWIZZLE4(zwzx, Z, W, Z, X)
        SWIZZLE4(zwzy, Z, W, Z, Y)
        SWIZZLE4(zwzz, Z, W, Z, Z)
        SWIZZLE4(zwzw, Z, W, Z, W)
        SWIZZLE4(zwwx, Z, W, W, X)
        SWIZZLE4(zwwy, Z, W, W, Y)
        SWIZZLE4(zwwz, Z, W, W, Z)
        SWIZZLE4(zwww, Z, W, W, W)
        SWIZZLE4(wxxx, W, X, X, X)
        SWIZZLE4(wxxy, W, X, X, Y)
        SWIZZLE4(wxxz, W, X, X, Z)
        SWIZZLE4(wxxw, W, X, X, W)
        SWIZZLE4(wxyx, W, X, Y, X)
        SWIZZLE4(wxyy, W, X, Y, Y)
        SWIZZLE4(wxyz, W, X, Y, Z)
        SWIZZLE4(wxyw, W, X, Y, W)
        SWIZZLE4(wxzx, W, X, Z, X)
        SWIZZLE4(wxzy, W, X, Z, Y)
        SWIZZLE4(wxzz, W, X, Z, Z)
        SWIZZLE4(wxzw, W, X, Z, W)
        SWIZZLE4(wxwx, W, X, W, X)
        SWIZZLE4(wxwy, W, X, W, Y)
        SWIZZLE4(wxwz, W, X, W, Z)
        SWIZZLE4(wxww, W, X, W, W)
        SWIZZLE4(wyxx, W, Y, X, X)
        SWIZZLE4(wyxy, W, Y, X, Y)
        SWIZZLE4(wyxz, W, Y, X, Z)
        SWIZZLE4(wyxw, W, Y, X, W)
        SWIZZLE4(wyyx, W, Y, Y, X)
        SWIZZLE4(wyyy, W, Y, Y, Y)
        SWIZZLE4(wyyz, W, Y, Y, Z)
        SWIZZLE4(wyyw, W, Y, Y, W)
        SWIZZLE4(wyzx, W, Y, Z, X)
        SWIZZLE4(wyzy, W, Y, Z, Y)
        SWIZZLE4(wyzz, W, Y, Z, Z)
        SWIZZLE4(wyzw, W, Y, Z, W)
        SWIZZLE4(wywx, W, Y, W, X)
        SWIZZLE4(wywy, W, Y, W, Y)
        SWIZZLE4(wywz, W, Y, W, Z)
        SWIZZLE4(wyww, W, Y, W, W)
        SWIZZLE4(wzxx, W, Z, X, X)
        SWIZZLE4(wzxy, W, Z, X, Y)
        SWIZZLE4(wzxz, W, Z, X, Z)
        SWIZZLE4(wzxw, W, Z, X, W)
        SWIZZLE4(wzyx, W, Z, Y, X)
        SWIZZLE4(wzyy, W, Z, Y, Y)
        SWIZZLE4(wzyz, W, Z, Y, Z)
        SWIZZLE4(wzyw, W, Z, Y, W)
        SWIZZLE4(wzzx, W, Z, Z, X)
        SWIZZLE4(wzzy, W, Z, Z, Y)
        SWIZZLE4(wzzz, W, Z, Z, Z)
        SWIZZLE4(wzzw, W, Z, Z, W)
        SWIZZLE4(wzwx, W, Z, W, X)
        SWIZZLE4(wzwy, W, Z, W, Y)
        SWIZZLE4(wzwz, W, Z, W, Z)
        SWIZZLE4(wzww, W, Z, W, W)
        SWIZZLE4(wwxx, W, W, X, X)
        SWIZZLE4(wwxy, W, W, X, Y)
        SWIZZLE4(wwxz, W, W, X, Z)
        SWIZZLE4(wwxw, W, W, X, W)
        SWIZZLE4(wwyx, W, W, Y, X)
        SWIZZLE4(wwyy, W, W, Y, Y)
        SWIZZLE4(wwyz, W, W, Y, Z)
        SWIZZLE4(wwyw, W, W, Y, W)
        SWIZZLE4(wwzx, W, W, Z, X)
        SWIZZLE4(wwzy, W, W, Z, Y)
        SWIZZLE4(wwzz, W, W, Z, Z)
        SWIZZLE4(wwzw, W, W, Z, W)
        SWIZZLE4(wwwx, W, W, W, X)
        SWIZZLE4(wwwy, W, W, W, Y)
        SWIZZLE4(wwwz, W, W, W, Z)
        SWIZZLE4(wwww, W, W, W, W)

    #undef SWIZZLE4

    //------------------------------------------------------------------------------
    // Operator overloads
    //------------------------------------------------------------------------------

    inline fvec4 fvec4::operator+(const fvec4& other) const noexcept
    {
        return fvec4(_mm_add_ps(m_XYZW, other.m_XYZW));
    }

    inline fvec4 fvec4::operator-(const fvec4& other) const noexcept
    {
        return fvec4(_mm_sub_ps(m_XYZW, other.m_XYZW));
    }

    inline fvec4 fvec4::operator*(float scalar) const noexcept
    {
        return fvec4(_mm_mul_ps(m_XYZW, _mm_set1_ps(scalar)));
    }

    inline fvec4 fvec4::operator/(float scalar) const noexcept
    {
        return fvec4(_mm_div_ps(m_XYZW, _mm_set1_ps(scalar)));
    }

    inline fvec4& fvec4::operator+=(const fvec4& other) noexcept
    {
        m_XYZW = _mm_add_ps(m_XYZW, other.m_XYZW);
        return *this;
    }

    inline fvec4& fvec4::operator-=(const fvec4& other) noexcept
    {
        m_XYZW = _mm_sub_ps(m_XYZW, other.m_XYZW);
        return *this;
    }

    inline fvec4& fvec4::operator*=(float scalar) noexcept
    {
        m_XYZW = _mm_mul_ps(m_XYZW, _mm_set1_ps(scalar));
        return *this;
    }

    inline fvec4& fvec4::operator/=(float scalar) noexcept
    {
        m_XYZW = _mm_div_ps(m_XYZW, _mm_set1_ps(scalar));
        return *this;
    }

    inline bool fvec4::operator==(const fvec4& other) const noexcept
    {
        return _mm_movemask_ps(_mm_cmpeq_ps(m_XYZW, other.m_XYZW)) == 0xF;
    }

    inline bool fvec4::operator!=(const fvec4& other) const noexcept
    {
        return !(*this == other);
    }

    inline float fvec4::operator[](std::int32_t index) const noexcept
    {
        assert(index >= 0 && index < 4);
        return m_Elements[index];
    }

    inline float& fvec4::operator[](std::int32_t index) noexcept
    {
        assert(index >= 0 && index < 4);
        return m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // Friend operators
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Multiplies a scalar by a vector.
    //
    // Parameters:
    //  scalar - The scalar.
    //  v - The vector.
    //
    // Returns:
    //  The scaled vector.
    //
    inline fvec4 operator*(float scalar, const fvec4& v) noexcept
    {
        return v * scalar;
    }

    //------------------------------------------------------------------------------
    // Negates a vector.
    //
    // Parameters:
    //  v - The vector.
    //
    // Returns:
    //  The negated vector.
    //
    inline fvec4 operator-(const fvec4& v) noexcept
    {
        return fvec4(_mm_mul_ps(v.m_XYZW, _mm_set1_ps(-1.f)));
    }
}

