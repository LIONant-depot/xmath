#pragma once
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from individual components.
    //
    // Params:
    //  m_X - X component.
    //  m_Y - Y component.
    //  m_Z - Z component.
    //
    // Notes:
    //  Sets m_W=1 for homogeneous coords; uses _mm_set_ps for SIMD.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(float x, float y, float z) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(0.0f, z, y, x);
        }
        else
        {
            this->m_X = x;
            this->m_Y = y;
            this->m_Z = z;
        }
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // This constructor is meant to add compatibility between both types
    //
    // Params:
    //  value - another type of fvec3_t
    //
    // Notes:
    //  Sets m_W=1 for homogeneous coords; uses _mm_set_ps for SIMD.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const fvec3_t<!V>& other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(1.0f, other.m_Z, other.m_Y, other.m_X);
        }
        else
        {
            this->m_X = other.m_X;
            this->m_Y = other.m_Y;
            this->m_Z = other.m_Z;
        }
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor filling all components with value (m_W=1).
    //
    // Params:
    //  value - Value for m_X/m_Y/z.
    //
    // Notes:
    //  Uses _mm_set1_ps for SIMD; overrides m_W to 1.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(float value) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set1_ps(value);
            this->m_W = 1.0f;
        }
        else
        {
            this->m_X = this->m_Y = this->m_Z = value;
        }
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from pitch and yaw angles.
    //
    // Params:
    //  pitch - Pitch angle.
    //  yaw - Yaw angle.
    //
    // Returns:
    //  Direction vector from angles.
    //
    // Notes:
    //  Uses SinCos; no SIMD here as trig ops scalar.
    //  Asserts valid angles.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(radian pitch, radian yaw) noexcept
    {
        assert(isValid(pitch.m_Value) && isValid(yaw.m_Value));
        float ps, pc, ys, yc;
        SinCos(pitch, ps, pc);
        SinCos(yaw, ys, yc);
        *this = fvec3_t(ys * pc, -ps, yc * pc);
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from SIMD register.
    //
    // Params:
    //  reg - __m128 register.
    //
    // Notes:
    //  Direct copy; assumes valid reg.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const floatx4& reg) noexcept requires V
        : details::f_vec3::simd_data{ .m_XYZW = reg }
    {}

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from fvec2 and z.
    //
    // Params:
    //  other - 2D vector for x/y.
    //  z - Z component.
    //
    // Notes:
    //  Copies x/y from fvec2; sets z; m_W=1 for SIMD.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(const fvec2& other, float z) noexcept
    {
        assert(other.isFinite() && isValid(z));
        this->m_X = other.m_X;
        this->m_Y = other.m_Y;
        this->m_Z = z;
        if constexpr (V) this->m_W = 1.0f;
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from x and fvec2.
    //
    // Params:
    //  x - X component.
    //  other - 2D vector for y/z.
    //
    // Notes:
    //  Sets x; copies y/z from fvec2; m_W=1 for SIMD.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(float x, const fvec2& other) noexcept
    {
        assert(isValid(x) && other.isFinite());
        this->m_X = x;
        this->m_Y = other.m_X;
        this->m_Z = other.m_Y;
        if constexpr (V) this->m_W = 1.0f;
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from span of floats.
    //
    // Params:
    //  Span - Span of 3 floats.
    //
    // Notes:
    //  Copies from span; sets m_W=1 for SIMD.
    //  Asserts span size=3; assumes finite.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(std::span<float> Span) noexcept
    {
        assert(Span.size() == 3);
        this->m_X = Span[0];
        this->m_Y = Span[1];
        this->m_Z = Span[2];
        if constexpr (V) this->m_W = 1.0f;
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from array of doubles.
    //
    // Params:
    //  Conversion - Array of 3 doubles.
    //
    // Notes:
    //  Static cast to float; sets m_W=1 for SIMD.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const std::array<double, 3>& Conversion) noexcept
    {
        this->m_X = static_cast<float>(Conversion[0]);
        this->m_Y = static_cast<float>(Conversion[1]);
        this->m_Z = static_cast<float>(Conversion[2]);
        if constexpr (V) this->m_W = 1.0f;
    }

    //------------------------------------------------------------------------------
    // operator std::array<double,3>
    //------------------------------------------------------------------------------
    //
    // Converts to array of doubles.
    //
    // Returns:
    //  {x, y, z} as doubles.
    //
    // Notes:
    //  Static cast; ignores m_W.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    constexpr fvec3_t<V>::operator std::array<double, 3>() const noexcept
    {
        return { static_cast<double>(this->m_X), static_cast<double>(this->m_Y), static_cast<double>(this->m_Z) };
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Returns zero vector (0,0,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromZero(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,0,0,0 }} };
        else             return { 0,0,0 };
    }

    //------------------------------------------------------------------------------
    // fromOne
    //------------------------------------------------------------------------------
    //
    // Returns one vector (1,1,1).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromOne(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 1,1,1,1 }} };
        else             return { 1,1,1 };
    }

    //------------------------------------------------------------------------------
    // Up
    //------------------------------------------------------------------------------
    //
    // Returns up vector (0,1,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromUp(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,1,0,0 }} };
        else             return { 0,1,0 };
    }

    //------------------------------------------------------------------------------
    // Down
    //------------------------------------------------------------------------------
    //
    // Returns down vector (0,-1,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromDown(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,-1,0,0 }} };
        else             return { 0,-1,0 };
    }

    //------------------------------------------------------------------------------
    // Left
    //------------------------------------------------------------------------------
    //
    // Returns left vector (-1,0,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromLeft(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ -1, 0, 0, 0 }} };
        else             return { -1,0,0 };
    }

    //------------------------------------------------------------------------------
    // Right
    //------------------------------------------------------------------------------
    //
    // Returns right vector (1,0,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromRight(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 1,0,0,0 }} };
        else             return { 1,0,0 };
    }

    //------------------------------------------------------------------------------
    // Forward
    //------------------------------------------------------------------------------
    //
    // Returns forward vector (0,0,1).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromForward(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,0,1,0 }} };
        else             return { 0,0,1 };
    }

    //------------------------------------------------------------------------------
    // Back
    //------------------------------------------------------------------------------
    //
    // Returns back vector (0,0,-1).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromBack(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,0,-1,0 }} };
        else             return { 0,0,-1 };
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
    //  The dot product scalar.
    //
    template <bool V >
    inline float fvec3_t<V>::Dot(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            return _mm_cvtss_f32(_mm_dp_ps(a.m_XYZW, b.m_XYZW, 0xF7));
        }
        else
        {
            return a.m_X * b.m_X + a.m_Y * b.m_Y + a.m_Z * b.m_Z;
        }
    }

    //------------------------------------------------------------------------------
    // Cross
    //------------------------------------------------------------------------------
    //
    // Computes the cross product of two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  The cross product vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Cross(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            __m128 a_yzx = _mm_shuffle_ps(a.m_XYZW, a.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1));
            __m128 b_yzx = _mm_shuffle_ps(b.m_XYZW, b.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1));
            __m128 res = _mm_sub_ps(_mm_mul_ps(a.m_XYZW, b_yzx), _mm_mul_ps(a_yzx, b.m_XYZW));
            return fvec3_t<V>(_mm_shuffle_ps(res, res, _MM_SHUFFLE(3, 0, 2, 1)));
        }
        else
        {
            return fvec3_t<V>(a.m_Y * b.m_Z - a.m_Z * b.m_Y, a.m_Z * b.m_X - a.m_X * b.m_Z, a.m_X * b.m_Y - a.m_Y * b.m_X);
        }
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
    //  The minimum vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Min(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_min_ps(a.m_XYZW, b.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Min(a.m_X, b.m_X), xmath::Min(a.m_Y, b.m_Y), xmath::Min(a.m_Z, b.m_Z));
        }
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
    //  The maximum vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Max(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_max_ps(a.m_XYZW, b.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Max(a.m_X, b.m_X), xmath::Max(a.m_Y, b.m_Y), xmath::Max(a.m_Z, b.m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Computes the linear interpolation between two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //  t - The interpolation factor (0 to 1).
    //
    // Returns:
    //  The interpolated vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Lerp(const fvec3_t& a, const fvec3_t& b, float t) noexcept
    {
        t = xmath::Clamp(t, 0.0f, 1.0f);
        if constexpr (V)
        {
            __m128 t_vec = _mm_set1_ps(t);
            return fvec3_t<V>(_mm_add_ps(a.m_XYZW, _mm_mul_ps(_mm_sub_ps(b.m_XYZW, a.m_XYZW), t_vec)));
        }
        else
        {
            return fvec3_t<V>(a.m_X + t * (b.m_X - a.m_X), a.m_Y + t * (b.m_Y - a.m_Y), a.m_Z + t * (b.m_Z - a.m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Computes the distance between two vectors.
    //
    // Parameters:
    //  a - The first vector.
    //  b - The second vector.
    //
    // Returns:
    //  The distance.
    //
    template <bool V >
    inline float fvec3_t<V>::Distance(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        return (a - b).Length();
    }

    //------------------------------------------------------------------------------
    // Static methods as members
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes the dot product with another vector (static method as member).
    //
    // Parameters:
    //  a - The second vector.
    //
    // Returns:
    //  The dot product scalar.
    //
    template <bool V >
    inline float fvec3_t<V>::Dot(const fvec3_t& a) const noexcept
    {
        return Dot(*this, a);
    }

    //------------------------------------------------------------------------------
    // Cross
    //------------------------------------------------------------------------------
    //
    // Computes the cross product with another vector (static method as member).
    //
    // Parameters:
    //  a - The second vector.
    //
    // Returns:
    //  The cross product vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Cross(const fvec3_t& a) const noexcept
    {
        return Cross(*this, a);
    }

    //------------------------------------------------------------------------------
    // Min
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise minimum with another vector (static method as member).
    //
    // Parameters:
    //  a - The second vector.
    //
    // Returns:
    //  The minimum vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Min(const fvec3_t& a) const noexcept
    {
        return Min(*this, a);
    }

    //------------------------------------------------------------------------------
    // Max
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise maximum with another vector (static method as member).
    //
    // Parameters:
    //  a - The second vector.
    //
    // Returns:
    //  The maximum vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Max(const fvec3_t& a) const noexcept
    {
        return Max(*this, a);
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Computes the linear interpolation with another vector (static method as member).
    //
    // Parameters:
    //  a - The end vector.
    //  t - The interpolation factor (0 to 1).
    //
    // Returns:
    //  The interpolated vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Lerp(const fvec3_t& a, float t) const noexcept
    {
        return Lerp(*this, a, t);
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Computes the distance to another vector (static method as member).
    //
    // Parameters:
    //  a - The target vector.
    //
    // Returns:
    //  The distance.
    //
    template <bool V >
    inline float fvec3_t<V>::Distance(const fvec3_t& a) const noexcept
    {
        return Distance(*this, a);
    }

    //------------------------------------------------------------------------------
    // Instance methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Length
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean length (magnitude) of the vector.
    //
    // Returns:
    //  sqrt(m_X^2 + m_Y^2 + m_Z^2).
    //
    // Notes:
    //  Uses LengthSq and xmath::Sqrt; scalar operation.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline float fvec3_t<V>::Length() const noexcept
    {
        return xmath::Sqrt(LengthSq());
    }

    //------------------------------------------------------------------------------
    // LengthSq
    //------------------------------------------------------------------------------
    //
    // Computes the squared Euclidean length of the vector.
    //
    // Returns:
    //  m_X^2 + m_Y^2 + m_Z^2.
    //
    // Notes:
    //  Faster than Length for comparisons; uses Dot product.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline float fvec3_t<V>::LengthSq() const noexcept
    {
        return Dot(*this, *this);
    }

    //------------------------------------------------------------------------------
    // LimitLengthCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy of the vector with length capped to MaxLength.
    //
    // Parameters:
    //  MaxLength - Maximum allowed length.
    //
    // Returns:
    //  Vector with length <= MaxLength.
    //
    // Notes:
    //  Uses LengthSq and InvSqrt; returns original if length <= MaxLength.
    //  Assumes finite vector and positive MaxLength; asserts in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::LimitLengthCopy(float MaxLength) const noexcept
    {
        float lenSq = LengthSq();
        if (lenSq <= MaxLength * MaxLength) return *this;
        return *this * (MaxLength * xmath::InvSqrt(lenSq));
    }

    //------------------------------------------------------------------------------
    // NormalizeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy of the vector (unit length).
    //
    // Returns:
    //  Vector / Length, or (0,0,0) if length < 1e-6.
    //
    // Notes:
    //  Uses Length and division; scalar operation.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::NormalizeCopy() const noexcept
    {
        float len = Length();
        if (len < 1e-6f) return fvec3_t<V>(0.f, 0.f, 0.f);
        return *this / len;
    }

    //------------------------------------------------------------------------------
    // Normalize
    //------------------------------------------------------------------------------
    //
    // Normalizes the vector in-place to unit length.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Uses NormalizeCopy; sets to (0,0,0) if length < 1e-6.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Normalize() noexcept
    {
        *this = NormalizeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // NormalizeSafeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy or Right() if length is near zero.
    //
    // Returns:
    //  Normalized vector or Right() if length < 1e-6.
    //
    // Notes:
    //  Avoids division by zero; uses Length and division.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::NormalizeSafeCopy() const noexcept
    {
        float len = Length();
        if (len < 1e-6f) return fromRight();
        return *this / len;
    }

    //------------------------------------------------------------------------------
    // NormalizeSafe
    //------------------------------------------------------------------------------
    //
    // Normalizes the vector in-place or sets to Right() if length is near zero.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Uses NormalizeSafeCopy; avoids division by zero.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::NormalizeSafe() noexcept
    {
        *this = NormalizeSafeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all vector components are finite (no NaN or infinity).
    //
    // Returns:
    //  True if all components are finite.
    //
    // Notes:
    //  Uses xmath::isFinite on x, y, z components; ignores m_W.
    //  No assertions as this is a check function.
    //
    template <bool V >
    inline bool fvec3_t<V>::isFinite() const noexcept
    {
        return xmath::isFinite(this->m_X) && xmath::isFinite(this->m_Y) && xmath::isFinite(this->m_Z);
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if all components are within the specified range.
    //
    // Parameters:
    //  min - Lower bound of range.
    //  max - Upper bound of range.
    //
    // Returns:
    //  True if all components are in [min, max].
    //
    // Notes:
    //  Uses xmath::isInRange on x, y, z components; ignores m_W.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline bool fvec3_t<V>::isInRange(float min, float max) const noexcept
    {
        return xmath::isInRange(this->m_X, min, max) && xmath::isInRange(this->m_Y, min, max) && xmath::isInRange(this->m_Z, min, max);
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
    template <bool V >
    inline bool fvec3_t<V>::Equals(const fvec3_t& other, float tolerance) const noexcept
    {
        return xmath::Abs(this->m_X - other.m_X) <= tolerance && xmath::Abs(this->m_Y - other.m_Y) <= tolerance && xmath::Abs(this->m_Z - other.m_Z) <= tolerance;
    }

    //------------------------------------------------------------------------------
    // OneOverCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with the reciprocal of each component.
    //
    // Returns:
    //  Vector with 1/m_X, 1/m_Y, 1/m_Z (m_W unchanged).
    //
    // Notes:
    //  Scalar division; no SIMD optimization for reciprocals.
    //  Assumes non-zero components; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::OneOverCopy() const noexcept
    {
        return fvec3_t<V>(1.f / this->m_X, 1.f / this->m_Y, 1.f / this->m_Z);
    }

    //------------------------------------------------------------------------------
    // OneOver
    //------------------------------------------------------------------------------
    //
    // Computes the reciprocal of each component in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Uses OneOverCopy; scalar division.
    //  Assumes non-zero components; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::OneOver() noexcept
    {
        *this = OneOverCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // AbsCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with the absolute value of each component.
    //
    // Returns:
    //  Vector with |m_X|, |m_Y|, |m_Z| (m_W unchanged).
    //
    // Notes:
    //  SIMD uses _mm_andnot_ps to clear sign bit; scalar uses xmath::Abs.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AbsCopy() const noexcept
    {
        if constexpr (V)
        {
            // 0x80000000 is the sign bit in IEEE 754 float representation
            // We use _mm_andnot_ps to clear the sign bit (i.e., compute absolute value)
            return fvec3_t<V>(_mm_andnot_ps(_mm_castsi128_ps(_mm_set1_epi32(0x80000000)), this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Abs(this->m_X), xmath::Abs(this->m_Y), xmath::Abs(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Abs
    //------------------------------------------------------------------------------
    //
    // Computes the absolute value of each component in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    // Notes:
    //  Uses AbsCopy; SIMD or scalar path.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Abs() noexcept
    {
        *this = AbsCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Reflection
    //------------------------------------------------------------------------------
    //
    // Computes the reflection of the vector over a normal.
    //
    // Parameters:
    //  normal - The normal vector (assumed normalized).
    //
    // Returns:
    //  Reflected vector.
    //
    // Notes:
    //  Formula: v - 2 * (v · n) * n; scalar operation.
    //  Assumes finite vector and normalized normal; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Reflection(const fvec3_t & normal) const noexcept
    {
        float dot = Dot(normal);
        return *this - 2.f * dot * normal;
    }

    //------------------------------------------------------------------------------
    // DistanceSqr
    //------------------------------------------------------------------------------
    //
    // Computes the squared Euclidean distance to another vector.
    //
    // Parameters:
    //  v - The target vector.
    //
    // Returns:
    //  Squared distance.
    //
    // Notes:
    //  Uses LengthSq of difference; faster than Distance for comparisons.
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline float fvec3_t<V>::DistanceSqr(const fvec3_t & v) const noexcept
    {
        return (*this - v).LengthSq();
    }

    //------------------------------------------------------------------------------
    // Pitch
    //------------------------------------------------------------------------------
    //
    // Computes the pitch angle of the vector.
    //
    // Returns:
    //  Pitch angle in radians.
    //
    // Notes:
    //  Uses Atan2 for pitch calculation; scalar operation.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline radian fvec3_t<V>::Pitch() const noexcept
    {
        return radian(xmath::Atan2(-this->m_Y, xmath::Sqrt(this->m_X * this->m_X + this->m_Z * this->m_Z)));
    }

    //------------------------------------------------------------------------------
    // Yaw
    //------------------------------------------------------------------------------
    //
    // Computes the yaw angle of the vector.
    //
    // Returns:
    //  Yaw angle in radians.
    //
    // Notes:
    //  Uses Atan2 for yaw calculation; scalar operation.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline radian fvec3_t<V>::Yaw() const noexcept
    {
        return radian(xmath::Atan2(this->m_X, this->m_Z));
    }

    //------------------------------------------------------------------------------
    // PitchYaw
    //------------------------------------------------------------------------------
    //
    // Computes both pitch and yaw angles of the vector.
    //
    // Returns:
    //  Pair of pitch and yaw angles in radians.
    //
    // Notes:
    //  Calls Pitch and Yaw; scalar operation.
    //  Assumes finite vector; asserts validity in inline.
    //
    template <bool V >
    inline std::pair<radian, radian> fvec3_t<V>::PitchYaw() const noexcept
    {
        return { Pitch(), Yaw() };
    }

    //------------------------------------------------------------------------------
    // RotationTowards
    //------------------------------------------------------------------------------
    //
    // Computes the rotation axis and angle to align with a target vector.
    //
    // Parameters:
    //  dest - The target vector.
    //
    // Returns:
    //  Pair of rotation axis (unit vector) and angle in radians.
    //
    // Notes:
    //  Uses Dot, Cross, and Acos; handles near-zero cases with default axis.
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline std::pair<fvec3_t<V>, radian> fvec3_t<V>::RotationTowards(const fvec3_t & dest) const noexcept
    {
        float cos = xmath::Clamp(Dot(dest) / (Length() * dest.Length()), -1.f, 1.f);
        fvec3_t cross = Cross(dest);
        float sin = cross.Length();
        radian angle(xmath::Acos(cos));
        if (sin < 1e-6f) return { fvec3_t<V>(0.f, 0.f, 1.f), angle };
        return { cross.Normalize(), angle };
    }

    //------------------------------------------------------------------------------
    // AngleBetween
    //------------------------------------------------------------------------------
    //
    // Computes the unsigned angle between this vector and another.
    //
    // Parameters:
    //  v - The other vector.
    //
    // Returns:
    //  Angle in radians [0, pi].
    //
    // Notes:
    //  Uses Dot, Length, and Acos; clamps cosine to [-1, 1].
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline radian fvec3_t<V>::AngleBetween(const fvec3_t & v) const noexcept
    {
        return radian(xmath::Acos(xmath::Clamp(Dot(v) / (Length() * v.Length()), -1.f, 1.f)));
    }

    //------------------------------------------------------------------------------
    // VectorToLineSegment
    //------------------------------------------------------------------------------
    //
    // Computes the vector from this point to the closest point on a line segment.
    //
    // Parameters:
    //  start - The start of the line segment.
    //  end - The end of the line segment.
    //
    // Returns:
    //  Vector from this point to the closest point on the segment.
    //
    // Notes:
    //  Projects onto segment and clamps parameter; scalar operation.
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::VectorToLineSegment(const fvec3_t & start, const fvec3_t & end) const noexcept
    {
        fvec3_t dir = end - start;
        float t = Dot(*this - start, dir) / dir.LengthSq();
        t = xmath::Clamp(t, 0.f, 1.f);
        return (start + t * dir) - *this;
    }

    //------------------------------------------------------------------------------
    // SquareDistToLineSeg
    //------------------------------------------------------------------------------
    //
    // Computes the squared distance to the closest point on a line segment.
    //
    // Parameters:
    //  start - The start of the line segment.
    //  end - The end of the line segment.
    //
    // Returns:
    //  Squared distance to the closest point.
    //
    // Notes:
    //  Uses VectorToLineSegment and LengthSq; faster than Distance for comparisons.
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline float fvec3_t<V>::SquareDistToLineSeg(const fvec3_t & start, const fvec3_t & end) const noexcept
    {
        return VectorToLineSegment(start, end).LengthSq();
    }

    //------------------------------------------------------------------------------
    // ClosestPointInLineSegment
    //------------------------------------------------------------------------------
    //
    // Computes the closest point on a line segment to this vector.
    //
    // Parameters:
    //  start - The start of the line segment.
    //  end - The end of the line segment.
    //
    // Returns:
    //  The closest point on the segment.
    //
    // Notes:
    //  Uses VectorToLineSegment; scalar operation.
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClosestPointInLineSegment(const fvec3_t & start, const fvec3_t & end) const noexcept
    {
        return *this + VectorToLineSegment(start, end);
    }

    //------------------------------------------------------------------------------
    // ClosestPointToRectangle
    //------------------------------------------------------------------------------
    //
    // Computes the closest point and squared distance to a rectangle.
    //
    // Parameters:
    //  p0 - The origin point of the rectangle.
    //  e0 - First edge vector of the rectangle.
    //  e1 - Second edge vector of the rectangle.
    //  outClosestPoint - Output parameter for the closest point.
    //
    // Returns:
    //  Squared distance to the closest point.
    //
    // Notes:
    //  Projects onto rectangle and clamps parameters; scalar operation.
    //  Assumes finite vectors; asserts validity in inline.
    //
    template <bool V >
    inline float fvec3_t<V>::ClosestPointToRectangle(const fvec3_t & p0, const fvec3_t & e0, const fvec3_t & e1, fvec3_t & outClosestPoint) const noexcept
    {
        fvec3_t diff = *this - p0;
        float t0 = Dot(diff, e0) / e0.LengthSq();
        float t1 = Dot(diff, e1) / e1.LengthSq();
        t0 = xmath::Clamp(t0, 0.f, 1.f);
        t1 = xmath::Clamp(t1, 0.f, 1.f);
        outClosestPoint = p0 + t0 * e0 + t1 * e1;
        return (*this - outClosestPoint).LengthSq();
    }
    
    //------------------------------------------------------------------------------
    // RotateXCopy
    //------------------------------------------------------------------------------
    //
    // Rotates vector around X axis in-place.
    //
    // Params:
    //  rx - Rotation angle.
    //
    // Returns:
    //  new rotated vector instance
    //
    // Notes:
    //  Uses SinCos; scalar trig (no SIMD rot matrix yet).
    //  Skips if rx=0; assumes finite/rx valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RotateXCopy(radian rx) const noexcept
    {
        assert(isFinite() && isValid(rx.m_Value));
        if (rx.m_Value == 0.0f) return *this;
        float s, c;
        SinCos(rx, s, c);
        const float y_val = this->m_Y;
        const float z_val = this->m_Z;
        return fvec3_t{ this->m_X, c * y_val - s * z_val, c * z_val + s * y_val };
    }

    //------------------------------------------------------------------------------
    // RotateX
    //------------------------------------------------------------------------------
    //
    // Rotates vector around X axis in-place.
    //
    // Params:
    //  rx - Rotation angle.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses SinCos; scalar trig (no SIMD rot matrix yet).
    //  Skips if rx=0; assumes finite/rx valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateX(radian rx) noexcept
    {
        return *this = RotateXCopy(rx);
    }

    //------------------------------------------------------------------------------
    // RotateYCopy
    //------------------------------------------------------------------------------
    //
    // Rotates vector around Y axis in-place.
    //
    // Params:
    //  ry - Rotation angle.
    //
    // Returns:
    //  new rotated vector instance
    //
    // Notes:
    //  Uses SinCos; scalar.
    //  Skips if ry=0; assumes finite/ry valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RotateYCopy(radian ry) const noexcept
    {
        assert(isFinite() && isValid(ry.m_Value));
        if (ry.m_Value == 0.0f) return *this;
        float s, c;
        SinCos(ry, s, c);
        const float x_val = this->m_X;
        const float z_val = this->m_Z;
        return fvec3_t{ c * x_val + s * z_val, this->m_Y, c * z_val - s * x_val };
    }

    //------------------------------------------------------------------------------
    // RotateY
    //------------------------------------------------------------------------------
    //
    // Rotates vector around Y axis in-place.
    //
    // Params:
    //  ry - Rotation angle.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses SinCos; scalar.
    //  Skips if ry=0; assumes finite/ry valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateY(radian ry) noexcept
    {
        return *this = RotateYCopy(ry);
    }

    //------------------------------------------------------------------------------
    // RotateZCopy
    //------------------------------------------------------------------------------
    //
    // Rotates vector around Z axis in-place.
    //
    // Params:
    //  rz - Rotation angle.
    //
    // Returns:
    //  new rotated vector instance
    //
    // Notes:
    //  Uses SinCos; scalar.
    //  Skips if rz=0; assumes finite/rz valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RotateZCopy(radian rz) const noexcept
    {
        assert(isFinite() && isValid(rz.m_Value));
        if (rz.m_Value == 0.0f) return *this;
        float s, c;
        SinCos(rz, s, c);
        const float x_val = this->m_X;
        const float y_val = this->m_Y;
        return fvec3_t<V>{c* x_val - s * y_val, c* y_val + s * x_val, this->m_Z };
    }

    //------------------------------------------------------------------------------
    // rotateZ
    //------------------------------------------------------------------------------
    //
    // Rotates vector around Z axis in-place.
    //
    // Params:
    //  rz - Rotation angle.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses SinCos; scalar.
    //  Skips if rz=0; assumes finite/rz valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateZ(radian rz) noexcept
    {
        return *this = RotateZCopy(rz);
    }
    //------------------------------------------------------------------------------
    // gridSnap
    //------------------------------------------------------------------------------
    //
    // Snaps components to grid multiples in-place.
    //
    // Params:
    //  gridX - X grid step.
    //  gridY - Y grid step.
    //  gridZ - Z grid step.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses Round; assumes finite/grids >0; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::GridSnap(float gridX, float gridY, float gridZ) noexcept
    {
        assert(isFinite() && gridX > 0.0f && gridY > 0.0f && gridZ > 0.0f);
        this->m_X = Round(this->m_X, gridX);
        this->m_Y = Round(this->m_Y, gridY);
        this->m_Z = Round(this->m_Z, gridZ);
        return *this;
    }

    //------------------------------------------------------------------------------
    // isRightHanded
    //------------------------------------------------------------------------------
    //
    // Checks if points form right-handed 2D triangle (z ignored).
    //
    // Params:
    //  p1 - Second point.
    //  p2 - Third point.
    //
    // Returns:
    //  True if right-handed (negative cross z).
    //
    // Notes:
    //  2D projection; assumes finite; asserts validity.
    //
    template <bool V >
    inline bool fvec3_t<V>::isRightHanded(const fvec3_t& p1, const fvec3_t& p2) const noexcept
    {
        assert(isFinite() && p1.isFinite() && p2.isFinite());
        return ((p1.m_X - this->m_X) * (p2.m_Y - this->m_Y) - (p1.m_Y - this->m_Y) * (p2.m_X - this->m_X)) < 0.0f;
    }

    //------------------------------------------------------------------------------
    // Swizzle methods
    //------------------------------------------------------------------------------
    template <bool V> inline float  fvec3_t<V>::x(void) const noexcept { return this->m_X; }
    template <bool V> inline float  fvec3_t<V>::y(void) const noexcept { return this->m_Y; }
    template <bool V> inline float  fvec3_t<V>::z(void) const noexcept { return this->m_Z; }


    #define SWIZZLE2(func, a, b)                                                                \
        template <bool V> inline fvec2 fvec3_t<V>::func() const noexcept { using namespace xmath; \
            if constexpr (V) return { _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(simde::W, simde::W, simde::##b, simde::##a)) };  \
            else             return fvec2(this->m_##a, this->m_##b);    \
        }

    SWIZZLE2(xx, X, X)
    SWIZZLE2(xy, X, Y)
    SWIZZLE2(xz, X, Z)
    SWIZZLE2(yx, Y, X)
    SWIZZLE2(yy, Y, Y)
    SWIZZLE2(yz, Y, Z)
    SWIZZLE2(zx, Z, X)
    SWIZZLE2(zy, Z, Y)
    SWIZZLE2(zz, Z, Z)

    #undef SWIZZLE2

    #define SWIZZLE3(func, a, b, c)                                                 \
        template <bool V> inline fvec3_t<V> fvec3_t<V>::func() const noexcept { using namespace xmath;           \
            if constexpr (V) return { _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(simde::W,simde::##c,simde::##b,simde::##a)) }; \
            else             return fvec3_t(this->m_##a, this->m_##b, this->m_##c); \
        }

    SWIZZLE3(xxx, X,X,X)
    SWIZZLE3(xxy, X,X,Y)
    SWIZZLE3(xxz, X,X,Z)
    SWIZZLE3(xyx, X,Y,X)
    SWIZZLE3(xyy, X,Y,Y)
    SWIZZLE3(xyz, X,Y,Z)
    SWIZZLE3(xzx, X,Z,X)
    SWIZZLE3(xzy, X,Z,Y)
    SWIZZLE3(xzz, X,Z,Z)
    SWIZZLE3(yxx, Y,X,X)
    SWIZZLE3(yxy, Y,X,Y)
    SWIZZLE3(yxz, Y,X,Z)
    SWIZZLE3(yyx, Y,Y,X)
    SWIZZLE3(yyy, Y,Y,Y)
    SWIZZLE3(yyz, Y,Y,Z)
    SWIZZLE3(yzx, Y,Z,X)
    SWIZZLE3(yzy, Y,Z,Y)
    SWIZZLE3(yzz, Y,Z,Z)
    SWIZZLE3(zxx, Z,X,X)
    SWIZZLE3(zxy, Z,X,Y)
    SWIZZLE3(zxz, Z,X,Z)
    SWIZZLE3(zyx, Z,Y,X)
    SWIZZLE3(zyy, Z,Y,Y)
    SWIZZLE3(zyz, Z,Y,Z)
    SWIZZLE3(zzx, Z,Z,X)
    SWIZZLE3(zzy, Z,Z,Y)
    SWIZZLE3(zzz, Z,Z,Z)

    #undef SWIZZLE3


    #define SWIZZLE4(func, a, b, c, d)                                                              \
        template <bool V> inline fvec4 fvec3_t<V>::func() const noexcept { using namespace xmath;   \
            if constexpr (V) return { _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(simde::##d, simde::##c, simde::##b, simde::##a)) }; \
            else             return fvec4(this->m_##a, this->m_##b, this->m_##c, this->m_##d);      \
        }

    SWIZZLE4(xxxx, X, X, X, X)
    SWIZZLE4(xxxy, X, X, X, Y)
    SWIZZLE4(xxxz, X, X, X, Z)
    SWIZZLE4(xxyx, X, X, Y, X)
    SWIZZLE4(xxyy, X, X, Y, Y)
    SWIZZLE4(xxyz, X, X, Y, Z)
    SWIZZLE4(xxzx, X, X, Z, X)
    SWIZZLE4(xxzy, X, X, Z, Y)
    SWIZZLE4(xxzz, X, X, Z, Z)
    SWIZZLE4(xyxx, X, Y, X, X)
    SWIZZLE4(xyxy, X, Y, X, Y)
    SWIZZLE4(xyxz, X, Y, X, Z)
    SWIZZLE4(xyyx, X, Y, Y, X)
    SWIZZLE4(xyyy, X, Y, Y, Y)
    SWIZZLE4(xyyz, X, Y, Y, Z)
    SWIZZLE4(xyzx, X, Y, Z, X)
    SWIZZLE4(xyzy, X, Y, Z, Y)
    SWIZZLE4(xyzz, X, Y, Z, Z)
    SWIZZLE4(xzxx, X, Z, X, X)
    SWIZZLE4(xzxy, X, Z, X, Y)
    SWIZZLE4(xzxz, X, Z, X, Z)
    SWIZZLE4(xzyx, X, Z, Y, X)
    SWIZZLE4(xzyy, X, Z, Y, Y)
    SWIZZLE4(xzyz, X, Z, Y, Z)
    SWIZZLE4(xzzx, X, Z, Z, X)
    SWIZZLE4(xzzy, X, Z, Z, Y)
    SWIZZLE4(xzzz, X, Z, Z, Z)
    SWIZZLE4(yxxx, Y, X, X, X)
    SWIZZLE4(yxxy, Y, X, X, Y)
    SWIZZLE4(yxxz, Y, X, X, Z)
    SWIZZLE4(yxyx, Y, X, Y, X)
    SWIZZLE4(yxyy, Y, X, Y, Y)
    SWIZZLE4(yxyz, Y, X, Y, Z)
    SWIZZLE4(yxzx, Y, X, Z, X)
    SWIZZLE4(yxzy, Y, X, Z, Y)
    SWIZZLE4(yxzz, Y, X, Z, Z)
    SWIZZLE4(yyxx, Y, Y, X, X)
    SWIZZLE4(yyxy, Y, Y, X, Y)
    SWIZZLE4(yyxz, Y, Y, X, Z)
    SWIZZLE4(yyyx, Y, Y, Y, X)
    SWIZZLE4(yyyy, Y, Y, Y, Y)
    SWIZZLE4(yyyz, Y, Y, Y, Z)
    SWIZZLE4(yyzx, Y, Y, Z, X)
    SWIZZLE4(yyzy, Y, Y, Z, Y)
    SWIZZLE4(yyzz, Y, Y, Z, Z)
    SWIZZLE4(yzxx, Y, Z, X, X)
    SWIZZLE4(yzxy, Y, Z, X, Y)
    SWIZZLE4(yzxz, Y, Z, X, Z)
    SWIZZLE4(yzyx, Y, Z, Y, X)
    SWIZZLE4(yzyy, Y, Z, Y, Y)
    SWIZZLE4(yzyz, Y, Z, Y, Z)
    SWIZZLE4(yzzx, Y, Z, Z, X)
    SWIZZLE4(yzzy, Y, Z, Z, Y)
    SWIZZLE4(yzzz, Y, Z, Z, Z)
    SWIZZLE4(zxxx, Z, X, X, X)
    SWIZZLE4(zxxy, Z, X, X, Y)
    SWIZZLE4(zxxz, Z, X, X, Z)
    SWIZZLE4(zxyx, Z, X, Y, X)
    SWIZZLE4(zxyy, Z, X, Y, Y)
    SWIZZLE4(zxyz, Z, X, Y, Z)
    SWIZZLE4(zxzx, Z, X, Z, X)
    SWIZZLE4(zxzy, Z, X, Z, Y)
    SWIZZLE4(zxzz, Z, X, Z, Z)
    SWIZZLE4(zyxx, Z, Y, X, X)
    SWIZZLE4(zyxy, Z, Y, X, Y)
    SWIZZLE4(zyxz, Z, Y, X, Z)
    SWIZZLE4(zyyx, Z, Y, Y, X)
    SWIZZLE4(zyyy, Z, Y, Y, Y)
    SWIZZLE4(zyyz, Z, Y, Y, Z)
    SWIZZLE4(zyzx, Z, Y, Z, X)
    SWIZZLE4(zyzy, Z, Y, Z, Y)
    SWIZZLE4(zyzz, Z, Y, Z, Z)
    SWIZZLE4(zzxx, Z, Z, X, X)
    SWIZZLE4(zzxy, Z, Z, X, Y)
    SWIZZLE4(zzxz, Z, Z, X, Z)
    SWIZZLE4(zzyx, Z, Z, Y, X)
    SWIZZLE4(zzyy, Z, Z, Y, Y)
    SWIZZLE4(zzyz, Z, Z, Y, Z)
    SWIZZLE4(zzzx, Z, Z, Z, X)
    SWIZZLE4(zzzy, Z, Z, Z, Y)
    SWIZZLE4(zzzz, Z, Z, Z, Z)

    #undef SWIZZLE4


    //------------------------------------------------------------------------------
    // Operator overloads
    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator+(const fvec3_t& other) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_add_ps(this->m_XYZW, other.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(this->m_X + other.m_X, this->m_Y + other.m_Y, this->m_Z + other.m_Z);
        }
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator-(const fvec3_t& other) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_sub_ps(this->m_XYZW, other.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(this->m_X - other.m_X, this->m_Y - other.m_Y, this->m_Z - other.m_Z);
        }
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator*(float scalar) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_mul_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else
        {
            return fvec3_t<V>(this->m_X * scalar, this->m_Y * scalar, this->m_Z * scalar);
        }
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fvec3_t<V> fvec3_t<V>::operator*(const fvec3_t& other) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_mul_ps(this->m_XYZW, other.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(
                this->m_X * other.m_X,
                this->m_Y * other.m_Y,
                this->m_Z * other.m_Z
            );
        }
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator/(float scalar) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else
        {
            return fvec3_t<V>(this->m_X / scalar, this->m_Y / scalar, this->m_Z / scalar);
        }
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator+=(const fvec3_t& other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_add_ps(this->m_XYZW, other.m_XYZW);
        }
        else
        {
            this->m_X += other.m_X; this->m_Y += other.m_Y; this->m_Z += other.m_Z;
        }
        return *this;
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator-=(const fvec3_t& other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_sub_ps(this->m_XYZW, other.m_XYZW);
        }
        else
        {
            this->m_X -= other.m_X; this->m_Y -= other.m_Y; this->m_Z -= other.m_Z;
        }
        return *this;
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator*=(float scalar) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_mul_ps(this->m_XYZW, _mm_set1_ps(scalar));
        }
        else
        {
            this->m_X *= scalar; this->m_Y *= scalar; this->m_Z *= scalar;
        }
        return *this;
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator/=(float scalar) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar));
        }
        else
        {
            this->m_X /= scalar; this->m_Y /= scalar; this->m_Z /= scalar;
        }
        return *this;
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline bool fvec3_t<V>::operator==(const fvec3_t& other) const noexcept
    {
        return this->m_X == other.m_X && this->m_Y == other.m_Y && this->m_Z == other.m_Z;
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline bool fvec3_t<V>::operator!=(const fvec3_t& other) const noexcept
    {
        return !(*this == other);
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline float fvec3_t<V>::operator[](std::int32_t index) const noexcept
    {
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline float& fvec3_t<V>::operator[](std::int32_t index) noexcept
    {
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // Friend operators
    //------------------------------------------------------------------------------

    template <bool V >
    inline fvec3_t<V> operator*(float scalar, const fvec3_t<V>& v) noexcept
    {
        return v * scalar;
    }

    template <bool V >
    inline fvec3_t<V> operator-(const fvec3_t<V>& v) noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_mul_ps(v.m_XYZW, _mm_set1_ps(-1.f)));
        }
        else
        {
            return fvec3_t<V>(-v.m_X, -v.m_Y, -v.m_Z);
        }
    }

    //------------------------------------------------------------------------------
    // Missing functions from fvec2 upgrades
    //------------------------------------------------------------------------------

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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::SmoothStep(float edge0, float edge1) const noexcept
    {
        float tx = xmath::Clamp((this->m_X - edge0) / (edge1 - edge0), 0.f, 1.f);
        float ty = xmath::Clamp((this->m_Y - edge0) / (edge1 - edge0), 0.f, 1.f);
        float tz = xmath::Clamp((this->m_Z - edge0) / (edge1 - edge0), 0.f, 1.f);
        tx = tx * tx * (3.f - 2.f * tx);
        ty = ty * ty * (3.f - 2.f * ty);
        tz = tz * tz * (3.f - 2.f * tz);
        return fvec3_t<V>(tx, ty, tz);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Step(float edge) const noexcept
    {
        return fvec3_t<V>((this->m_X < edge) ? 0.f : 1.f, (this->m_Y < edge) ? 0.f : 1.f, (this->m_Z < edge) ? 0.f : 1.f);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClampCopy(float min_val, float max_val) const noexcept
    {
        return fvec3_t<V>(xmath::Clamp(this->m_X, min_val, max_val), xmath::Clamp(this->m_Y, min_val, max_val), xmath::Clamp(this->m_Z, min_val, max_val));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Clamp(float min_val, float max_val) noexcept
    {
        return *this = ClampCopy(min_val, max_val);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ModCopy(float divisor) const noexcept
    {
        return fvec3_t<V>(xmath::FMod(this->m_X, divisor), xmath::FMod(this->m_Y, divisor), xmath::FMod(this->m_Z, divisor));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Mod(float divisor) noexcept
    {
        return *this = ModCopy(divisor);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::FractCopy() const noexcept
    {
        return fvec3_t<V>(this->m_X - xmath::Floor(this->m_X), this->m_Y - xmath::Floor(this->m_Y), this->m_Z - xmath::Floor(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Fract() noexcept
    {
       return *this = FractCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::CeilCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Ceil(this->m_X), xmath::Ceil(this->m_Y), xmath::Ceil(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Ceil() noexcept
    {
        return *this = CeilCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::FloorCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Floor(this->m_X), xmath::Floor(this->m_Y), xmath::Floor(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Floor() noexcept
    {
        return *this = FloorCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::SignCopy() const noexcept
    {
        return fvec3_t<V>((this->m_X > 0.f) ? 1.f : (this->m_X < 0.f) ? -1.f : 0.f,
            (this->m_Y > 0.f) ? 1.f : (this->m_Y < 0.f) ? -1.f : 0.f,
            (this->m_Z > 0.f) ? 1.f : (this->m_Z < 0.f) ? -1.f : 0.f);
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Sign() noexcept
    {
        return *this = SignCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::InvSqrtCopy() const noexcept
    {
        return fvec3_t<V>(xmath::InvSqrt(this->m_X), xmath::InvSqrt(this->m_Y), xmath::InvSqrt(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::InvSqrt() noexcept
    {
        return *this = InvSqrtCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::SqrtCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Sqrt(this->m_X), xmath::Sqrt(this->m_Y), xmath::Sqrt(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Sqrt() noexcept
    {
        return *this = SqrtCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::PowCopy(float exp) const noexcept
    {
        return fvec3_t<V>(xmath::Pow(this->m_X, exp), xmath::Pow(this->m_Y, exp), xmath::Pow(this->m_Z, exp));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Pow(float exp) noexcept
    {
        return *this = PowCopy(exp);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Log2Copy() const noexcept
    {
        return fvec3_t<V>(xmath::Log2(this->m_X), xmath::Log2(this->m_Y), xmath::Log2(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Log2() noexcept
    {
        return *this = Log2Copy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::LogCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Log(this->m_X), xmath::Log(this->m_Y), xmath::Log(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Log() noexcept
    {
        return *this = LogCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::TanCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Tan(this->m_X), xmath::Tan(this->m_Y), xmath::Tan(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Tan() noexcept
    {
        return *this = TanCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AsinCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Asin(this->m_X), xmath::Asin(this->m_Y), xmath::Asin(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Asin() noexcept
    {
        return *this = AsinCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AcosCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Acos(this->m_X), xmath::Acos(this->m_Y), xmath::Acos(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Acos() noexcept
    {
        return *this = AcosCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AtanCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Atan(this->m_X), xmath::Atan(this->m_Y), xmath::Atan(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Atan() noexcept
    {
        return *this = AtanCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Atan2Copy(const fvec3_t& x) const noexcept
    {
        return fvec3_t<V>(xmath::Atan2(this->m_X, x.m_X), xmath::Atan2(this->m_Y, x.m_Y), xmath::Atan2(this->m_Z, x.m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Atan2(const fvec3_t& x) noexcept
    {
        return *this = Atan2Copy(x);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ExpCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Exp(this->m_X), xmath::Exp(this->m_Y), xmath::Exp(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Exp() noexcept
    {
        return *this = ExpCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RoundCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Round(this->m_X), xmath::Round(this->m_Y), xmath::Round(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Round() noexcept
    {
        return *this = RoundCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::TruncCopy() const noexcept
    {
        return fvec3_t<V>(xmath::Trunc(this->m_X), xmath::Trunc(this->m_Y), xmath::Trunc(this->m_Z));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Trunc() noexcept
    {
        return *this = TruncCopy();
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
    template <bool V >
    inline radian fvec3_t<V>::SignedAngleBetween(const fvec3_t& v) const noexcept
    {
        return xmath::Atan2(Cross(v).Length(), Dot(v));
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RotateCopy(const radian3& r) const noexcept
    {
        assert(isFinite());
        RotateZ(r.m_Roll);
        RotateX(r.m_Pitch);
        RotateY(r.m_Yaw);
        return *this;
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Rotate(const radian3& r) noexcept
    {
        return *this = RotateCopy(r);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ProjectCopy(const fvec3_t& onto) const noexcept
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Project(const fvec3_t& onto) noexcept
    {
        return *this = ProjectCopy(onto);
    }

    //------------------------------------------------------------------------------
    // PerpCopy
    //------------------------------------------------------------------------------
    //
    // Computes a vector perpendicular to this one in the plane defined by the normal.
    //
    // Parameters:
    //  normal - The normal vector defining the plane.
    //
    // Returns:
    //  A perpendicular vector (cross product with normal, unnormalized).
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Perpendicular(const fvec3_t& normal) const noexcept
    {
        return Cross(normal);
    }

    //------------------------------------------------------------------------------
    // ProjectOntoPlane
    //
    // Projects this vector onto a plane defined by its normal.
    // Removes the component in the direction of the normal.
    //
    // Parameters:
    //    normal - The normal vector of the plane.
    //
    // Returns:
    //    A vector projected onto the plane.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ProjectOntoPlane(const fvec3_t& normal) const noexcept
    {
        return *this - normal * this->Dot(normal);
    }

    //------------------------------------------------------------------------------
    // ClampCopy
    //
    // Clamps each component of the vector between the corresponding min and max vectors.
    //
    // Parameters:
    //    min - The vector representing the minimum bound for each component.
    //    max - The vector representing the maximum bound for each component.
    //
    // Returns:
    //    A vector where each component is clamped between min and max.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClampCopy(const fvec3_t& min, const fvec3_t& max) const noexcept
    {
        return fvec3_t(
            xmath::Clamp(this->m_X, min.m_X, max.m_X),
            xmath::Clamp(this->m_Y, min.m_Y, max.m_Y),
            xmath::Clamp(this->m_Z, min.m_Z, max.m_Z)
            );
    }

    //------------------------------------------------------------------------------
    // Clamp
    //
    // Clamps each component of the vector between the corresponding min and max vectors.
    //
    // Parameters:
    //    min - The vector representing the minimum bound for each component.
    //    max - The vector representing the maximum bound for each component.
    //
    // Returns:
    //    the modified self
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Clamp(const fvec3_t& min, const fvec3_t& max) noexcept
    {
        return *this = ClampCopy(min,max);
    }

    //------------------------------------------------------------------------------
    // isNearlyZero
    //
    // Determines if all components of the vector are close to zero within a tolerance.
    //
    // Parameters:
    //    tolerance - The tolerance threshold. Defaults to 1e-6.
    //
    // Returns:
    //    True if all components are within [-tolerance, tolerance]; otherwise, false.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline bool fvec3_t<V>::isNearlyZero(float tolerance) const noexcept
    {
        return xmath::Abs(this->m_X) < tolerance &&
                xmath::Abs(this->m_Y) < tolerance &&
                xmath::Abs(this->m_Z) < tolerance;
    }

    //------------------------------------------------------------------------------
    // isNormalized
    //
    // Checks whether the vector is normalized (length approximately 1).
    //
    // Parameters:
    //    tolerance - The acceptable deviation from 1.0. Defaults to 1e-6.
    //
    // Returns:
    //    True if length is within [1 - tolerance, 1 + tolerance]; otherwise, false.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline bool fvec3_t<V>::isNormalized(float tolerance) const noexcept
    {
        const float lenSq = this->LengthSq();
        return xmath::Abs(lenSq - 1.0f) < (2.0f * tolerance);
    }

    //------------------------------------------------------------------------------
    // MoveTowardsCopy
    //
    // Moves this vector towards a target vector by a maximum distance delta.
    //
    // Parameters:
    //    target            - The destination vector.
    //    maxDistanceDelta  - The maximum distance to move.
    //
    // Returns:
    //    The vector moved towards target without overshooting.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::MoveTowardsCopy(const fvec3_t& target, float maxDistanceDelta) const noexcept
    {
        const fvec3_t delta = target - *this;
        const float distance = delta.Length();
        if (distance <= maxDistanceDelta || distance < 1e-6f)
            return target;
        return *this + delta * (maxDistanceDelta / distance);
    }

    //------------------------------------------------------------------------------
    // MoveTowards
    //
    // Moves this vector towards a target vector by a maximum distance delta.
    //
    // Parameters:
    //    target            - The destination vector.
    //    maxDistanceDelta  - The maximum distance to move.
    //
    // Returns:
    //    The vector moved towards target without overshooting.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::MoveTowards(const fvec3_t& target, float maxDistanceDelta) noexcept
    {
        return *this = MoveTowardsCopy(target, maxDistanceDelta );
    }

    //------------------------------------------------------------------------------
    // RandomFloat
    //
    // Generates a random floating-point value between min and max using std::mt19937.
    //
    // Parameters:
    //    min - The minimum value (inclusive).
    //    max - The maximum value (inclusive).
    //
    // Returns:
    //    A random float between min and max.
    //
    // Notes:
    //    Thread-safe if you use thread_local rng.
    //------------------------------------------------------------------------------
    inline static float RandomFloat(float min, float max) noexcept
    {
        static thread_local std::mt19937 generator(std::random_device{}());
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(generator);
    }

    //------------------------------------------------------------------------------
    // RandomUnitVector
    //
    // Generates a random unit vector uniformly distributed on the unit sphere.
    //
    // Returns:
    //    A random unit vector.
    //
    // Note:
    //    Requires access to a float random generator between [-1, 1].
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RandomUnitVector() noexcept
    {
        const float  z = RandomFloat(-1.0f, 1.0f);
        const radian t = xmath::radian{ RandomFloat(0.0f, xmath::pi2_v.m_Value) };
        const float  r = xmath::Sqrt(1.0f - z * z);
        return fvec3_t(xmath::Cos(t) * r, xmath::Sin(t) * r, z);
    }

    //------------------------------------------------------------------------------
    // ToString
    //
    // Converts the vector into a string representation in the form "(x, y, z)".
    //
    // Returns:
    //    A std::string representing the vector components.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    std::string fvec3_t<V>::ToString() const noexcept
    {
        return std::format("({}, {}, {})", this->m_X, this->m_Y, this->m_Z);
    }

    //------------------------------------------------------------------------------
    // operator std::string
    //
    // Allows implicit conversion of the vector into a string representation,
    // formatted as "(x, y, z)".
    //
    // Returns:
    //    A std::string representing the vector.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>::operator std::string() const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // operator<<
    //
    // Overloads the stream output operator to print the vector in "(x, y, z)" format.
    //
    // Parameters:
    //    os   - The output stream.
    //    vec  - The vector to print.
    //
    // Returns:
    //    Reference to the output stream.
    //
    //------------------------------------------------------------------------------
    template <bool V >
    inline std::ostream& operator<<(std::ostream& os, const fvec3_t<V>& vec) noexcept
    {
        return os << '(' << vec.x() << ", " << vec.y() << ", " << vec.z() << ')';
    }

} // namespace xmath



/*************
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
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from individual components.
    //
    // Params:
    //  m_X - X component.
    //  m_Y - Y component.
    //  m_Z - Z component.
    //
    // Notes:
    //  Sets m_W=1 for homogeneous coords; uses _mm_set_ps for SIMD.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(float x, float y, float z) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(0.0f, z, y, x);
        }
        else
        {
            this->m_X = x;
            this->m_Y = y;
            this->m_Z = z;
        }
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // This constructor is meant to add compatibility between both types
    //
    // Params:
    //  value - another type of fvec3_t
    //
    // Notes:
    //  Sets m_W=1 for homogeneous coords; uses _mm_set_ps for SIMD.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const fvec3_t<!V>& other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(1.0f, other.m_Z, other.m_Y, other.m_X);
        }
        else
        {
            this->m_X = other.m_X;
            this->m_Y = other.m_Y;
            this->m_Z = other.m_Z;
        }
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor filling all components with value (m_W=1).
    //
    // Params:
    //  value - Value for m_X/m_Y/z.
    //
    // Notes:
    //  Uses _mm_set1_ps for SIMD; overrides m_W to 1.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(float value) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set1_ps(value);
            this->m_W    = 1.0f;
        }
        else
        {
            this->m_X = this->m_Y = this->m_Z = value;
        }
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from pitch and yaw angles.
    //
    // Params:
    //  pitch - Pitch angle.
    //  yaw - Yaw angle.
    //
    // Returns:
    //  Direction vector from angles.
    //
    // Notes:
    //  Uses SinCos; no SIMD here as trig ops scalar.
    //  Asserts valid angles.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(radian pitch, radian yaw) noexcept
    {
        assert(isValid(pitch.m_Value) && isValid(yaw.m_Value));
        float ps, pc, ys, yc;
        SinCos(pitch, ps, pc);
        SinCos(yaw, ys, yc);
        *this = fvec3_t(ys * pc, -ps, yc * pc);
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from SIMD register.
    //
    // Params:
    //  reg - __m128 register.
    //
    //
    // Notes:
    //  Direct copy; assumes valid reg.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const floatx4& reg) noexcept requires V
    {
        this->m_XYZW = reg;
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Returns zero vector (0,0,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromZero(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,0,0,0 }} };
        else             return { 0,0,0 };
    }

    //------------------------------------------------------------------------------
    // fromOne
    //------------------------------------------------------------------------------
    //
    // Returns one vector (1,1,1).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::fromOne(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 1,1,1,1 }} };
        else             return { 1,1,1 };
    }

    //------------------------------------------------------------------------------
    // Up
    //------------------------------------------------------------------------------
    //
    // Returns up vector (0,1,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::Up(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,1,0,0 }} };
        else             return { 0,1,0 };
    }

    //------------------------------------------------------------------------------
    // Down
    //------------------------------------------------------------------------------
    //
    // Returns down vector (0,-1,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::Down(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,-1,0,0 }} };
        else             return { 0,-1,0 };
    }

    //------------------------------------------------------------------------------
    // Left
    //------------------------------------------------------------------------------
    //
    // Returns left vector (-1,0,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::Left(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ -1, 0, 0, 0 }} };
        else             return { -1,0,0 };
    }

    //------------------------------------------------------------------------------
    // Right
    //------------------------------------------------------------------------------
    //
    // Returns right vector (1,0,0).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::Right(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 1,0,0,0 }} };
        else             return { 1,0,0 };
    }

    //------------------------------------------------------------------------------
    // Forward
    //------------------------------------------------------------------------------
    //
    // Returns forward vector (0,0,1).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::Forward(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,0,1,0 }} };
        else             return { 0,0,1 };
    }

    //------------------------------------------------------------------------------
    // Back
    //------------------------------------------------------------------------------
    //
    // Returns back vector (0,0,-1).
    //
    template <bool V >
    consteval fvec3_t<V> fvec3_t<V>::Back(void) noexcept
    {
        if constexpr (V) return fvec3_t{ floatx4{.m128_f32{ 0,0,-1,0 }} };
        else             return { 0,0,-1 };
    }

    //------------------------------------------------------------------------------
    // Static methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // dot
    //------------------------------------------------------------------------------
    //
    // Computes dot product of two vectors.
    //
    // Params:
    //  a - First vector.
    //  b - Second vector.
    //
    // Returns:
    //  a · b (scalar).
    //
    // Notes:
    //  SIMD uses _mm_dp_ps (ignores m_W); scalar fallback.
    //  Assumes finite vectors; asserts validity.
    //
    template <bool V >
    inline float fvec3_t<V>::Dot(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        assert(a.isFinite() && b.isFinite());
        if constexpr (V) 
        {
            return _mm_cvtss_f32(_mm_dp_ps(a.m_XYZW, b.m_XYZW, 0x77));
        }
        else 
        {
            return a.m_X * b.m_X + a.m_Y * b.m_Y + a.m_Z * b.m_Z;
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline float fvec3_t<V>::Dot(const fvec3_t& a) const noexcept
    {
        return Dot(*this, a);
    }

    //------------------------------------------------------------------------------
    // cross
    //------------------------------------------------------------------------------
    //
    // Computes cross product of two vectors.
    //
    // Params:
    //  a - First vector.
    //  b - Second vector.
    //
    // Returns:
    //  a × b (perpendicular vector).
    //
    // Notes:
    //  SIMD uses shuffles/sub/mul; scalar fallback.
    //  Assumes finite vectors; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Cross(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        assert(a.isFinite() && b.isFinite());
        if constexpr (V) 
        {
            __m128 a_yzx = _mm_shuffle_ps(a.m_XYZW, a.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1));
            __m128 b_yzx = _mm_shuffle_ps(b.m_XYZW, b.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1));
            __m128 res = _mm_sub_ps(_mm_mul_ps(a.m_XYZW, b_yzx), _mm_mul_ps(a_yzx, b.m_XYZW));
            return fvec3_t(_mm_shuffle_ps(res, res, _MM_SHUFFLE(3, 0, 2, 1)));
        }
        else 
        {
            return fvec3_t
            ( a.m_Y * b.m_Z - a.m_Z * b.m_Y
            , a.m_Z * b.m_X - a.m_X * b.m_Z
            , a.m_X * b.m_Y - a.m_Y * b.m_X
            );
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Cross(const fvec3_t& a) const noexcept
    {
        return Cross(*this, a);
    }

    //------------------------------------------------------------------------------
    // min
    //------------------------------------------------------------------------------
    //
    // Component-wise minimum of two vectors.
    //
    // Params:
    //  a - First vector.
    //  b - Second vector.
    //
    // Returns:
    //  Vector with min components.
    //
    // Notes:
    //  SIMD uses _mm_min_ps; scalar fallback.
    //  Assumes finite vectors; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Min(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        assert(a.isFinite() && b.isFinite());
        if constexpr (V) 
        {
            return fvec3_t(_mm_min_ps(a.m_XYZW, b.m_XYZW));
        }
        else 
        {
            return fvec3_t
            ( xmath::Min(a.m_X, b.m_X)
            , xmath::Min(a.m_Y, b.m_Y)
            , xmath::Min(a.m_Z, b.m_Z)
            );
        }
    }

    //------------------------------------------------------------------------------
    // max
    //------------------------------------------------------------------------------
    //
    // Component-wise maximum of two vectors.
    //
    // Params:
    //  a - First vector.
    //  b - Second vector.
    //
    // Returns:
    //  Vector with max components.
    //
    // Notes:
    //  SIMD uses _mm_max_ps; scalar fallback.
    //  Assumes finite vectors; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Max(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        assert(a.isFinite() && b.isFinite());
        if constexpr (V) 
        {
            return fvec3_t(_mm_max_ps(a.m_XYZW, b.m_XYZW));
        }
        else 
        {
            return fvec3_t
            ( xmath::Max(a.m_X, b.m_X)
            , xmath::Max(a.m_Y, b.m_Y)
            , xmath::Max(a.m_Z, b.m_Z)
            );
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Max(const fvec3_t& a) const noexcept
    {
        return Max(*this, a );
    }

    //------------------------------------------------------------------------------
    // lerp
    //------------------------------------------------------------------------------
    //
    // Linear interpolation between two vectors.
    //
    // Params:
    //  a - Start vector.
    //  b - End vector.
    //  t - Factor [0,1].
    //
    // Returns:
    //  Interpolated vector.
    //
    // Notes:
    //  SIMD uses mul/add; scalar fallback.
    //  Clamps t; asserts finite inputs and t range.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Lerp(const fvec3_t& a, const fvec3_t& b, float t) noexcept
    {
        assert(a.isFinite() && b.isFinite() && isValid(t) && xmath::isInRange(t, 0.0f, 1.0f));
        if constexpr (V) 
        {
            t = Range(t, 0.0f, 1.0f);
            __m128 t_vec = _mm_set1_ps(t);
            __m128 one_minus_t = _mm_set1_ps(1.0f - t);
            return fvec3_t(_mm_add_ps(_mm_mul_ps(a.m_XYZW, one_minus_t), _mm_mul_ps(b.m_XYZW, t_vec)));
        }
        else 
        {
            t = Range(t, 0.0f, 1.0f);
            return fvec3_t
            ( a.m_X + (b.m_X - a.m_X) * t
            , a.m_Y + (b.m_Y - a.m_Y) * t
            , a.m_Z + (b.m_Z - a.m_Z) * t
            );
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Lerp(const fvec3_t& a, float t) const noexcept
    {
        return Lerp(*this, a, t);
    }

    //------------------------------------------------------------------------------
    // distance
    //------------------------------------------------------------------------------
    //
    // Computes Euclidean distance between two vectors.
    //
    // Params:
    //  a - First vector.
    //  b - Second vector.
    //
    // Returns:
    //  Distance scalar >=0.
    //
    // Notes:
    //  Uses length of difference; SIMD via length.
    //  Assumes finite vectors; asserts validity.
    //
    template <bool V >
    inline float fvec3_t<V>::Distance(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        assert(a.isFinite() && b.isFinite());
        return (a - b).Length();
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline float fvec3_t<V>::Distance(const fvec3_t& a) const noexcept
    {
        return Distance( *this, a );
    }

    //------------------------------------------------------------------------------
    // Instance methods
    //------------------------------------------------------------------------------

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
    template <bool V >
    inline float fvec3_t<V>::Length(void) const noexcept
    {
        assert(isFinite());
        if constexpr (V) 
        {
            return std::sqrt(_mm_cvtss_f32(_mm_dp_ps(this->m_XYZW, this->m_XYZW, 0x77)));
        }
        else 
        {
            return std::sqrt( this->m_X * this->m_X 
                            + this->m_Y * this->m_Y 
                            + this->m_Z * this->m_Z
                            );
        }
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::LimitLengthCopy(float MaxLength) const noexcept
    {
        auto l = LengthSq();
        if (l <= (MaxLength * MaxLength)) return *this;
        return (*this) * (MaxLength * xmath::InvSqrt(l));
    }

    //------------------------------------------------------------------------------
    // lengthsq
    //------------------------------------------------------------------------------
    //
    // Computes squared magnitude of vector.
    //
    // Returns:
    //  m_X^2 + m_Y^2 + z^2 >=0.
    //
    // Notes:
    //  Faster than length for comparisons; SIMD dot.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline float fvec3_t<V>::LengthSq(void) const noexcept
    {
        assert(isFinite());
        if constexpr (V) 
        {
            return _mm_cvtss_f32(_mm_dp_ps(this->m_XYZW, this->m_XYZW, 0x77));
        }
        else 
        {
            return this->m_X * this->m_X
                 + this->m_Y * this->m_Y
                 + this->m_Z * this->m_Z;
        }
    }

    //------------------------------------------------------------------------------
    // normalize_copy
    //------------------------------------------------------------------------------
    //
    // Returns normalized copy (unit length).
    //
    // Returns:
    //  Vector / length (zero if length <1e-6).
    //
    // Notes:
    //  SIMD div; scalar fallback.
    //  Assumes finite/non-zero; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::NormalizeCopy(void) const noexcept
    {
        assert(isFinite());
        const float len = Length();
        if (len < 1e-6f) return fromZero();
        if constexpr (V) 
        {
            return fvec3_t(_mm_div_ps(this->m_XYZW, _mm_set1_ps(len)));
        }
        else 
        {
            return fvec3_t(this->m_X / len, this->m_Y / len, this->m_Z / len);
        }
    }

    //------------------------------------------------------------------------------
    // normalize
    //------------------------------------------------------------------------------
    //
    // Normalizes vector in-place (unit length).
    //
    // Returns:
    //  Reference to self.
    //
    // Notes:
    //  Uses normalize_copy; chainable.
    //  Assumes finite/non-zero; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Normalize(void) noexcept
    {
        assert(isFinite());
        *this = NormalizeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // normalize_safe_copy
    //------------------------------------------------------------------------------
    //
    // Returns normalized copy or safe default if near zero.
    //
    // Returns:
    //  Normalized or right() if length <1e-6.
    //
    // Notes:
    //  Avoids div-by-zero; SIMD/scalar.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::NormalizeSafeCopy(void) const noexcept
    {
        assert(isFinite());
        const float len = Length();
        if (len < 1e-6f) return Right();
        if constexpr (V) 
        {
            return fvec3_t(_mm_div_ps(this->m_XYZW, _mm_set1_ps(len)));
        }
        else 
        {
            return fvec3_t(this->m_X / len, this->m_Y / len, this->m_Z / len);
        }
    }

    //------------------------------------------------------------------------------
    // normalize_safe
    //------------------------------------------------------------------------------
    //
    // Normalizes in-place or sets to safe default if near zero.
    //
    // Returns:
    //  Reference to self.
    //
    // Notes:
    //  Uses normalize_safe_copy; chainable.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::NormalizeSafe(void) noexcept
    {
        assert(isFinite());
        *this = NormalizeSafeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if vector components are finite (no NaN/inf).
    //
    // Returns:
    //  True if all finite.
    //
    // Notes:
    //  Uses isValid on components; ignores m_W.
    //  No assert as check func.
    //
    template <bool V >
    inline bool fvec3_t<V>::isFinite(void) const noexcept
    {
        return isValid(this->m_X) && isValid(this->m_Y) && isValid(this->m_Z);
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if all components in [min, max].
    //
    // Params:
    //  min - Lower bound.
    //  max - Upper bound.
    //
    // Returns:
    //  True if all in range.
    //
    // Notes:
    //  Uses isInRange on components; ignores m_W.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline bool fvec3_t<V>::isInRange(float min, float max) const noexcept
    {
        assert(isFinite());
        return xmath::isInRange(this->m_X, min, max) && xmath::isInRange(this->m_Y, min, max) && xmath::isInRange(this->m_Z, min, max);
    }

    //------------------------------------------------------------------------------
    // oneOver_copy
    //------------------------------------------------------------------------------
    //
    // Returns reciprocal copy (1/v component-wise).
    //
    // Returns:
    //  Vector with 1/m_X, 1/m_Y, 1/z (m_W unchanged).
    //
    // Notes:
    //  SIMD _mm_rcp_ps; scalar fallback.
    //  Assumes non-zero components; asserts validity/no zero.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::OneOverCopy(void) const noexcept
    {
        assert(isFinite() && this->m_X != 0.0f && this->m_Y != 0.0f && this->m_Z != 0.0f);
        if constexpr (V) 
        {
            return fvec3_t(_mm_rcp_ps(this->m_XYZW));
        }
        else 
        {
            return fvec3_t(1.0f / this->m_X, 1.0f / this->m_Y, 1.0f / this->m_Z);
        }
    }

    //------------------------------------------------------------------------------
    // oneOver
    //------------------------------------------------------------------------------
    //
    // Computes reciprocal in-place (1/v component-wise).
    //
    // Returns:
    //  Reference to self.
    //
    // Notes:
    //  Uses oneOver_copy; chainable.
    //  Assumes non-zero; asserts validity/no zero.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::OneOver(void) noexcept
    {
        assert(isFinite() && this->m_X != 0.0f && this->m_Y != 0.0f && this->m_Z != 0.0f);
        *this = OneOverCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // abs_copy
    //------------------------------------------------------------------------------
    //
    // Returns absolute value copy component-wise.
    //
    // Returns:
    //  |m_X|, |m_Y|, |z| (m_W unchanged).
    //
    // Notes:
    //  SIMD mask; scalar Abs.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AbsCopy(void) const noexcept
    {
        assert(isFinite());
        if constexpr (V) 
        {
            return fvec3_t(_mm_andnot_ps(_mm_set1_ps(-0.0f), this->m_XYZW));
        }
        else 
        {
            return fvec3_t(xmath::Abs(this->m_X), xmath::Abs(this->m_Y), xmath::Abs(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // abs
    //------------------------------------------------------------------------------
    //
    // Computes absolute value in-place component-wise.
    //
    // Returns:
    //  Reference to self.
    //
    // Notes:
    //  Uses abs_copy; chainable.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Abs(void) noexcept
    {
        assert(isFinite());
        return *this = AbsCopy();
    }

    //------------------------------------------------------------------------------
    // reflection
    //------------------------------------------------------------------------------
    //
    // Computes reflection over normal.
    //
    // Params:
    //  normal - Reflection normal (unit vector assumed).
    //
    // Returns:
    //  Reflected vector.
    //
    // Notes:
    //  Formula: v - 2 * (v · n) * n.
    //  Assumes normalized normal; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Reflection(const fvec3_t& normal) const noexcept
    {
        assert(isFinite() && normal.isFinite() && FEqual(normal.LengthSq(), 1.0f));
        if constexpr (V) 
        {
            __m128 dot_val = _mm_dp_ps(this->m_XYZW, normal.m_XYZW, 0x77 );
            __m128 proj    = _mm_mul_ps(_mm_set1_ps(2.0f), _mm_mul_ps(dot_val, normal.m_XYZW));
            return fvec3_t(_mm_sub_ps(this->m_XYZW, proj));
        }
        else 
        {
            const float d = Dot(*this, normal);
            return *this - 2.0f * d * normal;
        }
    }

    //------------------------------------------------------------------------------
    // distanceSquare
    //------------------------------------------------------------------------------
    //
    // Computes squared distance to v.
    //
    // Params:
    //  v - Target vector.
    //
    // Returns:
    //  Squared Euclidean distance.
    //
    // Notes:
    //  Faster than distance for comparisons.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline float fvec3_t<V>::DistanceSquare(const fvec3_t& v) const noexcept
    {
        assert(isFinite() && v.isFinite());
        return (*this - v).LengthSq();
    }

    //------------------------------------------------------------------------------
    // pitch
    //------------------------------------------------------------------------------
    //
    // Computes pitch angle from vector.
    //
    // Returns:
    //  Pitch in radians.
    //
    // Notes:
    //  Assumes direction vector; negative atan2.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline radian fvec3_t<V>::Pitch(void) const noexcept
    {
        assert(isFinite());
        return -xmath::Atan2(this->m_Y, xmath::Sqrt(this->m_X * this->m_X + this->m_Z * this->m_Z));
    }

    //------------------------------------------------------------------------------
    // yaw
    //------------------------------------------------------------------------------
    //
    // Computes yaw angle from vector.
    //
    // Returns:
    //  Yaw in radians.
    //
    // Notes:
    //  Assumes direction vector.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline radian fvec3_t<V>::Yaw(void) const noexcept
    {
        assert(isFinite());
        return xmath::Atan2(this->m_X, this->m_Z);
    }

    //------------------------------------------------------------------------------
    // pitchYaw
    //------------------------------------------------------------------------------
    //
    // Computes pitch and yaw from vector.
    //
    // Returns:
    //  Pair of pitch and yaw radians.
    //
    // Notes:
    //  Calls pitch/yaw; for direction vectors.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline std::pair<radian, radian> fvec3_t<V>::PitchYaw(void) const noexcept
    {
        assert(isFinite());
        return { Pitch(), Yaw() };
    }

    //------------------------------------------------------------------------------
    // angleBetween
    //------------------------------------------------------------------------------
    //
    // Computes angle between this and v.
    //
    // Params:
    //  v - Other vector.
    //
    // Returns:
    //  Angle in [0, PI] radians.
    //
    // Notes:
    //  acos(dot / lengths); clamps cos.
    //  Assumes finite/non-zero; asserts validity.
    //
    template <bool V >
    inline radian fvec3_t<V>::AngleBetween(const fvec3_t& v) const noexcept
    {
        assert(isFinite() && v.isFinite());
        float d = Dot(*this, v);
        float c = Cross(*this, v).Length();
        auto e = xmath::Atan2(c, d);  
        return e;
    }

    //------------------------------------------------------------------------------
    // rotationTowards
    //------------------------------------------------------------------------------
    //
    // Computes rotation axis/angle to align with dest.
    //
    // Params:
    //  dest - Target direction.
    //
    // Returns:
    //  Pair of axis (unit) and angle radians.
    //
    // Notes:
    //  Handles parallel/opposite cases; default axis for 180 deg.
    //  Assumes finite/non-zero; asserts validity.
    //
    template <bool V >
    inline std::pair<fvec3_t<V>, radian> fvec3_t<V>::RotationTowards(const fvec3_t& dest) const noexcept
    {
        assert(isFinite() && dest.isFinite());
        float d = Length() * dest.Length();
        if (d < 1e-4f) 
        {
            return { Right(), radian(0.0f) };
        }
        float cos = Dot(*this, dest) / d;
        cos = Range(cos, -1.0f, 1.0f);
        fvec3_t rotAxis;
        radian rotAngle;
        if (cos <= -0.999f) 
        {
            fvec3_t absU = dest.AbsCopy();
            rotAxis = (absU.m_X < absU.m_Y) ? ((absU.m_Y < absU.m_Z) ? Forward() : Up()) : Right();
            rotAxis = Cross(rotAxis).NormalizeSafeCopy();
            rotAngle = radian(3.1415926535f);
        }
        else 
        {
            rotAxis  = Cross(dest).NormalizeSafeCopy();
            rotAngle = xmath::Acos(cos);
        }
        return { rotAxis, rotAngle };
    }

    //------------------------------------------------------------------------------
    // vectorToLineSegment
    //------------------------------------------------------------------------------
    //
    // Computes vector from this to closest point on line segment [start, end].
    //
    // Params:
    //  start - Line start.
    //  end - Line end.
    //
    // Returns:
    //  Vector to closest point.
    //
    // Notes:
    //  Projects/clamps t; for distance/closest ops.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::VectorToLineSegment(const fvec3_t& start, const fvec3_t& end) const noexcept
    {
        assert(isFinite() && start.isFinite() && end.isFinite());
        fvec3_t diff = *this - start;
        fvec3_t dir = end - start;
        float t = Dot(diff, dir);
        if (t <= 0.0f) return -diff;
        float sqrLen = dir.LengthSq();
        if (t >= sqrLen) return -(diff - dir);
        t /= sqrLen;
        return -(diff - t * dir);
    }

    //------------------------------------------------------------------------------
    // squareDistToLineSeg
    //------------------------------------------------------------------------------
    //
    // Computes squared distance to line segment [start, end].
    //
    // Params:
    //  start - Line start.
    //  end - Line end.
    //
    // Returns:
    //  Squared distance >=0.
    //
    // Notes:
    //  Uses vectorToLineSegment.lengthsq; faster for comparisons.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline float fvec3_t<V>::SquareDistToLineSeg(const fvec3_t& start, const fvec3_t& end) const noexcept
    {
        assert(isFinite() && start.isFinite() && end.isFinite());
        return VectorToLineSegment(start, end).LengthSq();
    }

    //------------------------------------------------------------------------------
    // closestPointInLineSegment
    //------------------------------------------------------------------------------
    //
    // Computes closest point on line segment [start, end].
    //
    // Params:
    //  start - Line start.
    //  end - Line end.
    //
    // Returns:
    //  Closest point on segment.
    //
    // Notes:
    //  this + vectorToLineSegment.
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClosestPointInLineSegment(const fvec3_t& start, const fvec3_t& end) const noexcept
    {
        assert(isFinite() && start.isFinite() && end.isFinite());
        return *this + VectorToLineSegment(start, end);
    }

    //------------------------------------------------------------------------------
    // closestPointToRectangle
    //------------------------------------------------------------------------------
    //
    // Computes closest point and squared distance to rectangle.
    //
    // Params:
    //  p0 - Rectangle origin.
    //  e0 - First edge vector.
    //  e1 - Second edge vector.
    //  outClosestPoint - Output closest point.
    //
    // Returns:
    //  Squared distance to closest point.
    //
    // Notes:
    //  Parametrized rectangle; clamps params [0,1].
    //  Assumes finite; asserts validity.
    //
    template <bool V >
    inline float fvec3_t<V>::ClosestPointToRectangle(const fvec3_t& p0, const fvec3_t& e0, const fvec3_t& e1, fvec3_t& outClosestPoint) const noexcept
    {
        assert(isFinite() && p0.isFinite() && e0.isFinite() && e1.isFinite());
        const fvec3_t kDiff = p0 - *this;
        const float fA00 = e0.LengthSq();
        const float fA11 = e1.LengthSq();
        float fS = -Dot(kDiff, e0);
        float fT = -Dot(kDiff, e1);
        float fSqrDist = kDiff.LengthSq();
        if (fS < 0.0f) fS = 0.0f;
        else if (fS > fA00) { fS = 1.0f; fSqrDist += fA00 + 2.0f * Dot(kDiff, e0); }
        else { fS /= fA00; fSqrDist += Dot(kDiff, e0) * fS; }
        if (fT < 0.0f) fT = 0.0f;
        else if (fT > fA11) { fT = 1.0f; fSqrDist += fA11 + 2.0f * Dot(kDiff, e1); }
        else { fT /= fA11; fSqrDist += Dot(kDiff, e1) * fT; }
        outClosestPoint = p0 + fS * e0 + fT * e1;
        return std::abs(fSqrDist);
    }

    //------------------------------------------------------------------------------
    // rotateX
    //------------------------------------------------------------------------------
    //
    // Rotates vector around X axis in-place.
    //
    // Params:
    //  rx - Rotation angle.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses SinCos; scalar trig (no SIMD rot matrix yet).
    //  Skips if rx=0; assumes finite/rx valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateX(radian rx) noexcept
    {
        assert(isFinite() && isValid(rx.m_Value));
        if (rx.m_Value == 0.0f) return *this;
        float s, c;
        SinCos(rx, s, c);
        const float y_val = this->m_Y;
        const float z_val = this->m_Z;
        this->m_Y = c * y_val - s * z_val;
        this->m_Z = c * z_val + s * y_val;
        return *this;
    }

    //------------------------------------------------------------------------------
    // rotateY
    //------------------------------------------------------------------------------
    //
    // Rotates vector around Y axis in-place.
    //
    // Params:
    //  ry - Rotation angle.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses SinCos; scalar.
    //  Skips if ry=0; assumes finite/ry valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateY(radian ry) noexcept
    {
        assert(isFinite() && isValid(ry.m_Value));
        if (ry.m_Value == 0.0f) return *this;
        float s, c;
        SinCos(ry, s, c);
        const float x_val = this->m_X;
        const float z_val = this->m_Z;
        this->m_X = c * x_val + s * z_val;
        this->m_Z = c * z_val - s * x_val;
        return *this;
    }

    //------------------------------------------------------------------------------
    // rotateZ
    //------------------------------------------------------------------------------
    //
    // Rotates vector around Z axis in-place.
    //
    // Params:
    //  rz - Rotation angle.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses SinCos; scalar.
    //  Skips if rz=0; assumes finite/rz valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateZ(radian rz) noexcept
    {
        assert(isFinite() && isValid(rz.m_Value));
        if (rz.m_Value == 0.0f) return *this;
        float s, c;
        SinCos(rz, s, c);
        const float x_val = this->m_X;
        const float y_val = this->m_Y;
        this->m_X = c * x_val - s * y_val;
        this->m_Y = c * y_val + s * x_val;
        return *this;
    }

    //------------------------------------------------------------------------------
    // rotate
    //------------------------------------------------------------------------------
    //
    // Applies Euler rotation (Z,X,Y order) in-place.
    //
    // Params:
    //  r - Euler angles (roll, pitch, yaw).
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Calls rotateZ/X/Y; assumes finite/r valid; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Rotate(const radian3& r) noexcept
    {
        assert(isFinite());
        RotateZ(r.m_Roll);
        RotateX(r.m_Pitch);
        RotateY(r.m_Yaw);
        return *this;
    }

    //------------------------------------------------------------------------------
    // gridSnap
    //------------------------------------------------------------------------------
    //
    // Snaps components to grid multiples in-place.
    //
    // Params:
    //  gridX - X grid step.
    //  gridY - Y grid step.
    //  gridZ - Z grid step.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses Round; assumes finite/grids >0; asserts.
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::GridSnap(float gridX, float gridY, float gridZ) noexcept
    {
        assert(isFinite() && gridX > 0.0f && gridY > 0.0f && gridZ > 0.0f);
        this->m_X = Round(this->m_X, gridX);
        this->m_Y = Round(this->m_Y, gridY);
        this->m_Z = Round(this->m_Z, gridZ);
        return *this;
    }

    //------------------------------------------------------------------------------
    // isRightHanded
    //------------------------------------------------------------------------------
    //
    // Checks if points form right-handed 2D triangle (z ignored).
    //
    // Params:
    //  p1 - Second point.
    //  p2 - Third point.
    //
    // Returns:
    //  True if right-handed (negative cross z).
    //
    // Notes:
    //  2D projection; assumes finite; asserts validity.
    //
    template <bool V >
    inline bool fvec3_t<V>::isRightHanded(const fvec3_t& p1, const fvec3_t& p2) const noexcept
    {
        assert(isFinite() && p1.isFinite() && p2.isFinite());
        return ((p1.m_X - this->m_X) * (p2.m_Y - this->m_Y) - (p1.m_Y - this->m_Y) * (p2.m_X - this->m_X)) < 0.0f;
    }

    //------------------------------------------------------------------------------
    // Swizzle methods
    //------------------------------------------------------------------------------

    template <bool V> inline float  fvec3_t<V>::x(void) const noexcept { return this->m_X; }
    template <bool V> inline float  fvec3_t<V>::y(void) const noexcept { return this->m_Y; }
    template <bool V> inline float  fvec3_t<V>::z(void) const noexcept { return this->m_Z; }


    #define SWIZZLE2(func, a, b)                                                               \
        template <bool V> inline fvec2 fvec3_t<V>::func() const noexcept {                   \
            if constexpr (V) { enum : int { X = 0, Y = 1, Z = 2, W = 3 };                      \
               return { _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(W, W, b, a)) }; \
            } else return fvec2(this->m_##a, this->m_##b);                                    \
        }

    SWIZZLE2(xx, X, X)
    SWIZZLE2(xy, X, Y)
    SWIZZLE2(xz, X, Z)
    SWIZZLE2(yx, Y, X)
    SWIZZLE2(yy, Y, Y)
    SWIZZLE2(yz, Y, Z)
    SWIZZLE2(zx, Z, X)
    SWIZZLE2(zy, Z, Y)
    SWIZZLE2(zz, Z, Z)

    #undef SWIZZLE2

    #define SWIZZLE3(func, a, b, c)                                                          \
        template <bool V> inline fvec3_t<V> fvec3_t<V>::func() const noexcept {            \
            if constexpr (V) { enum : int{ X = 0, Y = 1, Z = 2, W = 3 };                     \
                return { _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(W,c,b,a)) }; \
            } else return fvec3_t(this->m_##a, this->m_##b, this->m_##c);                   \
        }

    SWIZZLE3(xxx, X,X,X)
    SWIZZLE3(xxy, X,X,Y)
    SWIZZLE3(xxz, X,X,Z)
    SWIZZLE3(xyx, X,Y,X)
    SWIZZLE3(xyy, X,Y,Y)
    SWIZZLE3(xyz, X,Y,Z)
    SWIZZLE3(xzx, X,Z,X)
    SWIZZLE3(xzy, X,Z,Y)
    SWIZZLE3(xzz, X,Z,Z)
    SWIZZLE3(yxx, Y,X,X)
    SWIZZLE3(yxy, Y,X,Y)
    SWIZZLE3(yxz, Y,X,Z)
    SWIZZLE3(yyx, Y,Y,X)
    SWIZZLE3(yyy, Y,Y,Y)
    SWIZZLE3(yyz, Y,Y,Z)
    SWIZZLE3(yzx, Y,Z,X)
    SWIZZLE3(yzy, Y,Z,Y)
    SWIZZLE3(yzz, Y,Z,Z)
    SWIZZLE3(zxx, Z,X,X)
    SWIZZLE3(zxy, Z,X,Y)
    SWIZZLE3(zxz, Z,X,Z)
    SWIZZLE3(zyx, Z,Y,X)
    SWIZZLE3(zyy, Z,Y,Y)
    SWIZZLE3(zyz, Z,Y,Z)
    SWIZZLE3(zzx, Z,Z,X)
    SWIZZLE3(zzy, Z,Z,Y)
    SWIZZLE3(zzz, Z,Z,Z)

    #undef SWIZZLE3


    #define SWIZZLE4(func, a, b, c, d)                                                              \
        template <bool V> inline fvec4 fvec3_t<V>::func() const noexcept {                        \
            if constexpr (V) { enum : int { X = 0, Y = 1, Z = 2, W = 3 };                           \
                    return { _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(d, c, b, a)) }; \
            } else return fvec4(this->m_##a, this->m_##b, this->m_##c, this->m_##d);               \
        }

    SWIZZLE4(xxxx, X, X, X, X)
    SWIZZLE4(xxxy, X, X, X, Y)
    SWIZZLE4(xxxz, X, X, X, Z)
    SWIZZLE4(xxyx, X, X, Y, X)
    SWIZZLE4(xxyy, X, X, Y, Y)
    SWIZZLE4(xxyz, X, X, Y, Z)
    SWIZZLE4(xxzx, X, X, Z, X)
    SWIZZLE4(xxzy, X, X, Z, Y)
    SWIZZLE4(xxzz, X, X, Z, Z)
    SWIZZLE4(xyxx, X, Y, X, X)
    SWIZZLE4(xyxy, X, Y, X, Y)
    SWIZZLE4(xyxz, X, Y, X, Z)
    SWIZZLE4(xyyx, X, Y, Y, X)
    SWIZZLE4(xyyy, X, Y, Y, Y)
    SWIZZLE4(xyyz, X, Y, Y, Z)
    SWIZZLE4(xyzx, X, Y, Z, X)
    SWIZZLE4(xyzy, X, Y, Z, Y)
    SWIZZLE4(xyzz, X, Y, Z, Z)
    SWIZZLE4(xzxx, X, Z, X, X)
    SWIZZLE4(xzxy, X, Z, X, Y)
    SWIZZLE4(xzxz, X, Z, X, Z)
    SWIZZLE4(xzyx, X, Z, Y, X)
    SWIZZLE4(xzyy, X, Z, Y, Y)
    SWIZZLE4(xzyz, X, Z, Y, Z)
    SWIZZLE4(xzzx, X, Z, Z, X)
    SWIZZLE4(xzzy, X, Z, Z, Y)
    SWIZZLE4(xzzz, X, Z, Z, Z)
    SWIZZLE4(yxxx, Y, X, X, X)
    SWIZZLE4(yxxy, Y, X, X, Y)
    SWIZZLE4(yxxz, Y, X, X, Z)
    SWIZZLE4(yxyx, Y, X, Y, X)
    SWIZZLE4(yxyy, Y, X, Y, Y)
    SWIZZLE4(yxyz, Y, X, Y, Z)
    SWIZZLE4(yxzx, Y, X, Z, X)
    SWIZZLE4(yxzy, Y, X, Z, Y)
    SWIZZLE4(yxzz, Y, X, Z, Z)
    SWIZZLE4(yyxx, Y, Y, X, X)
    SWIZZLE4(yyxy, Y, Y, X, Y)
    SWIZZLE4(yyxz, Y, Y, X, Z)
    SWIZZLE4(yyyx, Y, Y, Y, X)
    SWIZZLE4(yyyy, Y, Y, Y, Y)
    SWIZZLE4(yyyz, Y, Y, Y, Z)
    SWIZZLE4(yyzx, Y, Y, Z, X)
    SWIZZLE4(yyzy, Y, Y, Z, Y)
    SWIZZLE4(yyzz, Y, Y, Z, Z)
    SWIZZLE4(yzxx, Y, Z, X, X)
    SWIZZLE4(yzxy, Y, Z, X, Y)
    SWIZZLE4(yzxz, Y, Z, X, Z)
    SWIZZLE4(yzyx, Y, Z, Y, X)
    SWIZZLE4(yzyy, Y, Z, Y, Y)
    SWIZZLE4(yzyz, Y, Z, Y, Z)
    SWIZZLE4(yzzx, Y, Z, Z, X)
    SWIZZLE4(yzzy, Y, Z, Z, Y)
    SWIZZLE4(yzzz, Y, Z, Z, Z)
    SWIZZLE4(zxxx, Z, X, X, X)
    SWIZZLE4(zxxy, Z, X, X, Y)
    SWIZZLE4(zxxz, Z, X, X, Z)
    SWIZZLE4(zxyx, Z, X, Y, X)
    SWIZZLE4(zxyy, Z, X, Y, Y)
    SWIZZLE4(zxyz, Z, X, Y, Z)
    SWIZZLE4(zxzx, Z, X, Z, X)
    SWIZZLE4(zxzy, Z, X, Z, Y)
    SWIZZLE4(zxzz, Z, X, Z, Z)
    SWIZZLE4(zyxx, Z, Y, X, X)
    SWIZZLE4(zyxy, Z, Y, X, Y)
    SWIZZLE4(zyxz, Z, Y, X, Z)
    SWIZZLE4(zyyx, Z, Y, Y, X)
    SWIZZLE4(zyyy, Z, Y, Y, Y)
    SWIZZLE4(zyyz, Z, Y, Y, Z)
    SWIZZLE4(zyzx, Z, Y, Z, X)
    SWIZZLE4(zyzy, Z, Y, Z, Y)
    SWIZZLE4(zyzz, Z, Y, Z, Z)
    SWIZZLE4(zzxx, Z, Z, X, X)
    SWIZZLE4(zzxy, Z, Z, X, Y)
    SWIZZLE4(zzxz, Z, Z, X, Z)
    SWIZZLE4(zzyx, Z, Z, Y, X)
    SWIZZLE4(zzyy, Z, Z, Y, Y)
    SWIZZLE4(zzyz, Z, Z, Y, Z)
    SWIZZLE4(zzzx, Z, Z, Z, X)
    SWIZZLE4(zzzy, Z, Z, Z, Y)
    SWIZZLE4(zzzz, Z, Z, Z, Z)

    #undef SWIZZLE4

    //------------------------------------------------------------------------------
    // Operator overloads
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator+(const fvec3_t& other) const noexcept
    {
        assert(isFinite() && other.isFinite());
        if constexpr (V) 
        {
            return fvec3_t(_mm_add_ps(this->m_XYZW, other.m_XYZW));
        }
        else 
        {
            return fvec3_t(this->m_X + other.m_X, this->m_Y + other.m_Y, this->m_Z + other.m_Z);
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator-(const fvec3_t& other) const noexcept
    {
        assert(isFinite() && other.isFinite());
        if constexpr (V) 
        {
            return fvec3_t(_mm_sub_ps(this->m_XYZW, other.m_XYZW));
        }
        else 
        {
            return fvec3_t(this->m_X - other.m_X, this->m_Y - other.m_Y, this->m_Z - other.m_Z);
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator*(float scalar) const noexcept
    {
        assert(isFinite() && isValid(scalar));
        if constexpr (V) 
        {
            return fvec3_t(_mm_mul_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else 
        {
            return fvec3_t(this->m_X * scalar, this->m_Y * scalar, this->m_Z * scalar);
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator/(float scalar) const noexcept
    {
        assert(isFinite() && isValid(scalar) && scalar != 0.0f);
        if constexpr (V) 
        {
            return fvec3_t(_mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else 
        {
            return fvec3_t(this->m_X / scalar, this->m_Y / scalar, this->m_Z / scalar);
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator+=(const fvec3_t& other) noexcept
    {
        assert(isFinite() && other.isFinite());
        if constexpr (V) 
        {
            this->m_XYZW = _mm_add_ps(this->m_XYZW, other.m_XYZW);
        }
        else 
        {
            this->m_X += other.m_X;
            this->m_Y += other.m_Y;
            this->m_Z += other.m_Z;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator-=(const fvec3_t& other) noexcept
    {
        assert(isFinite() && other.isFinite());
        if constexpr (V) 
        {
            this->m_XYZW = _mm_sub_ps(this->m_XYZW, other.m_XYZW);
        }
        else 
        {
            this->m_X -= other.m_X;
            this->m_Y -= other.m_Y;
            this->m_Z -= other.m_Z;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator*=(float scalar) noexcept
    {
        assert(isFinite() && isValid(scalar));
        if constexpr (V) 
        {
            this->m_XYZW = _mm_mul_ps(this->m_XYZW, _mm_set1_ps(scalar));
        }
        else 
        {
            this->m_X *= scalar;
            this->m_Y *= scalar;
            this->m_Z *= scalar;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator/=(float scalar) noexcept
    {
        assert(isFinite() && isValid(scalar) && scalar != 0.0f);
        if constexpr (V) 
        {
            this->m_XYZW = _mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar));
        }
        else 
        {
            this->m_X /= scalar;
            this->m_Y /= scalar;
            this->m_Z /= scalar;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline bool fvec3_t<V>::operator==(const fvec3_t& other) const noexcept
    {
        assert(isFinite() && other.isFinite());
        if constexpr (V) 
        {
            return(_mm_movemask_ps(_mm_cmpeq_ps(this->m_XYZW, other.m_XYZW)) & 0x7) == 0x7; // Ignore m_W
        }
        else 
        {
            return this->m_X == other.m_X && this->m_Y == other.m_Y && this->m_Z == other.m_Z;
        }
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline bool fvec3_t<V>::operator!=(const fvec3_t& other) const noexcept
    {
        return !(*this == other);
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline float fvec3_t<V>::operator[](std::int32_t index) const noexcept
    {
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline float& fvec3_t<V>::operator[](std::int32_t index) noexcept
    {
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // Friend operators
    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> operator*(float scalar, const fvec3_t<V>& v) noexcept
    {
        assert(isValid(scalar) && v.isFinite());
        return v * scalar;
    }

    //------------------------------------------------------------------------------
    template <bool V >
    inline fvec3_t<V> operator-(const fvec3_t<V>& v) noexcept
    {
        assert(v.isFinite());
        if constexpr (V) 
        {
            return fvec3_t<V>(_mm_sub_ps(_mm_setzero_ps(), v.m_XYZW));
        }
        else 
        {
            return fvec3_t<V>(-v.m_X, -v.m_Y, -v.m_Z);
        }
    }
}

*****/