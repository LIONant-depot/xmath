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
    // Constructs a vector with specified x, y, z components.
    //
    // Parameters:
    //  x - The x component.
    //  y - The y component.
    //  z - The z component.
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
    // Constructs a vector from another fvec3_t with opposite SIMD setting.
    //
    // Parameters:
    //  other - The source vector.
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
    // Constructs a vector with all components set to the same value.
    //
    // Parameters:
    //  value - The value for x, y, z components.
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
    // Constructs a direction vector from pitch and yaw angles.
    //
    // Parameters:
    //  pitch - The pitch angle in radians.
    //  yaw - The yaw angle in radians.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(radian pitch, radian yaw) noexcept
    {
        float ps, pc, ys, yc;
        xmath::SinCos(pitch, ps, pc);
        xmath::SinCos(yaw, ys, yc);

        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(1.0f, yc * pc, -ps, ys * pc);
        }
        else
        {
            this->m_X = ys * pc;
            this->m_Y = -ps;
            this->m_Z = yc * pc;
        }
    }

    //------------------------------------------------------------------------------
    // Constructs from a SIMD register (SIMD only).
    //
    // Parameters:
    //  reg - The __m128 register.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const floatx4& reg) noexcept requires V
        : details::f_vec3::simd_data{ .m_XYZW = reg }
    {
    }

    //------------------------------------------------------------------------------
    // Constructs from a fvec2 and z component.
    //
    // Parameters:
    //  other - The fvec2 for x and y.
    //  z - The z component.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(const fvec2& other, float z) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(1.0f, z, other.m_Y, other.m_X);
        }
        else
        {
            this->m_X = other.m_X;
            this->m_Y = other.m_Y;
            this->m_Z = z;
        }
    }

    //------------------------------------------------------------------------------
    // Constructs from x and a fvec2 for y and z.
    //
    // Parameters:
    //  x - The x component.
    //  other - The fvec2 for y and z.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(float x, const fvec2& other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(1.0f, other.m_Y, other.m_X, x);
        }
        else
        {
            this->m_X = x;
            this->m_Y = other.m_X;
            this->m_Z = other.m_Y;
        }
    }

    //------------------------------------------------------------------------------
    // Constructs from a span of floats (at least 3 elements).
    //
    // Parameters:
    //  Span - The span containing x, y, z.
    //
    template <bool V >
    inline fvec3_t<V>::fvec3_t(std::span<float> Span) noexcept
    {
        if (Span.size() >= 3)
        {
            if constexpr (V)
            {
                this->m_XYZW = _mm_set_ps(1.0f, Span[2], Span[1], Span[0]);
            }
            else
            {
                this->m_X = Span[0];
                this->m_Y = Span[1];
                this->m_Z = Span[2];
            }
        }
    }

    //------------------------------------------------------------------------------
    // Constructs from an array of doubles, casting to float.
    //
    // Parameters:
    //  Conversion - The array with x, y, z values.
    //
    template <bool V >
    constexpr fvec3_t<V>::fvec3_t(const std::array<double, 3>& Conversion) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(1.0f, static_cast<float>(Conversion[2]), static_cast<float>(Conversion[1]), static_cast<float>(Conversion[0]));
        }
        else
        {
            this->m_X = static_cast<float>(Conversion[0]);
            this->m_Y = static_cast<float>(Conversion[1]);
            this->m_Z = static_cast<float>(Conversion[2]);
        }
    }

    //------------------------------------------------------------------------------
    // Conversion to array of doubles.
    //
    // Returns:
    //  An array with x, y, z as doubles.
    //
    template <bool V >
    constexpr fvec3_t<V>::operator std::array<double, 3>() const noexcept
    {
        return { static_cast<double>(this->m_X), static_cast<double>(this->m_Y), static_cast<double>(this->m_Z) };
    }

    //------------------------------------------------------------------------------
    // Conversion to string representation.
    //
    // Returns:
    //  A string in "(x, y, z)" format.
    //
    template <bool V >
    inline fvec3_t<V>::operator std::string() const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // Returns a string representation of the vector.
    //
    // Returns:
    //  A string in "(x, y, z)" format.
    //
    template <bool V >
    std::string fvec3_t<V>::ToString(void) const noexcept
    {
        return std::format("({}, {}, {})", this->m_X, this->m_Y, this->m_Z);
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
    template <bool V >
    inline std::ostream& operator<< (std::ostream& os, const fvec3_t<V>& vec) noexcept
    {
        return os << '(' << vec.m_X << ", " << vec.m_Y << ", " << vec.m_Z << ')';
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
    // Returns a random unit vector.
    //
    // Returns:
    //  A normalized random vector.
    //
    // Notes:
    //  Uses random angles; scalar operation.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::fromRandomUnitVector(void) noexcept
    {
        radian theta{ static_cast<float>(rand()) / RAND_MAX * 2.f * xmath::pi_v.m_Value };
        radian phi  { static_cast<float>(rand()) / RAND_MAX * xmath::pi_v.m_Value };
        float sin_phi = xmath::Sin(phi);
        return fvec3_t<V>(sin_phi * xmath::Cos(theta), sin_phi * xmath::Sin(theta), xmath::Cos(phi)).Normalize();
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
    template <bool V >
    inline float fvec3_t<V>::Dot(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            return _mm_cvtss_f32(_mm_dp_ps(a.m_XYZW, b.m_XYZW, 0x71));
        }
        else
        {
            return a.m_X * b.m_X + a.m_Y * b.m_Y + a.m_Z * b.m_Z;
        }
    }

    //------------------------------------------------------------------------------
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
            __m128 tmp0 = _mm_shuffle_ps(a.m_XYZW, a.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1));
            __m128 tmp1 = _mm_shuffle_ps(b.m_XYZW, b.m_XYZW, _MM_SHUFFLE(3, 1, 0, 2));
            __m128 tmp2 = _mm_mul_ps(tmp0, tmp1);
            __m128 tmp3 = _mm_shuffle_ps(a.m_XYZW, a.m_XYZW, _MM_SHUFFLE(3, 1, 0, 2));
            __m128 tmp4 = _mm_shuffle_ps(b.m_XYZW, b.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1));
            __m128 tmp5 = _mm_mul_ps(tmp3, tmp4);
            return fvec3_t<V>(_mm_sub_ps(tmp2, tmp5));
        }
        else
        {
            return fvec3_t<V>(
                a.m_Y * b.m_Z - a.m_Z * b.m_Y,
                a.m_Z * b.m_X - a.m_X * b.m_Z,
                a.m_X * b.m_Y - a.m_Y * b.m_X
            );
        }
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Min(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_min_ps(a.m_XYZW, b.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(
                std::min(a.m_X, b.m_X),
                std::min(a.m_Y, b.m_Y),
                std::min(a.m_Z, b.m_Z)
            );
        }
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Max(const fvec3_t& a, const fvec3_t& b) noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_max_ps(a.m_XYZW, b.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(
                std::max(a.m_X, b.m_X),
                std::max(a.m_Y, b.m_Y),
                std::max(a.m_Z, b.m_Z)
            );
        }
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Lerp(const fvec3_t& a, const fvec3_t& b, float t) noexcept
    {
        if constexpr (V)
        {
            __m128 t_vec = _mm_set1_ps(t);
            return fvec3_t<V>(_mm_add_ps(a.m_XYZW, _mm_mul_ps(_mm_sub_ps(b.m_XYZW, a.m_XYZW), t_vec)));
        }
        else
        {
            return fvec3_t<V>(
                a.m_X + t * (b.m_X - a.m_X),
                a.m_Y + t * (b.m_Y - a.m_Y),
                a.m_Z + t * (b.m_Z - a.m_Z)
            );
        }
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
    template <bool V >
    inline float fvec3_t<V>::Distance(const fvec3_t& a, const fvec3_t& b) noexcept
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
    template <bool V >
    inline float fvec3_t<V>::Dot(const fvec3_t& a) const noexcept
    {
        return Dot(*this, a);
    }

    //------------------------------------------------------------------------------
    // Computes the cross product with another vector.
    //
    // Parameters:
    //  a - The other vector.
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
    // Computes the component-wise minimum with another vector.
    //
    // Parameters:
    //  a - The other vector.
    //
    // Returns:
    //  A vector with the minimum components.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Min(const fvec3_t& a) const noexcept
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Max(const fvec3_t& a) const noexcept
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Lerp(const fvec3_t& a, float t) const noexcept
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
    template <bool V >
    inline float fvec3_t<V>::Distance(const fvec3_t& a) const noexcept
    {
        return Distance(*this, a);
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
    template <bool V >
    inline float fvec3_t<V>::DistanceSquare(const fvec3_t& v) const noexcept
    {
        return (*this - v).LengthSq();
    }

    //------------------------------------------------------------------------------
    // Instance methods - Basic operations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Computes the Euclidean length (magnitude) of the vector.
    //
    // Returns:
    //  The length.
    //
    template <bool V >
    inline float fvec3_t<V>::Length(void) const noexcept
    {
        return xmath::Sqrt(LengthSq());
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
    template <bool V >
    inline float fvec3_t<V>::LengthSq(void) const noexcept
    {
        return Dot(*this);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::NormalizeCopy(void) const noexcept
    {
        float len = Length();
        return len > 0.f ? *this / len : fvec3_t<V>(0.f);
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Normalize(void) noexcept
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::NormalizeSafeCopy(void) const noexcept
    {
        if (!isFinite()) return fvec3_t<V>(0.f);
        float len = Length();
        return len > 0.f ? *this / len : fvec3_t<V>(0,0,1.f);
    }

    //------------------------------------------------------------------------------
    // Safely normalizes the vector in-place, handling infinite/NaN and zero.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::NormalizeSafe(void) noexcept
    {
        if (!isFinite())
        {
            *this = fvec3_t<V>(0.f);
            return *this;
        }
        float len = Length();
        if (len > 0.f) *this /= len;
        else *this = fvec3_t<V>(0.f);
        return *this;
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::LimitLengthCopy(float MaxLength) const noexcept
    {
        float l = LengthSq();
        if (l <= (MaxLength * MaxLength)) return *this;
        return (*this) * (MaxLength * xmath::InvSqrt(l));
    }

    //------------------------------------------------------------------------------
    // Checks if all components are finite (not NaN or infinite).
    //
    // Returns:
    //  True if finite, false otherwise.
    //
    template <bool V >
    inline bool fvec3_t<V>::isFinite(void) const noexcept
    {
        return std::isfinite(this->m_X) && std::isfinite(this->m_Y) && std::isfinite(this->m_Z);
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
    template <bool V >
    inline bool fvec3_t<V>::isInRange(float min, float max) const noexcept
    {
        return this->m_X >= min && this->m_X <= max && this->m_Y >= min && this->m_Y <= max && this->m_Z >= min && this->m_Z <= max;
    }

    //------------------------------------------------------------------------------
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
    // Instance methods - Component-wise math
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Returns a copy with absolute value components.
    //
    // Returns:
    //  The absolute vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AbsCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_andnot_ps(_mm_set1_ps(-0.f), this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Abs(this->m_X), xmath::Abs(this->m_Y), xmath::Abs(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes absolute values in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Abs(void) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_andnot_ps(_mm_set1_ps(-0.f), this->m_XYZW);
        }
        else
        {
            this->m_X = xmath::Abs(this->m_X);
            this->m_Y = xmath::Abs(this->m_Y);
            this->m_Z = xmath::Abs(this->m_Z);
        }
        return *this;
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::OneOverCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_div_ps(_mm_set1_ps(1.f), this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(1.f / this->m_X, 1.f / this->m_Y, 1.f / this->m_Z);
        }
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::OneOver(void) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_div_ps(_mm_set1_ps(1.f), this->m_XYZW);
        }
        else
        {
            this->m_X = 1.f / this->m_X;
            this->m_Y = 1.f / this->m_Y;
            this->m_Z = 1.f / this->m_Z;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise square root of the vector.
    //
    // Returns:
    //  A vector with the square root of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::SqrtCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_sqrt_ps(this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Sqrt(this->m_X), xmath::Sqrt(this->m_Y), xmath::Sqrt(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes the square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Sqrt(void) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_sqrt_ps(this->m_XYZW);
        }
        else
        {
            this->m_X = xmath::Sqrt(this->m_X);
            this->m_Y = xmath::Sqrt(this->m_Y);
            this->m_Z = xmath::Sqrt(this->m_Z);
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise inverse square root (1/sqrt(x)) of the vector.
    //
    // Returns:
    //  A vector with the inverse square root of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::InvSqrtCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_rsqrt_ps(this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::InvSqrt(this->m_X), xmath::InvSqrt(this->m_Y), xmath::InvSqrt(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes the inverse square root of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::InvSqrt(void) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_rsqrt_ps(this->m_XYZW);
        }
        else
        {
            this->m_X = xmath::InvSqrt(this->m_X);
            this->m_Y = xmath::InvSqrt(this->m_Y);
            this->m_Z = xmath::InvSqrt(this->m_Z);
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise sign of the vector.
    //
    // Returns:
    //  A vector with -1 (negative), 0 (zero), or 1 (positive) per component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::SignCopy(void) const noexcept
    {
        if constexpr (V)
        {
            __m128 zero = _mm_setzero_ps();
            __m128 one = _mm_set1_ps(1.f);
            __m128 neg_one = _mm_set1_ps(-1.f);
            __m128 gt_zero = _mm_cmpgt_ps(this->m_XYZW, zero);
            __m128 lt_zero = _mm_cmplt_ps(this->m_XYZW, zero);
            __m128 sign_pos = _mm_and_ps(gt_zero, one);
            __m128 sign_neg = _mm_and_ps(lt_zero, neg_one);
            return fvec3_t<V>(_mm_or_ps(sign_pos, sign_neg));
        }
        else
        {
            return fvec3_t<V>(
                (this->m_X > 0.f) ? 1.f : (this->m_X < 0.f) ? -1.f : 0.f,
                (this->m_Y > 0.f) ? 1.f : (this->m_Y < 0.f) ? -1.f : 0.f,
                (this->m_Z > 0.f) ? 1.f : (this->m_Z < 0.f) ? -1.f : 0.f
            );
        }
    }

    //------------------------------------------------------------------------------
    // Computes the sign of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Sign(void) noexcept
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::FloorCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_floor_ps(this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Floor(this->m_X), xmath::Floor(this->m_Y), xmath::Floor(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes the floor of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Floor(void) noexcept
    {
        *this = FloorCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise ceiling of the vector.
    //
    // Returns:
    //  A vector with the ceiling of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::CeilCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_ceil_ps(this->m_XYZW));
        }
        else
        {
            return fvec3_t<V>(xmath::Ceil(this->m_X), xmath::Ceil(this->m_Y), xmath::Ceil(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes the ceiling of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Ceil(void) noexcept
    {
        *this = CeilCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise fractional part of the vector.
    //
    // Returns:
    //  A vector with the fractional part of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::FractCopy(void) const noexcept
    {
        return *this - FloorCopy();
    }

    //------------------------------------------------------------------------------
    // Computes the fractional part of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Fract(void) noexcept
    {
        *this -= FloorCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise rounding of the vector to the nearest integer.
    //
    // Returns:
    //  A vector with the rounded value of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RoundCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_round_ps(this->m_XYZW, _MM_FROUND_TO_NEAREST_INT));
        }
        else
        {
            return fvec3_t<V>(xmath::Round(this->m_X), xmath::Round(this->m_Y), xmath::Round(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes the rounding of the vector in-place to the nearest integer.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Round(void) noexcept
    {
        *this = RoundCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise truncation of the vector (towards zero).
    //
    // Returns:
    //  A vector with the truncated value of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::TruncCopy(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_round_ps(this->m_XYZW, _MM_FROUND_TO_ZERO));
        }
        else
        {
            return fvec3_t<V>(xmath::Trunc(this->m_X), xmath::Trunc(this->m_Y), xmath::Trunc(this->m_Z));
        }
    }

    //------------------------------------------------------------------------------
    // Computes the truncation of the vector in-place (towards zero).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Trunc(void) noexcept
    {
        *this = TruncCopy();
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ModCopy(float divisor) const noexcept
    {
        return fvec3_t<V>(xmath::FMod(this->m_X, divisor), xmath::FMod(this->m_Y, divisor), xmath::FMod(this->m_Z, divisor));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Mod(float divisor) noexcept
    {
        this->m_X = xmath::FMod(this->m_X, divisor);
        this->m_Y = xmath::FMod(this->m_Y, divisor);
        this->m_Z = xmath::FMod(this->m_Z, divisor);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClampCopy(float min_val, float max_val) const noexcept
    {
        if constexpr (V)
        {
            __m128 min_vec = _mm_set1_ps(min_val);
            __m128 max_vec = _mm_set1_ps(max_val);
            return fvec3_t<V>(_mm_max_ps(min_vec, _mm_min_ps(this->m_XYZW, max_vec)));
        }
        else
        {
            return fvec3_t<V>(xmath::Clamp(this->m_X, min_val, max_val), xmath::Clamp(this->m_Y, min_val, max_val), xmath::Clamp(this->m_Z, min_val, max_val));
        }
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Clamp(float min_val, float max_val) noexcept
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClampCopy(const fvec3_t& min, const fvec3_t& max) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_max_ps(min.m_XYZW, _mm_min_ps(this->m_XYZW, max.m_XYZW)));
        }
        else
        {
            return fvec3_t<V>(xmath::Clamp(this->m_X, min.m_X, max.m_X), xmath::Clamp(this->m_Y, min.m_Y, max.m_Y), xmath::Clamp(this->m_Z, min.m_Z, max.m_Z));
        }
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Clamp(const fvec3_t& min, const fvec3_t& max) noexcept
    {
        *this = ClampCopy(min, max);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Step(float edge) const noexcept
    {
        if constexpr (V)
        {
            __m128 edge_vec = _mm_set1_ps(edge);
            __m128 cmp = _mm_cmpge_ps(this->m_XYZW, edge_vec);
            return fvec3_t<V>(_mm_and_ps(cmp, _mm_set1_ps(1.f)));
        }
        else
        {
            return fvec3_t<V>(
                (this->m_X >= edge) ? 1.f : 0.f,
                (this->m_Y >= edge) ? 1.f : 0.f,
                (this->m_Z >= edge) ? 1.f : 0.f
            );
        }
    }

    //------------------------------------------------------------------------------
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
        fvec3_t t = ((*this - edge0) / (edge1 - edge0)).ClampCopy(0.f, 1.f);
        return t * t * (3.f - 2.f * t);
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise natural logarithm of the vector.
    //
    // Returns:
    //  A vector with the natural log of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::LogCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Log(this->m_X), xmath::Log(this->m_Y), xmath::Log(this->m_Z));
    }

    //------------------------------------------------------------------------------
    // Computes the natural logarithm of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Log(void) noexcept
    {
        this->m_X = xmath::Log(this->m_X);
        this->m_Y = xmath::Log(this->m_Y);
        this->m_Z = xmath::Log(this->m_Z);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise base-2 logarithm of the vector.
    //
    // Returns:
    //  A vector with the base-2 log of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Log2Copy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Log2(this->m_X), xmath::Log2(this->m_Y), xmath::Log2(this->m_Z));
    }

    //------------------------------------------------------------------------------
    // Computes the base-2 logarithm of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Log2(void) noexcept
    {
        this->m_X = xmath::Log2(this->m_X);
        this->m_Y = xmath::Log2(this->m_Y);
        this->m_Z = xmath::Log2(this->m_Z);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::PowCopy(float exp) const noexcept
    {
        return fvec3_t<V>(xmath::Pow(this->m_X, exp), xmath::Pow(this->m_Y, exp), xmath::Pow(this->m_Z, exp));
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Pow(float exp) noexcept
    {
        this->m_X = xmath::Pow(this->m_X, exp);
        this->m_Y = xmath::Pow(this->m_Y, exp);
        this->m_Z = xmath::Pow(this->m_Z, exp);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise exponential (e^x) of the vector.
    //
    // Returns:
    //  A vector with the exponential of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ExpCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Exp(this->m_X), xmath::Exp(this->m_Y), xmath::Exp(this->m_Z));
    }

    //------------------------------------------------------------------------------
    // Computes the exponential of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Exp(void) noexcept
    {
        this->m_X = xmath::Exp(this->m_X);
        this->m_Y = xmath::Exp(this->m_Y);
        this->m_Z = xmath::Exp(this->m_Z);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise sine of the vector (in radians).
    //
    // Returns:
    //  A vector with the sine of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::SinCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Sin(radian{ this->m_X }), xmath::Sin(radian{ this->m_Y }), xmath::Sin(radian{ this->m_Z }));
    }

    //------------------------------------------------------------------------------
    // Computes the sine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Sin(void) noexcept
    {
        this->m_X = xmath::Sin(radian{ this->m_X });
        this->m_Y = xmath::Sin(radian{ this->m_Y });
        this->m_Z = xmath::Sin(radian{ this->m_Z });
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise cosine of the vector (in radians).
    //
    // Returns:
    //  A vector with the cosine of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::CosCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Cos(radian{ this->m_X }), xmath::Cos(radian{ this->m_Y }), xmath::Cos(radian{ this->m_Z }));
    }

    //------------------------------------------------------------------------------
    // Computes the cosine of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Cos(void) noexcept
    {
        this->m_X = xmath::Cos(radian{ this->m_X });
        this->m_Y = xmath::Cos(radian{ this->m_Y });
        this->m_Z = xmath::Cos(radian{ this->m_Z });
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise tangent of the vector (in radians).
    //
    // Returns:
    //  A vector with the tangent of each component.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::TanCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Tan(radian{ this->m_X }), xmath::Tan(radian{ this->m_Y }), xmath::Tan(radian{ this->m_Z }));
    }

    //------------------------------------------------------------------------------
    // Computes the tangent of the vector in-place (in radians).
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Tan(void) noexcept
    {
        this->m_X = xmath::Tan(radian{ this->m_X });
        this->m_Y = xmath::Tan(radian{ this->m_Y });
        this->m_Z = xmath::Tan(radian{ this->m_Z });
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arcsine of the vector.
    //
    // Returns:
    //  A vector with the arcsine of each component (in radians).
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AsinCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Asin(this->m_X).m_Value, xmath::Asin(this->m_Y).m_Value, xmath::Asin(this->m_Z).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arcsine of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Asin(void) noexcept
    {
        this->m_X = xmath::Asin(this->m_X).m_Value;
        this->m_Y = xmath::Asin(this->m_Y).m_Value;
        this->m_Z = xmath::Asin(this->m_Z).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arccosine of the vector.
    //
    // Returns:
    //  A vector with the arccosine of each component (in radians).
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AcosCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Acos(this->m_X).m_Value, xmath::Acos(this->m_Y).m_Value, xmath::Acos(this->m_Z).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arccosine of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Acos(void) noexcept
    {
        this->m_X = xmath::Acos(this->m_X).m_Value;
        this->m_Y = xmath::Acos(this->m_Y).m_Value;
        this->m_Z = xmath::Acos(this->m_Z).m_Value;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the component-wise arctangent of the vector.
    //
    // Returns:
    //  A vector with the arctangent of each component (in radians).
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::AtanCopy(void) const noexcept
    {
        return fvec3_t<V>(xmath::Atan(this->m_X).m_Value, xmath::Atan(this->m_Y).m_Value, xmath::Atan(this->m_Z).m_Value);
    }

    //------------------------------------------------------------------------------
    // Computes the arctangent of the vector in-place.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Atan(void) noexcept
    {
        this->m_X = xmath::Atan(this->m_X).m_Value;
        this->m_Y = xmath::Atan(this->m_Y).m_Value;
        this->m_Z = xmath::Atan(this->m_Z).m_Value;
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Atan2Copy(const fvec3_t& x) const noexcept
    {
        return fvec3_t<V>(xmath::Atan2(this->m_X, x.m_X).m_Value, xmath::Atan2(this->m_Y, x.m_Y).m_Value, xmath::Atan2(this->m_Z, x.m_Z).m_Value);
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Atan2(const fvec3_t& x) noexcept
    {
        *this = Atan2Copy(x);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Instance methods - Geometry
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Computes the pitch angle of the vector.
    //
    // Returns:
    //  Pitch angle in radians.
    //
    template <bool V >
    inline radian fvec3_t<V>::Pitch(void) const noexcept
    {
        return xmath::Atan2(-this->m_Y, xmath::Sqrt(this->m_X * this->m_X + this->m_Z * this->m_Z));
    }

    //------------------------------------------------------------------------------
    // Computes the yaw angle of the vector.
    //
    // Returns:
    //  Yaw angle in radians.
    //
    template <bool V >
    inline radian fvec3_t<V>::Yaw(void) const noexcept
    {
        return xmath::Atan2(this->m_X, this->m_Z);
    }

    //------------------------------------------------------------------------------
    // Computes both pitch and yaw angles of the vector.
    //
    // Returns:
    //  Pair of pitch and yaw angles in radians.
    //
    template <bool V >
    inline std::pair<radian, radian> fvec3_t<V>::PitchYaw(void) const noexcept
    {
        return { Pitch(), Yaw() };
    }

    //------------------------------------------------------------------------------
    // Computes the rotation axis and angle to align with a target vector.
    //
    // Parameters:
    //  dest - The target vector.
    //
    // Returns:
    //  Pair of rotation axis (unit vector) and angle in radians.
    //
    template <bool V >
    inline std::pair<fvec3_t<V>, radian> fvec3_t<V>::RotationTowards(const fvec3_t& dest) const noexcept
    {
        float cos = xmath::Clamp(Dot(dest) / (Length() * dest.Length()), -1.f, 1.f);
        fvec3_t cross = Cross(dest);
        float sin = cross.Length();
        if (sin < 1e-6f)
            return { fvec3_t<V>(0.f, 0.f, 1.f), radian(0.f) };
        return { cross / sin, xmath::Acos(cos) };
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
    template <bool V >
    inline radian fvec3_t<V>::AngleBetween(const fvec3_t& v) const noexcept
    {
        float len = Length() * v.Length();
        if (len == 0.f) return radian(0.f);
        float cos = Dot(v) / len;
        cos = xmath::Clamp(cos, -1.f, 1.f);
        return xmath::Acos(cos);
    }

    //------------------------------------------------------------------------------
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
    // Computes the vector from this point to the closest point on a line segment.
    //
    // Parameters:
    //  start - The start of the line segment.
    //  end - The end of the line segment.
    //
    // Returns:
    //  Vector from this point to the closest point on the segment.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::VectorToLineSegment(const fvec3_t& start, const fvec3_t& end) const noexcept
    {
        fvec3_t dir = end - start;
        float t = Dot(*this - start, dir) / dir.LengthSq();
        t = xmath::Clamp(t, 0.f, 1.f);
        return start + dir * t - *this;
    }

    //------------------------------------------------------------------------------
    // Computes the squared distance to the closest point on a line segment.
    //
    // Parameters:
    //  start - The start of the line segment.
    //  end - The end of the line segment.
    //
    // Returns:
    //  Squared distance to the closest point.
    //
    template <bool V >
    inline float fvec3_t<V>::SquareDistToLineSeg(const fvec3_t& start, const fvec3_t& end) const noexcept
    {
        return VectorToLineSegment(start, end).LengthSq();
    }

    //------------------------------------------------------------------------------
    // Computes the closest point on a line segment to this vector.
    //
    // Parameters:
    //  start - The start of the line segment.
    //  end - The end of the line segment.
    //
    // Returns:
    //  The closest point on the segment.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ClosestPointInLineSegment(const fvec3_t& start, const fvec3_t& end) const noexcept
    {
        return *this + VectorToLineSegment(start, end);
    }

    //------------------------------------------------------------------------------
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
    template <bool V >
    inline float fvec3_t<V>::ClosestPointToRectangle(const fvec3_t& p0, const fvec3_t& e0, const fvec3_t& e1, fvec3_t& outClosestPoint) const noexcept
    {
        fvec3_t diff = *this - p0;
        float t0 = Dot(diff, e0) / e0.LengthSq();
        float t1 = Dot(diff, e1) / e1.LengthSq();
        t0 = xmath::Clamp(t0, 0.f, 1.f);
        t1 = xmath::Clamp(t1, 0.f, 1.f);
        outClosestPoint = p0 + e0 * t0 + e1 * t1;
        return (*this - outClosestPoint).LengthSq();
    }

    //------------------------------------------------------------------------------
    // Computes the rotation of the vector by an angle around the X axis (in radians).
    //
    // Parameters:
    //  rx - The rotation angle around X.
    //
    // Returns:
    //  A rotated vector.
    //
    template <bool V>
    inline fvec3_t<V> fvec3_t<V>::RotateXCopy(radian rx) const noexcept
    {
        float s, c;
        xmath::SinCos(rx, s, c);
        return fvec3_t<V>(this->m_X, this->m_Y * c + this->m_Z * s, -this->m_Y * s + this->m_Z * c);
    }

    //------------------------------------------------------------------------------
    // Rotates the vector in-place by an angle around the X axis (in radians).
    //
    // Parameters:
    //  rx - The rotation angle around X.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateX(radian rx) noexcept
    {
        *this = RotateXCopy(rx);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the rotation of the vector by an angle around the Y axis (in radians).
    //
    // Parameters:
    //  ry - The rotation angle around Y.
    //
    // Returns:
    //  A rotated vector.
    //
    template <bool V>
    inline fvec3_t<V> fvec3_t<V>::RotateYCopy(radian ry) const noexcept
    {
        float s, c;
        xmath::SinCos(ry, s, c);
        return fvec3_t<V>(c * this->m_X + s * this->m_Z, this->m_Y, -s * this->m_X + c * this->m_Z);
    }

    //------------------------------------------------------------------------------
    // Rotates the vector in-place by an angle around the Y axis (in radians).
    //
    // Parameters:
    //  ry - The rotation angle around Y.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateY(radian ry) noexcept
    {
        *this = RotateYCopy(ry);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the rotation of the vector by an angle around the Z axis (in radians).
    //
    // Parameters:
    //  rz - The rotation angle around Z.
    //
    // Returns:
    //  A rotated vector.
    //
    template <bool V>
    inline fvec3_t<V> fvec3_t<V>::RotateZCopy(radian rz) const noexcept
    {
        float s, c;
        xmath::SinCos(rz, s, c);
        return fvec3_t<V>(c * this->m_X - s * this->m_Y, s * this->m_X + c * this->m_Y, this->m_Z);
    }

    //------------------------------------------------------------------------------
    // Rotates the vector in-place by an angle around the Z axis (in radians).
    //
    // Parameters:
    //  rz - The rotation angle around Z.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::RotateZ(radian rz) noexcept
    {
        *this = RotateZCopy(rz);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes the rotation of the vector by Euler angles (in radians).
    //
    // Parameters:
    //  r - The radian3 containing roll, pitch, yaw.
    //
    // Returns:
    //  A rotated vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::RotateCopy(const radian3& r) const noexcept
    {
        return RotateZCopy(r.m_Roll).RotateXCopy(r.m_Pitch).RotateYCopy(r.m_Yaw);
    }

    //------------------------------------------------------------------------------
    // Rotates the vector in-place by Euler angles (in radians).
    //
    // Parameters:
    //  r - The radian3 containing roll, pitch, yaw.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::Rotate(const radian3& r) noexcept
    {
        *this = RotateCopy(r);
        return *this;
    }

    // Computes the inverse rotation of the vector by Euler angles in-place (in radians).
    // Follows right-hand rule with reverse order for inverse: Y (-yaw), X (-pitch), Z (-roll).
    //
    // Parameters:
    //  r - The radian3 containing pitch, yaw, roll.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V>
    inline fvec3_t<V>& fvec3_t<V>::RotateInverse(const radian3& r) noexcept
    {
        RotateY(-r.m_Yaw).RotateX(-r.m_Pitch).RotateZ(-r.m_Roll);
        return *this;
    }

    // Computes the inverse rotation of the vector by Euler angles (in radians).
    // Follows right-hand rule with reverse order for inverse: Y (-yaw), X (-pitch), Z (-roll).
    //
    // Parameters:
    //  r - The radian3 containing pitch, yaw, roll.
    //
    // Returns:
    //  A rotated vector (inverse of RotateCopy).
    //
    template <bool V>
    inline fvec3_t<V> fvec3_t<V>::RotateInverseCopy(const radian3& r) const noexcept
    {
        return RotateYCopy(-r.m_Yaw).RotateXCopy(-r.m_Pitch).RotateZCopy(-r.m_Roll);
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Reflection(const fvec3_t& normal) const noexcept
    {
        return *this - 2.f * Dot(normal) * normal;
    }

    //------------------------------------------------------------------------------
    // Snaps the vector components to a grid.
    //
    // Parameters:
    //  gridX - Grid size for X.
    //  gridY - Grid size for Y.
    //  gridZ - Grid size for Z.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::GridSnap(float gridX, float gridY, float gridZ) noexcept
    {
        this->m_X = xmath::Round(this->m_X / gridX) * gridX;
        this->m_Y = xmath::Round(this->m_Y / gridY) * gridY;
        this->m_Z = xmath::Round(this->m_Z / gridZ) * gridZ;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Checks if points form a right-handed triangle (z ignored).
    //
    // Parameters:
    //  p1 - Second point.
    //  p2 - Third point.
    //
    // Returns:
    //  True if right-handed.
    //
    template <bool V >
    inline bool fvec3_t<V>::isRightHanded(const fvec3_t& p1, const fvec3_t& p2) const noexcept
    {
        return ((p1.m_X - this->m_X) * (p2.m_Y - this->m_Y) - (p1.m_Y - this->m_Y) * (p2.m_X - this->m_X)) < 0.f;
    }

    //------------------------------------------------------------------------------
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
        *this = ProjectCopy(onto);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes a vector perpendicular to this one in the plane defined by the normal.
    //
    // Parameters:
    //  normal - The normal vector defining the plane.
    //
    // Returns:
    //  A perpendicular vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::Perpendicular(const fvec3_t& normal) const noexcept
    {
        return Cross(normal);
    }

    //------------------------------------------------------------------------------
    // Projects this vector onto a plane defined by its normal.
    //
    // Parameters:
    //  normal - The normal vector of the plane.
    //
    // Returns:
    //  The projected vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::ProjectOntoPlane(const fvec3_t& normal) const noexcept
    {
        return *this - ProjectCopy(normal);
    }

    //------------------------------------------------------------------------------
    // Checks if the vector is nearly zero within a tolerance.
    //
    // Parameters:
    //  tolerance - The tolerance (default 1e-6).
    //
    // Returns:
    //  True if nearly zero.
    //
    template <bool V >
    inline bool fvec3_t<V>::isNearlyZero(float tolerance) const noexcept
    {
        return xmath::Abs(this->m_X) <= tolerance && xmath::Abs(this->m_Y) <= tolerance && xmath::Abs(this->m_Z) <= tolerance;
    }

    //------------------------------------------------------------------------------
    // Checks if the vector is normalized (length ~1) within a tolerance.
    //
    // Parameters:
    //  tolerance - The tolerance (default 1e-6).
    //
    // Returns:
    //  True if normalized.
    //
    template <bool V >
    inline bool fvec3_t<V>::isNormalized(float tolerance) const noexcept
    {
        return xmath::Abs(LengthSq() - 1.f) <= tolerance;
    }

    //------------------------------------------------------------------------------
    // Moves towards a target by max distance delta.
    //
    // Parameters:
    //  target - The target vector.
    //  maxDistanceDelta - Max move distance.
    //
    // Returns:
    //  The moved vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::MoveTowardsCopy(const fvec3_t& target, float maxDistanceDelta) const noexcept
    {
        fvec3_t diff = target - *this;
        float dist = diff.Length();
        if (dist <= maxDistanceDelta || dist == 0.f) return target;
        return *this + diff * (maxDistanceDelta / dist);
    }

    //------------------------------------------------------------------------------
    // Moves in-place towards a target by max distance delta.
    //
    // Parameters:
    //  target - The target vector.
    //  maxDistanceDelta - Max move distance.
    //
    // Returns:
    //  Reference to this vector (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::MoveTowards(const fvec3_t& target, float maxDistanceDelta) noexcept
    {
        *this = MoveTowardsCopy(target, maxDistanceDelta);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Swizzle methods for float (HLSL-style, return copy with swizzled components)
    //------------------------------------------------------------------------------
    template <bool V> inline float  fvec3_t<V>::x(void) const noexcept { return this->m_X; }
    template <bool V> inline float  fvec3_t<V>::y(void) const noexcept { return this->m_Y; }
    template <bool V> inline float  fvec3_t<V>::z(void) const noexcept { return this->m_Z; }


    #define SWIZZLE2(func, a, b)                                                                \
        template <bool V> inline fvec2 fvec3_t<V>::func() const noexcept { using namespace xmath; \
            return fvec2(this->m_##a, this->m_##b);    \
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
            if constexpr (V) return fvec3_t{ _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(simde::W,simde::##c,simde::##b,simde::##a)) }; \
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
            if constexpr (V) return fvec4{ _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(simde::##d, simde::##c, simde::##b, simde::##a)) }; \
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

    //------------------------------------------------------------------------------
    // Adds two vectors.
    //
    // Parameters:
    //  other - The vector to add.
    //
    // Returns:
    //  The sum vector.
    //
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
    // Subtracts two vectors.
    //
    // Parameters:
    //  other - The vector to subtract.
    //
    // Returns:
    //  The difference vector.
    //
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
    // Multiplies by a scalar.
    //
    // Parameters:
    //  scalar - The scalar value.
    //
    // Returns:
    //  The scaled vector.
    //
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
    // Component-wise multiplication by another vector.
    //
    // Parameters:
    //  other - The other vector.
    //
    // Returns:
    //  The product vector.
    //
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator*(const fvec3_t& other) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_mul_ps(this->m_XYZW, other.m_XYZW));
        }
        else
        {
            return fvec3_t<V>(this->m_X * other.m_X, this->m_Y * other.m_Y, this->m_Z * other.m_Z);
        }
    }

    //------------------------------------------------------------------------------
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
    template <bool V >
    inline fvec3_t<V> fvec3_t<V>::operator/(float scalar) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(_mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else
        {
            float inv = 1.f / scalar;
            return *this * inv;
        }
    }

    //------------------------------------------------------------------------------
    // Adds another vector in-place.
    //
    // Parameters:
    //  other - The vector to add.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator+=(const fvec3_t& other) noexcept
    {
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
    // Subtracts another vector in-place.
    //
    // Parameters:
    //  other - The vector to subtract.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator-=(const fvec3_t& other) noexcept
    {
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
    // Multiplies by a scalar in-place.
    //
    // Parameters:
    //  scalar - The scalar value.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator*=(float scalar) noexcept
    {
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
    template <bool V >
    inline fvec3_t<V>& fvec3_t<V>::operator/=(float scalar) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar));
        }
        else
        {
            float inv = 1.f / scalar;
            this->m_X *= inv;
            this->m_Y *= inv;
            this->m_Z *= inv;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Checks for equality with another vector.
    //
    // Parameters:
    //  other - The other vector.
    //
    // Returns:
    //  True if equal, false otherwise.
    //
    template <bool V >
    inline bool fvec3_t<V>::operator==(const fvec3_t& other) const noexcept
    {
        return this->m_X == other.m_X && this->m_Y == other.m_Y && this->m_Z == other.m_Z;
    }

    //------------------------------------------------------------------------------
    // Checks for inequality with another vector.
    //
    // Parameters:
    //  other - The other vector.
    //
    // Returns:
    //  True if not equal, false otherwise.
    //
    template <bool V >
    inline bool fvec3_t<V>::operator!=(const fvec3_t& other) const noexcept
    {
        return !(*this == other);
    }

    //------------------------------------------------------------------------------
    // Accesses a component by index (const).
    //
    // Parameters:
    //  index - The index (0 for x, 1 for y, 2 for z).
    //
    // Returns:
    //  The component value.
    //
    // Notes:
    //  Asserts index in [0,2].
    //
    template <bool V >
    inline float fvec3_t<V>::operator[](std::int32_t index) const noexcept
    {
        assert(index >= 0 && index < 3);
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // Accesses a component by index (mutable).
    //
    // Parameters:
    //  index - The index (0 for x, 1 for y, 2 for z).
    //
    // Returns:
    //  Reference to the component.
    //
    // Notes:
    //  Asserts index in [0,2].
    //
    template <bool V >
    inline float& fvec3_t<V>::operator[](std::int32_t index) noexcept
    {
        assert(index >= 0 && index < 3);
        return this->m_Elements[index];
    }

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
    template <bool V >
    inline fvec3_t<V> operator*(float scalar, const fvec3_t<V>& v) noexcept
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
}
