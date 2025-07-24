#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
    #error "You must include xmath_vector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    inline fvec4::fvec4(float x, float y, float z, float w) noexcept
        : m_XYZW(_mm_set_ps(w, z, y, x))
    {
    }

    inline fvec4::fvec4(float value) noexcept
        : m_XYZW(_mm_set_ps1(value))
    {
    }

    inline fvec4::fvec4(const fvec3& other, float w) noexcept
        : m_XYZW(_mm_insert_ps(other.m_XYZW, _mm_set_ss(w), 0x30)) // SSE4.1, assume available or alternative
    {
        // Alternative without SSE4.1: _mm_shuffle_ps(other.m_XYZW, _mm_set_ps(w, 0.f, 0.f, 0.f), _MM_SHUFFLE(0, 2, 1, 0)) but adjust
        // But for simplicity, set members
        // m_X = other.m_X; m_Y = other.m_Y; m_Z = other.m_Z; m_W = w;
    }

    inline fvec4::fvec4(float x, const fvec3& other) noexcept
        : m_XYZW(_mm_set_ps(other.m_Z, other.m_Y, other.m_X, x))
    {
    }

    constexpr fvec4::fvec4(const floatx4& reg) noexcept
        : m_XYZW(reg)
    {
    }

    inline fvec4::fvec4(const fvec2& xy, const fvec2& zw) noexcept
        : m_X(xy.m_X), m_Y(xy.m_Y), m_Z(zw.m_X), m_W(zw.m_Y)
    {
    }

    inline fvec4::fvec4(std::span<float> Span) noexcept
    {
        assert(Span.size() >= 4);
        m_XYZW = _mm_loadu_ps(Span.data());
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    constexpr fvec4 fvec4::fromZero(void) noexcept
    {
        return fvec4(_mm_setzero_ps());
    }

    constexpr fvec4 fvec4::fromOne(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 1,1,1,1 }} };
    }

    constexpr fvec4 fvec4::fromUnitX(void) noexcept
    {
        return fromRight();
    }

    constexpr fvec4 fvec4::fromUnitY(void) noexcept
    {
        return fromUp();
    }

    constexpr fvec4 fvec4::fromUnitZ(void) noexcept
    {
        return fromForward();
    }

    constexpr fvec4 fvec4::fromUnitW(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,0,0,1 }} };
    }

    //------------------------------------------------------------------------------
    // Up
    //------------------------------------------------------------------------------
    //
    // Returns up vector (0,1,0).
    //
    constexpr fvec4 fvec4::fromUp(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,1,0,0 }} };
    }

    //------------------------------------------------------------------------------
    // Down
    //------------------------------------------------------------------------------
    //
    // Returns down vector (0,-1,0).
    //
    constexpr fvec4 fvec4::fromDown(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,-1,0,0 }} };
    }

    //------------------------------------------------------------------------------
    // Left
    //------------------------------------------------------------------------------
    //
    // Returns left vector (-1,0,0).
    //
    constexpr fvec4 fvec4::fromLeft(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ -1, 0, 0, 0 }} };
    }

    //------------------------------------------------------------------------------
    // Right
    //------------------------------------------------------------------------------
    //
    // Returns right vector (1,0,0).
    //
    constexpr fvec4 fvec4::fromRight(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 1,0,0,0 }} };
    }

    //------------------------------------------------------------------------------
    // Forward
    //------------------------------------------------------------------------------
    //
    // Returns forward vector (0,0,1).
    //
    constexpr fvec4 fvec4::fromForward(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,0,1,0 }} };
    }

    //------------------------------------------------------------------------------
    // Back
    //------------------------------------------------------------------------------
    //
    // Returns back vector (0,0,-1).
    //
    constexpr fvec4 fvec4::fromBack(void) noexcept
    {
        return fvec4{ floatx4{.m128_f32{ 0,0,-1,0 }} };
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
    // Params:
    //  a - First vector.
    //  b - Second vector.
    //
    // Returns:
    //  The dot product as a scalar.
    //
    // Notes:
    //  Uses SIMD multiplication and horizontal add.
    //
    inline float fvec4::Dot(const fvec4& a, const fvec4& b) noexcept
    {
        floatx4 mul = _mm_mul_ps(a.m_XYZW, b.m_XYZW);
        mul = _mm_add_ps(mul, _mm_shuffle_ps(mul, mul, _MM_SHUFFLE(2, 3, 0, 1)));
        mul = _mm_add_ps(mul, _mm_shuffle_ps(mul, mul, _MM_SHUFFLE(0, 1, 2, 3)));
        float result;
        _mm_store_ss(&result, mul);
        return result;
    }

    inline fvec4 fvec4::Min(const fvec4& a, const fvec4& b) noexcept
    {
        return fvec4{ _mm_min_ps(a.m_XYZW, b.m_XYZW) };
    }

    inline fvec4 fvec4::Max(const fvec4& a, const fvec4& b) noexcept
    {
        return fvec4{ _mm_max_ps(a.m_XYZW, b.m_XYZW) };
    }

    inline fvec4 fvec4::Lerp(const fvec4& a, const fvec4& b, float t) noexcept
    {
        floatx4 tt = _mm_set_ps1(t);
        return fvec4{ _mm_add_ps(a.m_XYZW, _mm_mul_ps(tt, _mm_sub_ps(b.m_XYZW, a.m_XYZW))) };
    }

    inline float fvec4::Distance(const fvec4& a, const fvec4& b) noexcept
    {
        return std::sqrt(Dot(a - b, a - b));
    }

    //------------------------------------------------------------------------------
    // Static methods as members
    //------------------------------------------------------------------------------

    inline float fvec4::Dot(const fvec4& a) const noexcept
    {
        return Dot(*this, a);
    }

    inline fvec4 fvec4::Min(const fvec4& a) const noexcept
    {
        return Min(*this, a);
    }

    inline fvec4 fvec4::Max(const fvec4& a) const noexcept
    {
        return Max(*this, a);
    }

    inline fvec4 fvec4::Lerp(const fvec4& a, float t) const noexcept
    {
        return Lerp(*this, a, t);
    }

    inline float fvec4::Distance(const fvec4& a) const noexcept
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
    // Computes the Euclidean length of the vector.
    //
    // Returns:
    //  The length as a scalar.
    //
    // Notes:
    //  Equivalent to sqrt(LengthSq()).
    //
    inline float fvec4::Length(void) const noexcept
    {
        return xmath::Sqrt(Dot(*this, *this));
    }

    //------------------------------------------------------------------------------
    // Length
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean length of the vector.
    //
    // Returns:
    //  The length as a scalar.
    //
    // Notes:
    //  Equivalent to sqrt(LengthSq()).
    //
    inline fvec4 fvec4::LimitLengthCopy(float MaxLength) const noexcept
    {
        auto l = LengthSq();
        if (l <= (MaxLength * MaxLength)) return *this;
        return (*this) * (MaxLength * xmath::InvSqrt(l));
    }

    //------------------------------------------------------------------------------
    // LengthSq
    //------------------------------------------------------------------------------
    //
    // Computes the squared Euclidean length of the vector.
    //
    // Returns:
    //  The squared length as a scalar.
    //
    // Notes:
    //  Faster than Length() for comparisons.
    //
    inline float fvec4::LengthSq(void) const noexcept
    {
        return Dot(*this, *this);
    }

    //------------------------------------------------------------------------------
    // NormalizeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy of the vector.
    //
    // Returns:
    //  A unit vector in the same direction.
    //
    // Notes:
    //  Assumes non-zero length; use NormalizeSafeCopy for zero-check.
    //
    inline fvec4 fvec4::NormalizeCopy(void) const noexcept
    {
        float len = Length();
        assert(len > 0.f);
        return *this / len;
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
    inline fvec4& fvec4::Normalize(void) noexcept
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
    // Returns a normalized copy or zero if length is zero.
    //
    // Returns:
    //  Unit vector or zero.
    //
    // Notes:
    //  Safe for zero-length vectors.
    //
    inline fvec4 fvec4::NormalizeSafeCopy(void) const noexcept
    {
        float len = Length();
        if (len > 0.f) return *this / len;
        return fromZero();
    }

    //------------------------------------------------------------------------------
    // NormalizeSafe
    //------------------------------------------------------------------------------
    //
    // Normalizes in-place or sets to zero if length is zero.
    //
    // Returns:
    //  Reference to self.
    //
    // Notes:
    //  Safe for zero-length vectors.
    //
    inline fvec4& fvec4::NormalizeSafe(void) noexcept
    {
        float len = Length();
        if (len > 0.f) *this /= len;
        else *this = fromZero();
        return *this;
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all components are finite (not NaN or inf).
    //
    // Returns:
    //  True if finite, false otherwise.
    //
    inline bool fvec4::isFinite(void) const noexcept
    {
        return std::isfinite(m_X) && std::isfinite(m_Y) && std::isfinite(m_Z) && std::isfinite(m_W);
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if all components are within [min, max].
    //
    // Params:
    //  min - Lower bound.
    //  max - Upper bound.
    //
    // Returns:
    //  True if all in range.
    //
    inline bool fvec4::isInRange(float min, float max) const noexcept
    {
        return (m_X >= min && m_X <= max) && (m_Y >= min && m_Y <= max) &&
            (m_Z >= min && m_Z <= max) && (m_W >= min && m_W <= max);
    }

    //------------------------------------------------------------------------------
    // OneOverCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with reciprocal components.
    //
    // Returns:
    //  Vector with 1/x, 1/y, 1/z, 1/w.
    //
    // Notes:
    //  Assumes non-zero components.
    //
    inline fvec4 fvec4::OneOverCopy(void) const noexcept
    {
        return { 1.f / m_X, 1.f / m_Y, 1.f / m_Z, 1.f / m_W };
    }

    //------------------------------------------------------------------------------
    // OneOver
    //------------------------------------------------------------------------------
    //
    // Sets components to their reciprocals in-place.
    //
    // Returns:
    //  Reference to self.
    //
    // Notes:
    //  Assumes non-zero components.
    //
    inline fvec4& fvec4::OneOver(void) noexcept
    {
        m_X = 1.f / m_X;
        m_Y = 1.f / m_Y;
        m_Z = 1.f / m_Z;
        m_W = 1.f / m_W;
        return *this;
    }

    //------------------------------------------------------------------------------
    // AbsCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with absolute value components.
    //
    // Returns:
    //  Vector with |x|, |y|, |z|, |w|.
    //
    inline fvec4 fvec4::AbsCopy(void) const noexcept
    {
        return { std::abs(m_X), std::abs(m_Y), std::abs(m_Z), std::abs(m_W) };
    }

    //------------------------------------------------------------------------------
    // Abs
    //------------------------------------------------------------------------------
    //
    // Sets components to their absolute values in-place.
    //
    // Returns:
    //  Reference to self.
    //
    inline fvec4& fvec4::Abs(void) noexcept
    {
        m_X = std::abs(m_X);
        m_Y = std::abs(m_Y);
        m_Z = std::abs(m_Z);
        m_W = std::abs(m_W);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Reflection
    //------------------------------------------------------------------------------
    //
    // Computes reflection over a normal.
    //
    // Params:
    //  normal - Unit normal vector.
    //
    // Returns:
    //  Reflected vector.
    //
    // Notes:
    //  Formula: v - 2 * dot(v, n) * n
    //
    inline fvec4 fvec4::Reflection(const fvec4& normal) const noexcept
    {
        float d = 2.f * Dot(*this, normal);
        return *this - d * normal;
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
    // DistanceSquare
    //------------------------------------------------------------------------------
    //
    // Computes squared distance to another vector.
    //
    // Params:
    //  v - Other vector.
    //
    // Returns:
    //  Squared distance.
    //
    inline float fvec4::DistanceSquare(const fvec4& v) const noexcept
    {
        return Dot(*this - v, *this - v);
    }

    //------------------------------------------------------------------------------
    // AngleBetween
    //------------------------------------------------------------------------------
    //
    // Computes angle between two vectors.
    //
    // Params:
    //  v - Other vector.
    //
    // Returns:
    //  Angle in radians.
    //
    // Notes:
    //  Uses acos(dot / (len1 * len2))
    //
    inline radian fvec4::AngleBetween(const fvec4& v) const noexcept
    {
        float d = Dot(*this, v);
        float lens = Length() * v.Length();
        return radian{ std::acos(d / lens) };
    }

    //------------------------------------------------------------------------------
    // GridSnap
    //------------------------------------------------------------------------------
    //
    // Snaps components to grid multiples in-place.
    //
    // Params:
    //  gridX - X grid step.
    //  gridY - Y grid step.
    //  gridZ - Z grid step.
    //  gridW - W grid step.
    //
    // Returns:
    //  Reference to self (chainable).
    //
    // Notes:
    //  Uses std::round; assumes finite/grids >0; asserts.
    //
    inline fvec4& fvec4::GridSnap(float gridX, float gridY, float gridZ, float gridW) noexcept
    {
        assert(isFinite() && gridX > 0.0f && gridY > 0.0f && gridZ > 0.0f && gridW > 0.0f);
        m_X = std::round(m_X / gridX) * gridX;
        m_Y = std::round(m_Y / gridY) * gridY;
        m_Z = std::round(m_Z / gridZ) * gridZ;
        m_W = std::round(m_W / gridW) * gridW;
        return *this;
    }

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
        return fvec4{ _mm_add_ps(m_XYZW, other.m_XYZW) };
    }

    inline fvec4 fvec4::operator-(const fvec4& other) const noexcept
    {
        return fvec4{ _mm_sub_ps(m_XYZW, other.m_XYZW) };
    }

    inline fvec4 fvec4::operator*(float scalar) const noexcept
    {
        return fvec4{ _mm_mul_ps(m_XYZW, _mm_set_ps1(scalar)) };
    }

    inline fvec4 fvec4::operator/(float scalar) const noexcept
    {
        return fvec4{ _mm_div_ps(m_XYZW, _mm_set_ps1(scalar)) };
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
        m_XYZW = _mm_mul_ps(m_XYZW, _mm_set_ps1(scalar));
        return *this;
    }

    inline fvec4& fvec4::operator/=(float scalar) noexcept
    {
        m_XYZW = _mm_div_ps(m_XYZW, _mm_set_ps1(scalar));
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

    inline fvec4 operator*(float scalar, const fvec4& v) noexcept
    {
        return v * scalar;
    }

    inline fvec4 operator-(const fvec4& v) noexcept
    {
        return fvec4{ _mm_xor_ps(v.m_XYZW, _mm_set1_ps(-0.f)) };
    }

}