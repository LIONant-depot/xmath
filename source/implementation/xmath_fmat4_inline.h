#pragma once

namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Constructor for diagonal matrix.
    //
    // Params:
    //  diagonal - Value to set on the diagonal; others zero.
    //
    // Notes:
    //  Asserts finite input.
    //
    template <bool V>
    constexpr fmat4_t<V>::fmat4_t(float diagonal) noexcept
    {
        assert(xmath::isFinite(diagonal));
        if constexpr (V)
        {
            floatx4 diag = _mm_set1_ps(diagonal);
            floatx4 zero = _mm_setzero_ps();
            this->m_Columns[0] = _mm_blend_ps(zero, diag, 0x1); // 1 0 0 0
            this->m_Columns[1] = _mm_blend_ps(zero, diag, 0x2); // 0 1 0 0
            this->m_Columns[2] = _mm_blend_ps(zero, diag, 0x4); // 0 0 1 0
            this->m_Columns[3] = _mm_blend_ps(zero, diag, 0x8); // 0 0 0 1
        }
        else
        {
            this->m_00 = diagonal;  this->m_01 = 0.0f;      this->m_02 = 0.0f;      this->m_03 = 0.0f;
            this->m_10 = 0.0f;      this->m_11 = diagonal;  this->m_12 = 0.0f;      this->m_13 = 0.0f;
            this->m_20 = 0.0f;      this->m_21 = 0.0f;      this->m_22 = diagonal;  this->m_23 = 0.0f;
            this->m_30 = 0.0f;      this->m_31 = 0.0f;      this->m_32 = 0.0f;      this->m_33 = diagonal;
        }
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Constructor from array.
    //
    // Params:
    //  arr - Array of 16 floats in column-major order.
    //
    // Notes:
    //  No asserts; assumes valid input.
    //
    template <bool V>
    constexpr fmat4_t<V>::fmat4_t(const std::array<float, 16>& arr) noexcept
    {
        if constexpr (V)
        {
            this->m_Columns[0] = _mm_loadu_ps(&arr[0]);
            this->m_Columns[1] = _mm_loadu_ps(&arr[4]);
            this->m_Columns[2] = _mm_loadu_ps(&arr[8]);
            this->m_Columns[3] = _mm_loadu_ps(&arr[12]);
        }
        else
        {
            this->m_Elements = arr;
        }
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Constructor from span.
    //
    // Params:
    //  span - Span of 16 floats.
    //
    // Notes:
    //  Converts to array internally.
    //
    template <bool V>
    constexpr fmat4_t<V>::fmat4_t(std::span<const float, 16> Span) noexcept
    {
        for ( auto& E : Span )
        {
            this->m_Elements[ static_cast<int>(&E - Span.data()) ] = E;
        }
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Constructor from quaternion.
    //
    // Params:
    //  q - Input quaternion.
    //
    // Notes:
    //  Asserts finite and normalized quaternion.
    //
    template <bool V>
    inline fmat4_t<V>::fmat4_t(const fquat& q) noexcept
    {
        assert(q.isFinite());
        assert(q.isNormalized());

        float xx = q.m_X * q.m_X;
        float xy = q.m_X * q.m_Y;
        float xz = q.m_X * q.m_Z;
        float xw = q.m_X * q.m_W;
        float yy = q.m_Y * q.m_Y;
        float yz = q.m_Y * q.m_Z;
        float yw = q.m_Y * q.m_W;
        float zz = q.m_Z * q.m_Z;
        float zw = q.m_Z * q.m_W;

        if constexpr (V)
        {
            this->m_Columns[0] = _mm_set_ps(0.0f, 2 * (xz - yw), 2 * (xy + zw), 1 - 2 * (yy + zz));
            this->m_Columns[1] = _mm_set_ps(0.0f, 2 * (yz + xw), 1 - 2 * (xx + zz), 2 * (xy - zw));
            this->m_Columns[2] = _mm_set_ps(0.0f, 1 - 2 * (xx + yy), 2 * (yz - xw), 2 * (xz + yw));
            this->m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            this->m_00 = 1 - 2 * (yy + zz); this->m_01 =     2 * (xy - zw); this->m_02 =     2 * (xz + yw); this->m_03 = 0.0f;
            this->m_10 =     2 * (xy + zw); this->m_11 = 1 - 2 * (xx + zz); this->m_12 =     2 * (yz - xw); this->m_13 = 0.0f;
            this->m_20 =     2 * (xz - yw); this->m_21 =     2 * (yz + xw); this->m_22 = 1 - 2 * (xx + yy); this->m_23 = 0.0f;
            this->m_30 =             0.0f;  this->m_31 =              0.0f; this->m_32 =              0.0f; this->m_33 = 1.0f;
        }
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Constructor from Euler angles.
    //
    // Params:
    //  euler - Euler angles in radians (pitch, yaw, roll).
    //
    // Notes:
    //  Delegates to quaternion constructor.
    //
    template <bool V>
    inline fmat4_t<V>::fmat4_t(const radian3& Euler) noexcept
        : fmat4_t( fromRotation(Euler))
    {
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Constructor from translation, rotation, scale (TRS).
    //
    // Params:
    //  translation - Position vector.
    //  rotation - fromRotation quaternion.
    //  scale - Scale vector.
    //
    // Notes:
    //  Asserts finite inputs, normalized rotation, non-zero scale.
    //  Order: scale * rotation * translation (post-multiply).
    //
    template <bool V>
    inline fmat4_t<V>::fmat4_t(const fvec3& translation, const fquat& rotation, const fvec3& scale) noexcept
    {
        setup(translation, rotation, scale);
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Explicit constructor from other SIMD variant.
    //
    // Params:
    //  other - Matrix with opposite SIMD flag.
    //
    // Notes:
    //  Copies elements directly.
    //
    template <bool V>
    constexpr fmat4_t<V>::fmat4_t(const fmat4_t<!V>& other) noexcept
    {
        this->m_Elements = other.m_Elements;
    }

    //------------------------------------------------------------------------------
    // fromIdentity
    //------------------------------------------------------------------------------
    //
    // Static identity matrix.
    //
    // Returns:
    //  fromIdentity matrix.
    //
    template <bool V>
    constexpr fmat4_t<V> fmat4_t<V>::fromIdentity(void) noexcept
    {
        fmat4_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
            m.m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = 1.0f;  m.m_01 = 0.0f;  m.m_02 = 0.0f;  m.m_03 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = 1.0f;  m.m_12 = 0.0f;  m.m_13 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = 1.0f;  m.m_23 = 0.0f;
            m.m_30 = 0.0f;  m.m_31 = 0.0f;  m.m_32 = 0.0f;  m.m_33 = 1.0f;
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Static zero matrix.
    //
    // Returns:
    //  fromZero matrix.
    //
    template <bool V>
    constexpr fmat4_t<V> fmat4_t<V>::fromZero(void) noexcept
    {
        fmat4_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_setzero_ps();
            m.m_Columns[1] = _mm_setzero_ps();
            m.m_Columns[2] = _mm_setzero_ps();
            m.m_Columns[3] = _mm_setzero_ps();
        }
        else
        {
            std::fill(m.m_Elements.begin(), m.m_Elements.end(), 0.0f);
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // fromTranslation
    //------------------------------------------------------------------------------
    //
    // Creates a translation matrix.
    //
    // Params:
    //  t - Translation vector.
    //
    // Returns:
    //  Translation matrix.
    //
    // Notes:
    //  Asserts finite t.
    //  Sets identity with translation in last column (m_03, m_13, m_23, 1.0).
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromTranslation(const fvec3& t) noexcept
    {
        assert(t.isFinite());
        fmat4_t<V> m;
        if constexpr (V)
        {
           m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
           m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
           m.m_Columns[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
           m.m_Columns[3] = _mm_set_ps(1.0f, t.m_Z, t.m_Y, t.m_X);
        }
        else
        {
            m.m_00 = 1.0f; m.m_01 = 0.0f; m.m_02 = 0.0f; m.m_03 = t.m_X;
            m.m_10 = 0.0f; m.m_11 = 1.0f; m.m_12 = 0.0f; m.m_13 = t.m_Y;
            m.m_20 = 0.0f; m.m_21 = 0.0f; m.m_22 = 1.0f; m.m_23 = t.m_Z;
            m.m_30 = 0.0f; m.m_31 = 0.0f; m.m_32 = 0.0f; m.m_33 = 1.0f;
        }
        return m;
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromRotationX(radian angle) noexcept
    {
        fmat4_t<V> m;
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
            m.m_Columns[1] = _mm_set_ps(0.0f, s, c, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, c, -s, 0.0f);
            m.m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = 1.0f;  m.m_01 = 0.0f;  m.m_02 = 0.0f;  m.m_03 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = c;     m.m_12 = -s;    m.m_13 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = s;     m.m_22 = c;     m.m_23 = 0.0f;
            m.m_30 = 0.0f;  m.m_31 = 0.0f;  m.m_32 = 0.0f;  m.m_33 = 1.0f;
        }
        return m;
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromRotationY(radian angle) noexcept
    {
        fmat4_t<V> m;
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, -s, 0.0f, c);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, c, 0.0f, s);
            m.m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = c;     m.m_01 = 0.0f;  m.m_02 = s;     m.m_03 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = 1.0f;  m.m_12 = 0.0f;  m.m_13 = 0.0f;
            m.m_20 = -s;    m.m_21 = 0.0f;  m.m_22 = c;     m.m_23 = 0.0f;
            m.m_30 = 0.0f;  m.m_31 = 0.0f;  m.m_32 = 0.0f;  m.m_33 = 1.0f;
        }
        return m;
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromRotationZ(radian angle) noexcept
    {
        fmat4_t<V> m;
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, s, c);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, c, -s);
            m.m_Columns[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
            m.m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = c;     m.m_01 = -s;    m.m_02 = 0.0f;  m.m_03 = 0.0f;
            m.m_10 = s;     m.m_11 = c;     m.m_12 = 0.0f;  m.m_13 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = 1.0f;  m.m_23 = 0.0f;
            m.m_30 = 0.0f;  m.m_31 = 0.0f;  m.m_32 = 0.0f;  m.m_33 = 1.0f;
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // fromRotation
    //------------------------------------------------------------------------------
    //
    // Static rotation matrix from quaternion.
    //
    // Params:
    //  q - fromRotation quaternion.
    //
    // Returns:
    //  fromRotation matrix.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromRotation(const fquat& q) noexcept
    {
        assert(q.isFinite());
        assert(q.isNormalized());

        const float xx = q.m_X * q.m_X;
        const float xy = q.m_X * q.m_Y;
        const float xz = q.m_X * q.m_Z;
        const float xw = q.m_X * q.m_W;
        const float yy = q.m_Y * q.m_Y;
        const float yz = q.m_Y * q.m_Z;
        const float yw = q.m_Y * q.m_W;
        const float zz = q.m_Z * q.m_Z;
        const float zw = q.m_Z * q.m_W;

        fmat4_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 2 * (xz - yw), 2 * (xy + zw), 1 - 2 * (yy + zz));
            m.m_Columns[1] = _mm_set_ps(0.0f, 2 * (yz + xw), 1 - 2 * (xx + zz), 2 * (xy - zw));
            m.m_Columns[2] = _mm_set_ps(0.0f, 1 - 2 * (xx + yy), 2 * (yz - xw), 2 * (xz + yw));
            m.m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = 1 - 2 * (yy + zz); m.m_01 =     2 * (xy - zw); m.m_02 =     2 * (xz + yw); m.m_03 = 0.0f;
            m.m_10 =     2 * (xy + zw); m.m_11 = 1 - 2 * (xx + zz); m.m_12 =     2 * (yz - xw); m.m_13 = 0.0f;
            m.m_20 =     2 * (xz - yw); m.m_21 =     2 * (yz + xw); m.m_22 = 1 - 2 * (xx + yy); m.m_23 = 0.0f;
            m.m_30 =             0.0f;  m.m_31 =              0.0f; m.m_32 =              0.0f; m.m_33 = 1.0f;
        }

        return m;
    }



    //------------------------------------------------------------------------------
    // fromRotation
    //------------------------------------------------------------------------------
    //
    // Static rotation matrix from axis-angle.
    //
    // Params:
    //  axis - fromRotation axis (normalized).
    //  angle - fromRotation angle in radians.
    //
    // Returns:
    //  fromRotation matrix.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromRotation(const fvec3& axis, radian angle) noexcept
    {
        return fmat4_t<V>(fquat(axis, angle));
    }

    //------------------------------------------------------------------------------
    // fromRotation
    //------------------------------------------------------------------------------
    //
    // Constructor from Euler angles.
    //
    // Params:
    //  euler - Euler angles in radians (pitch, yaw, roll).
    //
    // Returns:
    //  fromRotation matrix.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromRotation(const radian3& Euler) noexcept
    {
        return fromRotation(fquat(Euler));
    }

    //------------------------------------------------------------------------------
    // fromScale
    //------------------------------------------------------------------------------
    //
    // Creates a scale matrix.
    //
    // Params:
    //  s - Scale vector.
    //
    // Returns:
    //  Scale matrix.
    //
    // Notes:
    //  Asserts finite and non-zero s.
    //  Sets diagonal to scale components, translation to zero.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromScale(const fvec3& s) noexcept
    {
        assert(s.isFinite() && s.m_X != 0.0f && s.m_Y != 0.0f && s.m_Z != 0.0f);
        fmat4_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, s.m_X);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, s.m_Y, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, s.m_Z, 0.0f, 0.0f);
            m.m_Columns[3] = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = s.m_X; m.m_01 = 0.0f;  m.m_02 = 0.0f;  m.m_03 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = s.m_Y; m.m_12 = 0.0f;  m.m_13 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = s.m_Z; m.m_23 = 0.0f;
            m.m_30 = 0.0f;  m.m_31 = 0.0f;  m.m_32 = 0.0f;  m.m_33 = 1.0f;
        }
        return m;
    }
    //------------------------------------------------------------------------------
    // Perspective
    //------------------------------------------------------------------------------
    //
    // Static perspective projection matrix (FOV).
    //
    // Params:
    //  fov - Field of view in radians.
    //  aspect - Aspect ratio (width/height).
    //  near_plane - Near clip plane.
    //  far_plane - Far clip plane.
    //
    // Returns:
    //  Perspective matrix.
    //
    // Notes:
    //  Right-hand, zero-to-one depth.
    //  Asserts valid inputs.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromPerspective(radian fov, float aspect, float near_plane, float far_plane) noexcept
    {
        assert(xmath::isFinite(fov.m_Value) && fov.m_Value > 0.0f && fov.m_Value < 3.14159f);
        assert(xmath::isFinite(aspect) && aspect > 0.0f);
        assert(xmath::isFinite(near_plane) && near_plane > 0.0f);
        assert(xmath::isFinite(far_plane) && far_plane > near_plane);

        float tan_half = xmath::Tan(fov / 2.0f);
        fmat4_t<V> m = fromZero();
        m(0, 0) = 1.0f / (aspect * tan_half);
        m(1, 1) = 1.0f / tan_half;
        m(2, 2) = far_plane / (far_plane - near_plane);
        m(2, 3) = 1.0f;
        m(3, 2) = -(far_plane * near_plane) / (far_plane - near_plane);
        return m;
    }

    //------------------------------------------------------------------------------
    // Perspective
    //------------------------------------------------------------------------------
    //
    // Static perspective projection matrix (frustum).
    //
    // Params:
    //  left - Left plane.
    //  right - Right plane.
    //  bottom - Bottom plane.
    //  top - Top plane.
    //  near_plane - Near clip plane.
    //  far_plane - Far clip plane.
    //
    // Returns:
    //  Perspective matrix.
    //
    // Notes:
    //  Asserts valid inputs.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromPerspective(float left, float right, float bottom, float top, float near_plane, float far_plane) noexcept
    {
        assert(xmath::isFinite(left) && left < right);
        assert(xmath::isFinite(bottom) && bottom < top);
        assert(xmath::isFinite(near_plane) && near_plane > 0.0f);
        assert(xmath::isFinite(far_plane) && far_plane > near_plane);

        fmat4_t<V> m = fromZero();
        m(0, 0) = 2.0f * near_plane / (right - left);
        m(1, 1) = 2.0f * near_plane / (top - bottom);
        m(0, 2) = (right + left) / (right - left);
        m(1, 2) = (top + bottom) / (top - bottom);
        m(2, 2) = far_plane / (far_plane - near_plane);
        m(2, 3) = 1.0f;
        m(3, 2) = -(far_plane * near_plane) / (far_plane - near_plane);
        return m;
    }

    //------------------------------------------------------------------------------
    // Ortho
    //------------------------------------------------------------------------------
    //
    // Static orthographic projection matrix (frustum).
    //
    // Params:
    //  left - Left plane.
    //  right - Right plane.
    //  bottom - Bottom plane.
    //  top - Top plane.
    //  near_plane - Near clip plane.
    //  far_plane - Far clip plane.
    //
    // Returns:
    //  Orthographic matrix.
    //
    // Notes:
    //  Asserts valid inputs.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromOrtho(float left, float right, float bottom, float top, float near_plane, float far_plane) noexcept
    {
        assert(xmath::isFinite(left) && left < right);
        assert(xmath::isFinite(bottom) && bottom < top);
        assert(xmath::isFinite(near_plane) && near_plane < far_plane);

        fmat4_t<V> m = fromIdentity();
        m(0, 0) = 2.0f / (right - left);
        m(1, 1) = 2.0f / (top - bottom);
        m(2, 2) = 1.0f / (far_plane - near_plane);
        m(0, 3) = -(right + left) / (right - left);
        m(1, 3) = -(top + bottom) / (top - bottom);
        m(2, 3) = -near_plane / (far_plane - near_plane);
        return m;
    }

    //------------------------------------------------------------------------------
    // Ortho
    //------------------------------------------------------------------------------
    //
    // Static orthographic projection matrix (width/height).
    //
    // Params:
    //  width - View width.
    //  height - View height.
    //  near_plane - Near clip plane.
    //  far_plane - Far clip plane.
    //
    // Returns:
    //  Orthographic matrix.
    //
    // Notes:
    //  Centers at origin.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromOrtho(float width, float height, float near_plane, float far_plane) noexcept
    {
        return fromOrtho(-width / 2.0f, width / 2.0f, -height / 2.0f, height / 2.0f, near_plane, far_plane);
    }

    //------------------------------------------------------------------------------
    // fromLookAt
    //------------------------------------------------------------------------------
    //
    // Static look-at matrix.
    //
    // Params:
    //  eye - Eye position.
    //  target - Target position.
    //  up - Up vector (normalized).
    //
    // Returns:
    //  View matrix.
    //
    // Notes:
    //  Right-hand system; asserts finite and normalized inputs.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromLookAt(const fvec3& eye, const fvec3& target, const fvec3& up) noexcept
    {
        assert(eye.isFinite());
        assert(target.isFinite());
        assert(up.isFinite() && up.isNormalized());

        const fvec3 z = (target - eye).Normalize();
        const fvec3 x = up.Cross(z).Normalize();
        const fvec3 y = z.Cross(x).Normalize();

        fmat4_t<V> m;
        m(0, 0) = x.m_X; m(0, 1) = x.m_Y; m(0, 2) = x.m_Z; m(0, 3) = -x.Dot(eye);
        m(1, 0) = y.m_X; m(1, 1) = y.m_Y; m(1, 2) = y.m_Z; m(1, 3) = -y.Dot(eye);
        m(2, 0) = z.m_X; m(2, 1) = z.m_Y; m(2, 2) = z.m_Z; m(2, 3) = -z.Dot(eye);
        m(3, 0) = 0.0f;  m(3, 1) = 0.0f;  m(3, 2) = 0.0f;  m(3, 3) = 1.0f;
        return m;
    }

    //------------------------------------------------------------------------------
    // fromBillboard
    //------------------------------------------------------------------------------
    //
    // Static billboard matrix.
    //
    // Params:
    //  from - Position.
    //  to - Target direction.
    //  up - Up vector.
    //
    // Returns:
    //  fromBillboard matrix.
    //
    // Notes:
    //  Aligns to view.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::fromBillboard(const fvec3& from, const fvec3& to, const fvec3& up) noexcept
    {
        const fvec3 dir     = (to - from).Normalize();
        const fvec3 right   = up.Cross(dir).Normalize();
        const fvec3 bill_up = dir.Cross(right);
        fmat4_t<V> m;
        m(0, 0) = right.m_X;    m(0, 1) = right.m_Y;        m(0, 2) = right.m_Z;    m(0, 3) = from.m_X;
        m(1, 0) = bill_up.m_X;  m(1, 1) = bill_up.m_Y;      m(1, 2) = bill_up.m_Z;  m(1, 3) = from.m_Y;
        m(2, 0) = dir.m_X;      m(2, 1) = dir.m_Y;          m(2, 2) = dir.m_Z;      m(2, 3) = from.m_Z;
        m(3, 0) = 0.0f;         m(3, 1) = 0.0f;             m(3, 2) = 0.0f;         m(3, 3) = 1.0f;
        return m;
    }

    //------------------------------------------------------------------------------
    // Setup methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // setupIdentity
    //------------------------------------------------------------------------------
    //
    // Sets to identity matrix.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupIdentity(void) noexcept
    {
        *this = fromIdentity();
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupZero
    //------------------------------------------------------------------------------
    //
    // Sets to zero matrix.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupZero(void) noexcept
    {
        *this = fromZero();
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupTranslation
    //------------------------------------------------------------------------------
    //
    // Sets to translation matrix.
    //
    // Params:
    //  t - Translation vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupTranslation(const fvec3& t) noexcept
    {
        *this = fromTranslation(t);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupRotation
    //------------------------------------------------------------------------------
    //
    // Sets to rotation matrix from quaternion.
    //
    // Params:
    //  q - fromRotation quaternion.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupRotation(const fquat& q) noexcept
    {
        *this = fromRotation(q);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupRotation
    //------------------------------------------------------------------------------
    //
    // Sets to rotation matrix from Euler angles.
    //
    // Params:
    //  euler - Euler angles.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupRotation(const radian3& euler) noexcept
    {
        *this = fmat4_t<V>(euler);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupScale
    //------------------------------------------------------------------------------
    //
    // Sets to scale matrix.
    //
    // Params:
    //  s - Scale vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupScale(const fvec3& s) noexcept
    {
        *this = fromScale(s);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupScale
    //------------------------------------------------------------------------------
    //
    // Sets to uniform scale matrix.
    //
    // Params:
    //  s - Uniform scale.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setupScale(float s) noexcept
    {
        *this = fromScale(fvec3(s));
        return *this;
    }

    //------------------------------------------------------------------------------
    // setup
    //------------------------------------------------------------------------------
    //
    // Sets to translation-rotation-scale (TRS) matrix.
    //
    // Params:
    //  translation - Position vector.
    //  rotation - Rotation quaternion.
    //  scale - Scale vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Asserts finite inputs, normalized rotation, non-zero scale.
    //  Directly sets rotation, scales columns, and sets translation for efficiency.
    //  Order: scale * rotation * translation (post-multiply).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::setup(const fvec3& translation, const fquat& rotation, const fvec3& scale) noexcept
    {
        assert(translation.isFinite());
        assert(rotation.isNormalized());
        assert(scale.isFinite() && scale.m_X != 0.0f && scale.m_Y != 0.0f && scale.m_Z != 0.0f);

        // Set rotation
        setupRotation(rotation);

        // Set translation in last column
        this->m_03 = translation.m_X;
        this->m_13 = translation.m_Y;
        this->m_23 = translation.m_Z;
        this->m_33 = 1.0f;

        // Scale the rotation columns
        if constexpr (V)
        {
            floatx4 scale_wzyx = _mm_set_ps(0.0f, scale.m_Z, scale.m_Y, scale.m_X);
            this->m_Columns[0] = _mm_mul_ps(this->m_Columns[0], _mm_shuffle_ps(scale_wzyx, scale_wzyx, _MM_SHUFFLE(0, 0, 0, 0)));
            this->m_Columns[1] = _mm_mul_ps(this->m_Columns[1], _mm_shuffle_ps(scale_wzyx, scale_wzyx, _MM_SHUFFLE(1, 1, 1, 1)));
            this->m_Columns[2] = _mm_mul_ps(this->m_Columns[2], _mm_shuffle_ps(scale_wzyx, scale_wzyx, _MM_SHUFFLE(2, 2, 2, 2)));
        }
        else
        {
            this->m_00 *= scale.m_X; this->m_10 *= scale.m_X; this->m_20 *= scale.m_X;
            this->m_01 *= scale.m_Y; this->m_11 *= scale.m_Y; this->m_21 *= scale.m_Y;
            this->m_02 *= scale.m_Z; this->m_12 *= scale.m_Z; this->m_22 *= scale.m_Z;
        }

        return *this;
    }


    //------------------------------------------------------------------------------
    // Accessors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Access row as vector.
    //
    // Params:
    //  row - Row index (0-3).
    //
    // Returns:
    //  Row vector.
    //
    // Notes:
    //  Asserts valid row.
    //  Gathers from columns for row-major view.
    //
    template <bool V>
    constexpr fvec4 fmat4_t<V>::operator[](size_t row) const noexcept
    {
        assert(row < 4);
        if constexpr (V)
        {
            alignas(16) float col0[4], col1[4], col2[4], col3[4];
            _mm_store_ps(col0, this->m_Columns[0]);
            _mm_store_ps(col1, this->m_Columns[1]);
            _mm_store_ps(col2, this->m_Columns[2]);
            _mm_store_ps(col3, this->m_Columns[3]);
            return fvec4(col0[row], col1[row], col2[row], col3[row]);
        }
        else
        {
            return fvec4(this->m_Cells[0][row], this->m_Cells[1][row], this->m_Cells[2][row], this->m_Cells[3][row]);
        }
    }

    //------------------------------------------------------------------------------
    // operator()
    //------------------------------------------------------------------------------
    //
    // Access element (row-major API).
    //
    // Params:
    //  row - Row index (0-3).
    //  col - Column index (0-3).
    //
    // Returns:
    //  Reference to element.
    //
    // Notes:
    //  Asserts valid indices; transposes internally.
    //
    template <bool V>
    constexpr float& fmat4_t<V>::operator()(size_t row, size_t col) noexcept
    {
        assert(row < 4 && col < 4);
        return this->m_Cells[col][row]; // Transposed access
    }

    //------------------------------------------------------------------------------
    // operator()
    //------------------------------------------------------------------------------
    //
    // Const access element (row-major API).
    //
    // Params:
    //  row - Row index (0-3).
    //  col - Column index (0-3).
    //
    // Returns:
    //  Const reference to element.
    //
    // Notes:
    //  Asserts valid indices; transposes internally.
    //
    template <bool V>
    constexpr const float& fmat4_t<V>::operator()(size_t row, size_t col) const noexcept
    {
        assert(row < 4 && col < 4);
        return this->m_Cells[col][row]; // Transposed access
    }

    //------------------------------------------------------------------------------
/*
    template <bool V>
    constexpr fmat4_t<V>::operator std::span<float, 16>() noexcept
    {
        return this->m_Elements;
    }
*/
    //------------------------------------------------------------------------------
    template <bool V>
    constexpr fmat4_t<V>::operator std::span<const float, 16>() const noexcept
    {
        return this->m_Elements;
    }

    //------------------------------------------------------------------------------
    template <bool V>
    inline fmat4_t<V>::operator xmath::fquat() const noexcept
    {
        return ExtractRotation();
    }

    //------------------------------------------------------------------------------
    // Operations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator+
    //------------------------------------------------------------------------------
    //
    // Matrix addition.
    //
    // Params:
    //  other - Matrix to add.
    //
    // Returns:
    //  Sum matrix.
    //
    // Notes:
    //  Asserts finite inputs.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::operator+(const fmat4_t<V>& other) const noexcept
    {
        assert(this->isFinite());
        assert(other.isFinite());
        fmat4_t<V> result;
        if constexpr (V)
        {
            result.m_Columns[0] = _mm_add_ps(this->m_Columns[0], other.m_Columns[0]);
            result.m_Columns[1] = _mm_add_ps(this->m_Columns[1], other.m_Columns[1]);
            result.m_Columns[2] = _mm_add_ps(this->m_Columns[2], other.m_Columns[2]);
            result.m_Columns[3] = _mm_add_ps(this->m_Columns[3], other.m_Columns[3]);
        }
        else
        {
            for (size_t i = 0; i < 16; ++i) result.m_Elements[i] = this->m_Elements[i] + other.m_Elements[i];
        }
        return result;
    }

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Matrix subtraction.
    //
    // Params:
    //  other - Matrix to subtract.
    //
    // Returns:
    //  Difference matrix.
    //
    // Notes:
    //  Asserts finite inputs.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::operator-(const fmat4_t<V>& other) const noexcept
    {
        assert(this->isFinite());
        assert(other.isFinite());
        fmat4_t<V> result;
        if constexpr (V)
        {
            result.m_Columns[0] = _mm_sub_ps(this->m_Columns[0], other.m_Columns[0]);
            result.m_Columns[1] = _mm_sub_ps(this->m_Columns[1], other.m_Columns[1]);
            result.m_Columns[2] = _mm_sub_ps(this->m_Columns[2], other.m_Columns[2]);
            result.m_Columns[3] = _mm_sub_ps(this->m_Columns[3], other.m_Columns[3]);
        }
        else
        {
            for (size_t i = 0; i < 16; ++i) result.m_Elements[i] = this->m_Elements[i] - other.m_Elements[i];
        }
        return result;
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Matrix multiplication (this * other).
    //
    // Params:
    //  other - Matrix to multiply by.
    //
    // Returns:
    //  Product matrix.
    //
    // Notes:
    //  Asserts finite matrices.
    //  Optimized for column-major storage; computes result.column i = this * other.column i.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::operator*(const fmat4_t<V>& other) const noexcept
    {
        assert(this->isFinite());
        assert(other.isFinite());

        fmat4_t<V> result;
        if constexpr (V)
        {
            // Optimized SSE version, inspired by old SSE4 code (broadcast src1 elements * src2 columns)
            // Compute result.column0 = sum this.column k * other( k, 0 )
            floatx4 col0 = other.m_Columns[0]; // other.column0
            floatx4 bc0 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(0, 0, 0, 0)); // other(0,0) broadcast
            floatx4 bc1 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(1, 1, 1, 1)); // other(1,0) broadcast
            floatx4 bc2 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(2, 2, 2, 2)); // other(2,0) broadcast
            floatx4 bc3 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(3, 3, 3, 3)); // other(3,0) broadcast
            result.m_Columns[0] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_add_ps(_mm_mul_ps(this->m_Columns[2], bc2), _mm_mul_ps(this->m_Columns[3], bc3)));

            // Repeat for column1
            const floatx4 col1 = other.m_Columns[1];
            bc0 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(0, 0, 0, 0));
            bc1 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(1, 1, 1, 1));
            bc2 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(2, 2, 2, 2));
            bc3 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(3, 3, 3, 3));
            result.m_Columns[1] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_add_ps(_mm_mul_ps(this->m_Columns[2], bc2), _mm_mul_ps(this->m_Columns[3], bc3)));

            // Column2
            const floatx4 col2 = other.m_Columns[2];
            bc0 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(0, 0, 0, 0));
            bc1 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(1, 1, 1, 1));
            bc2 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(2, 2, 2, 2));
            bc3 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(3, 3, 3, 3));
            result.m_Columns[2] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_add_ps(_mm_mul_ps(this->m_Columns[2], bc2), _mm_mul_ps(this->m_Columns[3], bc3)));

            // Column3
            const floatx4 col3 = other.m_Columns[3];
            bc0 = _mm_shuffle_ps(col3, col3, _MM_SHUFFLE(0, 0, 0, 0));
            bc1 = _mm_shuffle_ps(col3, col3, _MM_SHUFFLE(1, 1, 1, 1));
            bc2 = _mm_shuffle_ps(col3, col3, _MM_SHUFFLE(2, 2, 2, 2));
            bc3 = _mm_shuffle_ps(col3, col3, _MM_SHUFFLE(3, 3, 3, 3));
            result.m_Columns[3] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_add_ps(_mm_mul_ps(this->m_Columns[2], bc2), _mm_mul_ps(this->m_Columns[3], bc3)));
        }
        else
        {
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    result(i, j) = 0.0f;
                    for (int k = 0; k < 4; ++k)
                        result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return result;
    }

    //------------------------------------------------------------------------------
    // operator+=
    //------------------------------------------------------------------------------
    //
    // In-place addition.
    //
    // Params:
    //  other - Matrix to add.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::operator+=(const fmat4_t<V>& other) noexcept
    {
        *this = *this + other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator-=
    //------------------------------------------------------------------------------
    //
    // In-place subtraction.
    //
    // Params:
    //  other - Matrix to subtract.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::operator-=(const fmat4_t<V>& other) noexcept
    {
        *this = *this - other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*=
    //------------------------------------------------------------------------------
    //
    // In-place multiplication.
    //
    // Params:
    //  other - Matrix to multiply.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::operator*=(const fmat4_t<V>& other) noexcept
    {
        *this = *this * other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Matrix-vector multiplication (vec4).
    //
    // Params:
    //  v - Vector to transform.
    //
    // Returns:
    //  Transformed vector.
    //
    // Notes:
    //  Asserts finite inputs.
    //
    template <bool V>
    inline fvec4 fmat4_t<V>::operator*(const fvec4& v) const noexcept
    {
        assert(this->isFinite());
        assert(v.isFinite());

        if constexpr (V)
        {
            floatx4 result = _mm_setzero_ps();
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[0], _mm_shuffle_ps(v.m_XYZW, v.m_XYZW, _MM_SHUFFLE(0, 0, 0, 0))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[1], _mm_shuffle_ps(v.m_XYZW, v.m_XYZW, _MM_SHUFFLE(1, 1, 1, 1))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[2], _mm_shuffle_ps(v.m_XYZW, v.m_XYZW, _MM_SHUFFLE(2, 2, 2, 2))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[3], _mm_shuffle_ps(v.m_XYZW, v.m_XYZW, _MM_SHUFFLE(3, 3, 3, 3))));
            return fvec4(result);
        }
        else
        {
            return fvec4
            ( this->m_00 * v.m_X + this->m_01 * v.m_Y + this->m_02 * v.m_Z + this->m_03 * v.m_W
            , this->m_10 * v.m_X + this->m_11 * v.m_Y + this->m_12 * v.m_Z + this->m_13 * v.m_W
            , this->m_20 * v.m_X + this->m_21 * v.m_Y + this->m_22 * v.m_Z + this->m_23 * v.m_W
            , this->m_30 * v.m_X + this->m_31 * v.m_Y + this->m_32 * v.m_Z + this->m_33 * v.m_W
            );
        }
    }


    //------------------------------------------------------------------------------
    // Equals
    //------------------------------------------------------------------------------
    //
    // Checks if this vector is approximately equal to another within a tolerance.
    //
    // Parameters:
    //  other - The matrix to compare.
    //  tolerance - The epsilon tolerance.
    //
    // Returns:
    //  True if equal within tolerance, false otherwise.
    //
    template <bool V>
    inline bool fmat4_t<V>::Equals(const fmat4_t& other, float tolerance) const noexcept
    {
        for( int i=0; i<16; ++i)
        {
            if ( false == xmath::Abs(this->m_Elements[i] - other.m_Elements[i]) <= tolerance) return false;
        }
        return true;
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Matrix-vector multiplication (vec3, assumes w=1).
    //
    // Params:
    //  v - Vector to transform.
    //
    // Returns:
    //  Transformed vector (xyz).
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::operator*(const fvec3& v) const noexcept
    {
        return (*this * fvec4(v, 1.0f)).xyz();
    }

    //------------------------------------------------------------------------------
    // Math functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Transpose
    //------------------------------------------------------------------------------
    //
    // Computes the transpose of the matrix.
    //
    // Returns:
    //  Transposed matrix.
    //
    // Notes:
    //  Asserts finite matrix.
    //  Optimized for column-major storage with correct unpack pairing.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::Transpose(void) const noexcept
    {
        assert(this->isFinite());
        fmat4_t<V> result;
        if constexpr (V)
        {
            floatx4 tmp0 = _mm_unpacklo_ps(this->m_Columns[0], this->m_Columns[2]);
            floatx4 tmp1 = _mm_unpackhi_ps(this->m_Columns[0], this->m_Columns[2]);
            floatx4 tmp2 = _mm_unpacklo_ps(this->m_Columns[1], this->m_Columns[3]);
            floatx4 tmp3 = _mm_unpackhi_ps(this->m_Columns[1], this->m_Columns[3]);

            result.m_Columns[0] = _mm_unpacklo_ps(tmp0, tmp2);
            result.m_Columns[1] = _mm_unpackhi_ps(tmp0, tmp2);
            result.m_Columns[2] = _mm_unpacklo_ps(tmp1, tmp3);
            result.m_Columns[3] = _mm_unpackhi_ps(tmp1, tmp3);
        }
        else
        {
            for (size_t i = 0; i < 4; ++i)
                for (size_t j = 0; j < 4; ++j)
                    result(i, j) = (*this)(j, i);
        }
        return result;
    }

    //------------------------------------------------------------------------------
    // Determinant
    //------------------------------------------------------------------------------
    //
    // Computes determinant.
    //
    // Returns:
    //  Determinant value.
    //
    // Notes:
    //  Asserts finite matrix.
    //
    template <bool V>
    inline float fmat4_t<V>::Determinant(void) const noexcept
    {
        assert(this->isFinite());

        // Scalar implementation for both (SIMD det is complex, not always faster)
        float det = 0.0f;
        det += this->m_00 * (this->m_11 * (this->m_22 * this->m_33 - this->m_23 * this->m_32) - this->m_12 * (this->m_21 * this->m_33 - this->m_23 * this->m_31) + this->m_13 * (this->m_21 * this->m_32 - this->m_22 * this->m_31));
        det -= this->m_01 * (this->m_10 * (this->m_22 * this->m_33 - this->m_23 * this->m_32) - this->m_12 * (this->m_20 * this->m_33 - this->m_23 * this->m_30) + this->m_13 * (this->m_20 * this->m_32 - this->m_22 * this->m_30));
        det += this->m_02 * (this->m_10 * (this->m_21 * this->m_33 - this->m_23 * this->m_31) - this->m_11 * (this->m_20 * this->m_33 - this->m_23 * this->m_30) + this->m_13 * (this->m_20 * this->m_31 - this->m_21 * this->m_30));
        det -= this->m_03 * (this->m_10 * (this->m_21 * this->m_32 - this->m_22 * this->m_31) - this->m_11 * (this->m_20 * this->m_32 - this->m_22 * this->m_30) + this->m_12 * (this->m_20 * this->m_31 - this->m_21 * this->m_30));
        return det;
    }

    //------------------------------------------------------------------------------
    // Inverse
    //------------------------------------------------------------------------------
    //
    // Computes full inverse of the matrix.
    //
    // Returns:
    //  Inverse matrix.
    //
    // Notes:
    //  Asserts finite and invertible matrix (|det| > 0.00001f).
    //  General-purpose for any 4x4; uses SIMD adjoint for V=true, Gaussian elimination for V=false.
    //  Optimized for performance and stability, ensuring precision for identity matrix.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::Inverse(void) const noexcept
    {
        assert(this->isFinite());
        fmat4_t<V> result;
        if constexpr (V)
        {
            // SIMD branch: Adjoint method for 4x4
            // Compute determinant using Laplace expansion along first row
            float det_val = this->m_00 * (this->m_11 * this->m_22 * this->m_33 + this->m_12 * this->m_23 * this->m_31 + this->m_13 * this->m_21 * this->m_32 -
                this->m_11 * this->m_23 * this->m_32 - this->m_12 * this->m_21 * this->m_33 - this->m_13 * this->m_22 * this->m_31) -
                this->m_01 * (this->m_10 * this->m_22 * this->m_33 + this->m_12 * this->m_23 * this->m_30 + this->m_13 * this->m_20 * this->m_32 -
                    this->m_10 * this->m_23 * this->m_32 - this->m_12 * this->m_20 * this->m_33 - this->m_13 * this->m_22 * this->m_30) +
                this->m_02 * (this->m_10 * this->m_21 * this->m_33 + this->m_11 * this->m_23 * this->m_30 + this->m_13 * this->m_20 * this->m_31 -
                    this->m_10 * this->m_23 * this->m_31 - this->m_11 * this->m_20 * this->m_33 - this->m_13 * this->m_21 * this->m_30) -
                this->m_03 * (this->m_10 * this->m_21 * this->m_32 + this->m_11 * this->m_22 * this->m_30 + this->m_12 * this->m_20 * this->m_31 -
                    this->m_10 * this->m_22 * this->m_31 - this->m_11 * this->m_20 * this->m_32 - this->m_12 * this->m_21 * this->m_30);
            assert(std::abs(det_val) >= 0.00001f);
            float inv_det = 1.0f / det_val;

            // Compute adjoint (transpose of cofactors)
            floatx4 inv_det_v = _mm_set1_ps(inv_det);

            // Cofactors for column 0 (result.m_00, m_10, m_20, m_30)
            result.m_Columns[0] = _mm_set_ps(
                -inv_det * (this->m_10 * this->m_21 * this->m_32 + this->m_11 * this->m_22 * this->m_30 + this->m_12 * this->m_20 * this->m_31 -
                    this->m_10 * this->m_22 * this->m_31 - this->m_11 * this->m_20 * this->m_32 - this->m_12 * this->m_21 * this->m_30), // m_30
                inv_det * (this->m_10 * this->m_21 * this->m_33 + this->m_11 * this->m_23 * this->m_30 + this->m_13 * this->m_20 * this->m_31 -
                    this->m_10 * this->m_23 * this->m_31 - this->m_11 * this->m_20 * this->m_33 - this->m_13 * this->m_21 * this->m_30), // m_20
                -inv_det * (this->m_10 * this->m_22 * this->m_33 + this->m_12 * this->m_23 * this->m_30 + this->m_13 * this->m_20 * this->m_32 -
                    this->m_10 * this->m_23 * this->m_32 - this->m_12 * this->m_20 * this->m_33 - this->m_13 * this->m_22 * this->m_30), // m_10
                inv_det * (this->m_11 * this->m_22 * this->m_33 + this->m_12 * this->m_23 * this->m_31 + this->m_13 * this->m_21 * this->m_32 -
                    this->m_11 * this->m_23 * this->m_32 - this->m_12 * this->m_21 * this->m_33 - this->m_13 * this->m_22 * this->m_31)); // m_00

            // Cofactors for column 1 (result.m_01, m_11, m_21, m_31)
            result.m_Columns[1] = _mm_set_ps(
                inv_det * (this->m_00 * this->m_21 * this->m_32 + this->m_01 * this->m_22 * this->m_30 + this->m_02 * this->m_20 * this->m_31 -
                    this->m_00 * this->m_22 * this->m_31 - this->m_01 * this->m_20 * this->m_32 - this->m_02 * this->m_21 * this->m_30), // m_31
                -inv_det * (this->m_00 * this->m_21 * this->m_33 + this->m_01 * this->m_23 * this->m_30 + this->m_03 * this->m_20 * this->m_31 -
                    this->m_00 * this->m_23 * this->m_31 - this->m_01 * this->m_20 * this->m_33 - this->m_03 * this->m_21 * this->m_30), // m_21
                inv_det * (this->m_00 * this->m_22 * this->m_33 + this->m_02 * this->m_23 * this->m_30 + this->m_03 * this->m_20 * this->m_32 -
                    this->m_00 * this->m_23 * this->m_32 - this->m_02 * this->m_20 * this->m_33 - this->m_03 * this->m_22 * this->m_30), // m_11
                -inv_det * (this->m_01 * this->m_22 * this->m_33 + this->m_02 * this->m_23 * this->m_31 + this->m_03 * this->m_21 * this->m_32 -
                    this->m_01 * this->m_23 * this->m_32 - this->m_02 * this->m_21 * this->m_33 - this->m_03 * this->m_22 * this->m_31)); // m_01

            // Cofactors for column 2 (result.m_02, m_12, m_22, m_32)
            result.m_Columns[2] = _mm_set_ps(
                -inv_det * (this->m_00 * this->m_11 * this->m_32 + this->m_01 * this->m_12 * this->m_30 + this->m_02 * this->m_10 * this->m_31 -
                    this->m_00 * this->m_12 * this->m_31 - this->m_01 * this->m_10 * this->m_32 - this->m_02 * this->m_11 * this->m_30), // m_32
                inv_det * (this->m_00 * this->m_11 * this->m_33 + this->m_01 * this->m_13 * this->m_30 + this->m_03 * this->m_10 * this->m_31 -
                    this->m_00 * this->m_13 * this->m_31 - this->m_01 * this->m_10 * this->m_33 - this->m_03 * this->m_11 * this->m_30), // m_22
                -inv_det * (this->m_00 * this->m_12 * this->m_33 + this->m_02 * this->m_13 * this->m_30 + this->m_03 * this->m_10 * this->m_32 -
                    this->m_00 * this->m_13 * this->m_32 - this->m_02 * this->m_10 * this->m_33 - this->m_03 * this->m_12 * this->m_30), // m_12
                inv_det * (this->m_01 * this->m_12 * this->m_33 + this->m_02 * this->m_13 * this->m_31 + this->m_03 * this->m_11 * this->m_32 -
                    this->m_01 * this->m_13 * this->m_32 - this->m_02 * this->m_11 * this->m_33 - this->m_03 * this->m_12 * this->m_31)); // m_02

            // Cofactors for column 3 (result.m_03, m_13, m_23, m_33)
            result.m_Columns[3] = _mm_set_ps(
                inv_det * (this->m_00 * this->m_11 * this->m_22 + this->m_01 * this->m_12 * this->m_20 + this->m_02 * this->m_10 * this->m_21 -
                    this->m_00 * this->m_12 * this->m_21 - this->m_01 * this->m_10 * this->m_22 - this->m_02 * this->m_11 * this->m_20), // m_33
                -inv_det * (this->m_00 * this->m_11 * this->m_23 + this->m_01 * this->m_13 * this->m_20 + this->m_03 * this->m_10 * this->m_21 -
                    this->m_00 * this->m_13 * this->m_21 - this->m_01 * this->m_10 * this->m_23 - this->m_03 * this->m_11 * this->m_20), // m_23
                inv_det * (this->m_00 * this->m_12 * this->m_23 + this->m_02 * this->m_13 * this->m_20 + this->m_03 * this->m_10 * this->m_22 -
                    this->m_00 * this->m_13 * this->m_22 - this->m_02 * this->m_10 * this->m_23 - this->m_03 * this->m_12 * this->m_20), // m_13
                -inv_det * (this->m_01 * this->m_12 * this->m_23 + this->m_02 * this->m_13 * this->m_21 + this->m_03 * this->m_11 * this->m_22 -
                    this->m_01 * this->m_13 * this->m_22 - this->m_02 * this->m_11 * this->m_23 - this->m_03 * this->m_12 * this->m_21)); // m_03
        }
        else
        {
            // Scalar Gaussian elimination from old code
            float Scratch[4][8];
            float a;
            int i, j, k, jr, Pivot;
            int Row[4];

            // Initialize augmented matrix
            for (j = 0; j < 4; j++)
            {
                for (k = 0; k < 4; k++)
                {
                    Scratch[j][k] = this->m_Cells[j][k];
                    Scratch[j][4 + k] = 0.0f;
                }
                Scratch[j][4 + j] = 1.0f;
                Row[j] = j;
            }

            // Eliminate columns
            for (i = 0; i < 4; i++)
            {
                // Find pivot
                k = i;
                a = xmath::Abs(Scratch[Row[k]][i]);
                for (j = i + 1; j < 4; j++)
                {
                    jr = Row[j];
                    if (a < xmath::Abs(Scratch[jr][i]))
                    {
                        k = j;
                        a = xmath::Abs(Scratch[jr][i]);
                    }
                }

                // Swap pivot row
                Pivot = Row[k];
                Row[k] = Row[i];
                Row[i] = Pivot;

                // Normalize pivot row
                a = Scratch[Pivot][i];
                assert(xmath::Abs(a) >= 0.00001f);
                Scratch[Pivot][i] = 1.0f;
                for (k = i + 1; k < 8; k++)
                    Scratch[Pivot][k] /= a;

                // Eliminate pivot from other rows
                for (j = i + 1; j < 4; j++)
                {
                    jr = Row[j];
                    a = -Scratch[jr][i];
                    if (a == 0.0f) continue;
                    Scratch[jr][i] = 0.0f;
                    for (k = i + 1; k < 8; k++)
                        Scratch[jr][k] += (a * Scratch[Pivot][k]);
                }
            }

            // Back solve
            for (i = 3; i >= 0; i--)
            {
                Pivot = Row[i];
                for (j = i - 1; j >= 0; j--)
                {
                    jr = Row[j];
                    a = Scratch[jr][i];
                    for (k = i; k < 8; k++)
                        Scratch[jr][k] -= (a * Scratch[Pivot][k]);
                }
            }

            // Copy inverse back
            for (j = 0; j < 4; j++)
            {
                jr = Row[j];
                for (k = 0; k < 4; k++)
                {
                    result.m_Cells[j][k] = Scratch[jr][k + 4];
                }
            }
        }

        return result;
    }

    //------------------------------------------------------------------------------
    // InverseSRT
    //------------------------------------------------------------------------------
    //
    // Computes inverse for SRT matrix (scale * rotation * translation).
    //
    // Returns:
    //  Inverse matrix.
    //
    // Notes:
    //  Computes 3x3 inverse directly (adjoint / det), updates translation as -trans * inverse_3x3.
    //  Assumes non-zero determinant; asserts |det| > 0.00001f.
    //  Optimized for SIMD and scalar; matches Inverse() for SRT matrices.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::InverseSRT(void) const noexcept
    {
        assert(this->isFinite());
        fmat4_t<V> inv;
        const float d = (this->m_00 * (this->m_11 * this->m_22 - this->m_12 * this->m_21) -
                         this->m_01 * (this->m_10 * this->m_22 - this->m_12 * this->m_20) +
                         this->m_02 * (this->m_10 * this->m_21 - this->m_11 * this->m_20));
        assert(xmath::Abs(d) >= 0.00001f);

        const float inv_d = 1.0f / d;

        if constexpr (V)
        {
            // Compute 3x3 adjoint (transpose of cofactors)
            floatx4 cofactor0 = _mm_set_ps(0.0f,
                 inv_d * (this->m_10 * this->m_21 - this->m_11 * this->m_20),  // m_20
                -inv_d * (this->m_10 * this->m_22 - this->m_12 * this->m_20),  // m_10
                 inv_d * (this->m_11 * this->m_22 - this->m_12 * this->m_21)); // m_00

            floatx4 cofactor1 = _mm_set_ps(0.0f,
                -inv_d * (this->m_00 * this->m_21 - this->m_01 * this->m_20),  // m_21
                 inv_d * (this->m_00 * this->m_22 - this->m_02 * this->m_20),  // m_11
                -inv_d * (this->m_01 * this->m_22 - this->m_02 * this->m_21)); // m_01

            floatx4 cofactor2 = _mm_set_ps(0.0f,
                 inv_d * (this->m_00 * this->m_11 - this->m_01 * this->m_10),  // m_22
                -inv_d * (this->m_00 * this->m_12 - this->m_02 * this->m_10),  // m_12
                 inv_d * (this->m_01 * this->m_12 - this->m_02 * this->m_11)); // m_02

            // Set 3x3 inverse (transpose of cofactors)
            inv.m_Columns[0] = cofactor0;
            inv.m_Columns[1] = cofactor1;
            inv.m_Columns[2] = cofactor2;

            // Compute translation: -trans * inv_3x3
            floatx4 trans = _mm_set_ps(0.0f, this->m_23, this->m_13, this->m_03);
            floatx4 new_trans = _mm_setzero_ps();
            new_trans = _mm_sub_ps(new_trans, _mm_mul_ps(cofactor0, _mm_shuffle_ps(trans, trans, _MM_SHUFFLE(0, 0, 0, 0))));
            new_trans = _mm_sub_ps(new_trans, _mm_mul_ps(cofactor1, _mm_shuffle_ps(trans, trans, _MM_SHUFFLE(1, 1, 1, 1))));
            new_trans = _mm_sub_ps(new_trans, _mm_mul_ps(cofactor2, _mm_shuffle_ps(trans, trans, _MM_SHUFFLE(2, 2, 2, 2))));
            inv.m_Columns[3] = new_trans;
            inv.m_33 = 1.0f;
        }
        else
        {
            // Compute 3x3 cofactors and translation directly (compact, matching SIMD operations)
            inv.m_00 = inv_d * (this->m_11 * this->m_22 - this->m_12 * this->m_21);
            inv.m_10 = -inv_d * (this->m_10 * this->m_22 - this->m_12 * this->m_20);
            inv.m_20 = inv_d * (this->m_10 * this->m_21 - this->m_11 * this->m_20);
            inv.m_30 = 0.0f;

            inv.m_01 = -inv_d * (this->m_00 * this->m_21 - this->m_01 * this->m_20);
            inv.m_11 = inv_d * (this->m_00 * this->m_22 - this->m_02 * this->m_20);
            inv.m_21 = -inv_d * (this->m_01 * this->m_22 - this->m_02 * this->m_21);
            inv.m_31 = 0.0f;

            inv.m_02 = inv_d * (this->m_01 * this->m_12 - this->m_02 * this->m_11);
            inv.m_12 = -inv_d * (this->m_00 * this->m_12 - this->m_02 * this->m_10);
            inv.m_22 = inv_d * (this->m_00 * this->m_11 - this->m_01 * this->m_10);
            inv.m_32 = 0.0f;

            // Compute translation directly
            inv.m_03 = -(inv.m_00 * this->m_03 + inv.m_01 * this->m_13 + inv.m_02 * this->m_23);
            inv.m_13 = -(inv.m_10 * this->m_03 + inv.m_11 * this->m_13 + inv.m_12 * this->m_23);
            inv.m_23 = -(inv.m_20 * this->m_03 + inv.m_21 * this->m_13 + inv.m_22 * this->m_23);
            inv.m_33 = 1.0f;

            //// Compute cofactors (exact translation from SIMD)
            //float c00 = inv_d * (this->m_11 * this->m_22 - this->m_12 * this->m_21);
            //float c10 = -inv_d * (this->m_10 * this->m_22 - this->m_12 * this->m_20);
            //float c20 = inv_d * (this->m_10 * this->m_21 - this->m_11 * this->m_20);

            //float c01 = -inv_d * (this->m_00 * this->m_21 - this->m_01 * this->m_20);
            //float c11 = inv_d * (this->m_00 * this->m_22 - this->m_02 * this->m_20);
            //float c21 = -inv_d * (this->m_01 * this->m_22 - this->m_02 * this->m_21);

            //float c02 = inv_d * (this->m_01 * this->m_12 - this->m_02 * this->m_11);
            //float c12 = -inv_d * (this->m_00 * this->m_12 - this->m_02 * this->m_10);
            //float c22 = inv_d * (this->m_00 * this->m_11 - this->m_01 * this->m_10);

            //// Assign to matrix (column-major)
            //inv.m_00 = c00; inv.m_10 = c10; inv.m_20 = c20; inv.m_30 = 0.0f;
            //inv.m_01 = c01; inv.m_11 = c11; inv.m_21 = c21; inv.m_31 = 0.0f;
            //inv.m_02 = c02; inv.m_12 = c12; inv.m_22 = c22; inv.m_32 = 0.0f;

            //// Compute translation: -trans * inv_3x3
            //float tx = this->m_03;
            //float ty = this->m_13;
            //float tz = this->m_23;
            //inv.m_03 = -(c00 * tx + c01 * ty + c02 * tz);
            //inv.m_13 = -(c10 * tx + c11 * ty + c12 * tz);
            //inv.m_23 = -(c20 * tx + c21 * ty + c22 * tz);
            //inv.m_33 = 1.0f;
        }
    
        return inv;
    }


    //------------------------------------------------------------------------------
    // InverseRT
    //------------------------------------------------------------------------------
    //
    // Computes inverse for RT matrix (rotation * translation).
    //
    // Returns:
    //  Inverse matrix.
    //
    // Notes:
    //  Asserts finite matrix.
    //  Computes R^-1 (transpose of 3x3) and -trans * R^-1 for translation.
    //  Optimized for SIMD and scalar; matches Inverse() for RT matrices.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::InverseRT(void) const noexcept
    {
        assert(this->isFinite());
        fmat4_t<V> result;
        if constexpr (V)
        {
            // Transpose 3x3 rotation part
            floatx4 tmp0 = _mm_shuffle_ps(this->m_Columns[0], this->m_Columns[1], _MM_SHUFFLE(1, 0, 1, 0)); // [m_10, m_00, m_11, m_01]
            floatx4 tmp1 = _mm_shuffle_ps(this->m_Columns[0], this->m_Columns[1], _MM_SHUFFLE(3, 2, 3, 2)); // [m_30, m_20, m_31, m_21]
            result.m_Columns[0] = _mm_shuffle_ps(tmp0, this->m_Columns[2], _MM_SHUFFLE(3, 0, 2, 0));        // [m_00, m_10, m_20, 0]
            result.m_Columns[1] = _mm_shuffle_ps(tmp0, this->m_Columns[2], _MM_SHUFFLE(3, 1, 3, 1));        // [m_01, m_11, m_21, 0]
            result.m_Columns[2] = _mm_shuffle_ps(tmp1, this->m_Columns[2], _MM_SHUFFLE(3, 2, 2, 0));        // [m_02, m_12, m_22, 0]

            // Compute -trans * R^-1
            floatx4 trans = _mm_set_ps(0.0f, this->m_23, this->m_13, this->m_03);
            floatx4 new_trans = _mm_setzero_ps();
            new_trans = _mm_sub_ps(new_trans, _mm_mul_ps(result.m_Columns[0], _mm_shuffle_ps(trans, trans, _MM_SHUFFLE(0, 0, 0, 0))));
            new_trans = _mm_sub_ps(new_trans, _mm_mul_ps(result.m_Columns[1], _mm_shuffle_ps(trans, trans, _MM_SHUFFLE(1, 1, 1, 1))));
            new_trans = _mm_sub_ps(new_trans, _mm_mul_ps(result.m_Columns[2], _mm_shuffle_ps(trans, trans, _MM_SHUFFLE(2, 2, 2, 2))));
            result.m_Columns[3] = new_trans;
            result.m_33 = 1.0f;
        }
        else
        {
            // Transpose 3x3 rotation part
            result.m_00 = this->m_00;
            result.m_10 = this->m_01;
            result.m_20 = this->m_02;
            result.m_30 = 0.0f;
            result.m_01 = this->m_10;
            result.m_11 = this->m_11;
            result.m_21 = this->m_12;
            result.m_31 = 0.0f;
            result.m_02 = this->m_20;
            result.m_12 = this->m_21;
            result.m_22 = this->m_22;
            result.m_32 = 0.0f;

            // Compute translation: -trans * R^-1
            result.m_03 = -(result.m_00 * this->m_03 + result.m_01 * this->m_13 + result.m_02 * this->m_23);
            result.m_13 = -(result.m_10 * this->m_03 + result.m_11 * this->m_13 + result.m_12 * this->m_23);
            result.m_23 = -(result.m_20 * this->m_03 + result.m_21 * this->m_13 + result.m_22 * this->m_23);
            result.m_33 = 1.0f;
        }

        return result;
    }

    //------------------------------------------------------------------------------
    // Orthogonalize
    //------------------------------------------------------------------------------
    //
    // Orthogonalizes basis vectors.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::Orthogonalize(void) noexcept
    {
        const fvec3 x = fvec3(this->m_00, this->m_10, this->m_20).Normalize();
        const fvec3 y = (fvec3(this->m_01, this->m_11, this->m_21) - x * x.Dot(fvec3(this->m_01, this->m_11, this->m_21))).Normalize();
        const fvec3 z = x.Cross(y);

        this->m_00 = x.m_X; this->m_10 = x.m_Y; this->m_20 = x.m_Z;
        this->m_01 = y.m_X; this->m_11 = y.m_Y; this->m_21 = y.m_Z;
        this->m_02 = z.m_X; this->m_12 = z.m_Y; this->m_22 = z.m_Z;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Geometry helpers
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // ExtractPosition
    //------------------------------------------------------------------------------
    //
    // Extracts translation component.
    //
    // Returns:
    //  Position vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::ExtractPosition(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3_t<V>(this->m_Columns[3].m128_f32[0], this->m_Columns[3].m128_f32[1], this->m_Columns[3].m128_f32[2]);
        }
        else
        {
            return fvec3_t<V>(this->m_03, this->m_13, this->m_23);
        }
    }

    //------------------------------------------------------------------------------
    // ExtractRotation
    //------------------------------------------------------------------------------
    //
    // Extracts rotation as quaternion.
    //
    // Returns:
    //  Rotation quaternion.
    //
    // Notes:
    //  Handles scale by normalizing basis vectors.
    //  Assumes orthogonal basis (TRS matrix); for stability, normalizes columns.
    //  Flips Z column for negative determinant to handle mirroring.
    //  Normalizes quaternion sign (W >= 0) for consistent representation.
    //  If input from quaternion, extracts back the same (up to sign).
    //
//------------------------------------------------------------------------------
// ExtractRotation
//------------------------------------------------------------------------------
//
// Extracts rotation as quaternion.
//
// Returns:
//  Rotation quaternion.
//
// Notes:
//  Handles scale by normalizing basis vectors.
//  Assumes orthogonal basis (TRS matrix); for stability, normalizes columns.
//  No determinant flip; preserves mirroring if present.
//  If input from quaternion, extracts back the same (up to sign).
//
    template <bool V>
    inline fquat fmat4_t<V>::ExtractRotation(void) const noexcept
    {
        fvec3 scale = ExtractScale();

        fmat4_t<V> normalized = *this;
        if constexpr (V)
        {
            floatx4 inv_scale_x = _mm_set1_ps(1.0f / scale.m_X);
            floatx4 inv_scale_y = _mm_set1_ps(1.0f / scale.m_Y);
            floatx4 inv_scale_z = _mm_set1_ps(1.0f / scale.m_Z);

            normalized.m_Columns[0] = _mm_mul_ps(this->m_Columns[0], inv_scale_x);
            normalized.m_Columns[1] = _mm_mul_ps(this->m_Columns[1], inv_scale_y);
            normalized.m_Columns[2] = _mm_mul_ps(this->m_Columns[2], inv_scale_z);
        }
        else
        {
            normalized.m_00 = this->m_00 / scale.m_X;
            normalized.m_10 = this->m_10 / scale.m_X;
            normalized.m_20 = this->m_20 / scale.m_X;

            normalized.m_01 = this->m_01 / scale.m_Y;
            normalized.m_11 = this->m_11 / scale.m_Y;
            normalized.m_21 = this->m_21 / scale.m_Y;

            normalized.m_02 = this->m_02 / scale.m_Z;
            normalized.m_12 = this->m_12 / scale.m_Z;
            normalized.m_22 = this->m_22 / scale.m_Z;
        }

        // Extract quaternion from normalized 3x3
        float trace = normalized.m_00 + normalized.m_11 + normalized.m_22 + 1.0f;
        fquat q;
        if (trace > 1.0f)
        {
            const float s = 0.5f / std::sqrt(trace);
            q = fquat((normalized.m_21 - normalized.m_12) * s, (normalized.m_02 - normalized.m_20) * s, (normalized.m_10 - normalized.m_01) * s, 0.25f / s);
        }
        else
        {
            if (normalized.m_00 > normalized.m_11 && normalized.m_00 > normalized.m_22)
            {
                const float s = 2.0f * std::sqrt(1.0f + normalized.m_00 - normalized.m_11 - normalized.m_22);
                q = fquat(0.25f * s, (normalized.m_01 + normalized.m_10) / s, (normalized.m_02 + normalized.m_20) / s, (normalized.m_12 - normalized.m_21) / s);
            }
            else if (normalized.m_11 > normalized.m_22)
            {
                const float s = 2.0f * std::sqrt(1.0f + normalized.m_11 - normalized.m_00 - normalized.m_22);
                q = fquat((normalized.m_01 + normalized.m_10) / s, 0.25f * s, (normalized.m_12 + normalized.m_21) / s, (normalized.m_02 - normalized.m_20) / s);
            }
            else
            {
                const float s = 2.0f * std::sqrt(1.0f + normalized.m_22 - normalized.m_00 - normalized.m_11);
                q = fquat((normalized.m_02 + normalized.m_20) / s, (normalized.m_12 + normalized.m_21) / s, 0.25f * s, (normalized.m_01 - normalized.m_10) / s);
            }
        }

        // Normalize sign (ensure W >= 0 for consistent representation)
        if (q.m_W < 0.0f)
        {
            q = -q;
        }

        return q;
    }

    //------------------------------------------------------------------------------
    // ExtractScale
    //------------------------------------------------------------------------------
    //
    // Extracts scale components.
    //
    // Returns:
    //  Scale vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::ExtractScale(void) const noexcept
    {
        if constexpr (V)
        {
            return fvec3( reinterpret_cast<const fvec3*>(&this->m_Columns[0])->Length()
                        , reinterpret_cast<const fvec3*>(&this->m_Columns[1])->Length()
                        , reinterpret_cast<const fvec3*>(&this->m_Columns[2])->Length()
                        );
        }
        else
        {
            return fvec3(fvec3(this->m_00, this->m_10, this->m_20).Length(),
                         fvec3(this->m_01, this->m_11, this->m_21).Length(),
                         fvec3(this->m_02, this->m_12, this->m_22).Length());
        }
    }

    //------------------------------------------------------------------------------
    // Forward
    //------------------------------------------------------------------------------
    //
    // Gets forward vector.
    //
    // Returns:
    //  Normalized forward vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::Forward(void) const noexcept
    {
        return fvec3(this->m_02, this->m_12, this->m_22).Normalize();
    }

    //------------------------------------------------------------------------------
    // Back
    //------------------------------------------------------------------------------
    //
    // Gets back vector.
    //
    // Returns:
    //  Normalized back vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::Back(void) const noexcept
    {
        return -Forward();
    }

    //------------------------------------------------------------------------------
    // Up
    //------------------------------------------------------------------------------
    //
    // Gets up vector.
    //
    // Returns:
    //  Normalized up vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::Up(void) const noexcept
    {
        return fvec3(this->m_01, this->m_11, this->m_21).Normalize();
    }

    //------------------------------------------------------------------------------
    // Down
    //------------------------------------------------------------------------------
    //
    // Gets down vector.
    //
    // Returns:
    //  Normalized down vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::Down(void) const noexcept
    {
        return -Up();
    }

    //------------------------------------------------------------------------------
    // Left
    //------------------------------------------------------------------------------
    //
    // Gets left vector.
    //
    // Returns:
    //  Normalized left vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::Left(void) const noexcept
    {
        return -Right();
    }

    //------------------------------------------------------------------------------
    // Right
    //------------------------------------------------------------------------------
    //
    // Gets right vector.
    //
    // Returns:
    //  Normalized right vector.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::Right(void) const noexcept
    {
        return fvec3(this->m_00, this->m_10, this->m_20).Normalize();
    }

    //------------------------------------------------------------------------------
    // RotateVector
    //------------------------------------------------------------------------------
    //
    // Transforms a vector by the 3x3 rotation/scale part of the matrix.
    //
    // Params:
    //  v - Vector to transform.
    //
    // Returns:
    //  Transformed vector (ignores translation).
    //
    // Notes:
    //  Uses only the 3x3 submatrix [m_00,m_01,m_02; m_10,m_11,m_12; m_20,m_21,m_22].
    //  Optimized for SIMD and scalar; matches element-wise scale for diagonal matrices.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::RotateVector(const fvec3& v) const noexcept
    {
        if constexpr (V)
        {
            floatx4 vec    = _mm_set_ps(0.0f, v.m_Z, v.m_Y, v.m_X);
            floatx4 result = _mm_setzero_ps();
            return fvec3{_mm_add_ps(_mm_add_ps(
                _mm_mul_ps(this->m_Columns[0], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0, 0, 0, 0))),
                _mm_mul_ps(this->m_Columns[1], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(1, 1, 1, 1)))),
                _mm_mul_ps(this->m_Columns[2], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(2, 2, 2, 2)))) };
        }
        else
        {
            float x = this->m_00 * v.m_X + this->m_01 * v.m_Y + this->m_02 * v.m_Z;
            float y = this->m_10 * v.m_X + this->m_11 * v.m_Y + this->m_12 * v.m_Z;
            float z = this->m_20 * v.m_X + this->m_21 * v.m_Y + this->m_22 * v.m_Z;
            return fvec3_t<V>(x, y, z);
        }
    }

    //------------------------------------------------------------------------------
    // InvRotateVector
    //------------------------------------------------------------------------------
    //
    // Transforms a vector by the inverse of the 3x3 rotation/scale part of the matrix.
    //
    // Params:
    //  v - Vector to transform.
    //
    // Returns:
    //  Transformed vector (ignores translation).
    //
    // Notes:
    //  Uses inverse of 3x3 submatrix [m_00,m_01,m_02; m_10,m_11,m_12; m_20,m_21,m_22].
    //  Optimized for SIMD and scalar; inverts scale/rotation for diagonal matrices.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::InvRotateVector(const fvec3& v) const noexcept
    {
        if constexpr (V)
        {
            // Compute inverse of 3x3 submatrix (adjoint / det)
            float det = this->m_00 * (this->m_11 * this->m_22 - this->m_12 * this->m_21) -
                this->m_01 * (this->m_10 * this->m_22 - this->m_12 * this->m_20) +
                this->m_02 * (this->m_10 * this->m_21 - this->m_11 * this->m_20);
            assert(std::abs(det) >= 0.00001f);
            float inv_det = 1.0f / det;

            // Adjoint (transpose of cofactors)
            floatx4 inv_row0 = _mm_set_ps(0.0f,
                inv_det * (this->m_10 * this->m_21 - this->m_11 * this->m_20), // m_20
                -inv_det * (this->m_10 * this->m_22 - this->m_12 * this->m_20), // m_10
                inv_det * (this->m_11 * this->m_22 - this->m_12 * this->m_21)); // m_00

            floatx4 inv_row1 = _mm_set_ps(0.0f,
                -inv_det * (this->m_00 * this->m_21 - this->m_01 * this->m_20), // m_21
                inv_det * (this->m_00 * this->m_22 - this->m_02 * this->m_20), // m_11
                -inv_det * (this->m_01 * this->m_22 - this->m_02 * this->m_21)); // m_01

            floatx4 inv_row2 = _mm_set_ps(0.0f,
                inv_det * (this->m_00 * this->m_11 - this->m_01 * this->m_10), // m_22
                -inv_det * (this->m_00 * this->m_12 - this->m_02 * this->m_10), // m_12
                inv_det * (this->m_01 * this->m_12 - this->m_02 * this->m_11)); // m_02

            // Apply inverse 3x3 to vector
            floatx4 vec = _mm_set_ps(0.0f, v.m_Z, v.m_Y, v.m_X);
            floatx4 result = _mm_setzero_ps();
            result = _mm_add_ps(_mm_add_ps(
                _mm_mul_ps(inv_row0, _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0, 0, 0, 0))),
                _mm_mul_ps(inv_row1, _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(1, 1, 1, 1)))),
                _mm_mul_ps(inv_row2, _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(2, 2, 2, 2))));
            return fvec3_t<V>(result.m128_f32[0], result.m128_f32[1], result.m128_f32[2]);
        }
        else
        {
            // Compute inverse of 3x3 submatrix (adjoint / det)
            float det = this->m_00 * (this->m_11 * this->m_22 - this->m_12 * this->m_21) -
                this->m_01 * (this->m_10 * this->m_22 - this->m_12 * this->m_20) +
                this->m_02 * (this->m_10 * this->m_21 - this->m_11 * this->m_20);
            assert(std::abs(det) >= 0.00001f);
            float inv_det = 1.0f / det;

            return fvec3_t<V>(
                inv_det * (this->m_11 * this->m_22 - this->m_12 * this->m_21) * v.m_X +
                -inv_det * (this->m_01 * this->m_22 - this->m_02 * this->m_21) * v.m_Y +
                inv_det * (this->m_01 * this->m_12 - this->m_02 * this->m_11) * v.m_Z,
                -inv_det * (this->m_10 * this->m_22 - this->m_12 * this->m_20) * v.m_X +
                inv_det * (this->m_00 * this->m_22 - this->m_02 * this->m_20) * v.m_Y +
                -inv_det * (this->m_00 * this->m_12 - this->m_02 * this->m_10) * v.m_Z,
                inv_det * (this->m_10 * this->m_21 - this->m_11 * this->m_20) * v.m_X +
                -inv_det * (this->m_00 * this->m_21 - this->m_01 * this->m_20) * v.m_Y +
                inv_det * (this->m_00 * this->m_11 - this->m_01 * this->m_10) * v.m_Z
            );
        }
    }
    //------------------------------------------------------------------------------
    // TransformPosition
    //------------------------------------------------------------------------------
    //
    // Transforms position (full matrix).
    //
    // Params:
    //  p - Position to transform.
    //
    // Returns:
    //  Transformed position.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::TransformPosition(const fvec3& p) const noexcept
    {
        return (*this * p);
    }

    //------------------------------------------------------------------------------
    // TransformDirection
    //------------------------------------------------------------------------------
    //
    // Transforms direction (rotation only).
    //
    // Params:
    //  d - Direction to transform.
    //
    // Returns:
    //  Transformed direction.
    //
    template <bool V>
    inline fvec3 fmat4_t<V>::TransformDirection(const fvec3& d) const noexcept
    {
        return RotateVector(d);
    }

    //------------------------------------------------------------------------------
    // Mutable chaining methods (post-multiply)
    //------------------------------------------------------------------------------

     //------------------------------------------------------------------------------
    // Translate
    //------------------------------------------------------------------------------
    //
    // Post-multiplies translation (local space).
    //
    // Params:
    //  t - Translation vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Asserts finite t.
    //  Adds t transformed by current rotation/scale to translation column.
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::Translate(const fvec3& t) noexcept
    {
        assert(t.isFinite());
        return *this = fmat4_t::fromTranslation(t) * (*this);

        /*
        if constexpr (V)
        {
            // Compute t' = R * S * t (rotation * scale applied to t)
            floatx4 t_vec = _mm_set_ps(0.0f, t.m_Z, t.m_Y, t.m_X);
            floatx4 result = _mm_setzero_ps();
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[0], _mm_shuffle_ps(t_vec, t_vec, _MM_SHUFFLE(0, 0, 0, 0))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[1], _mm_shuffle_ps(t_vec, t_vec, _MM_SHUFFLE(1, 1, 1, 1))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[2], _mm_shuffle_ps(t_vec, t_vec, _MM_SHUFFLE(2, 2, 2, 2))));
            // Add to translation column
            this->m_Columns[3] = _mm_add_ps(this->m_Columns[3], result);
        }
        else
        {
            // Scalar: t' = R * S * t
            this->m_03 += this->m_00 * t.m_X + this->m_01 * t.m_Y + this->m_02 * t.m_Z;
            this->m_13 += this->m_10 * t.m_X + this->m_11 * t.m_Y + this->m_12 * t.m_Z;
            this->m_23 += this->m_20 * t.m_X + this->m_21 * t.m_Y + this->m_22 * t.m_Z;
        }
        return *this;
        */
    }

    //------------------------------------------------------------------------------
    // Rotate
    //------------------------------------------------------------------------------
    //
    // Post-multiplies rotation from quaternion.
    //
    // Params:
    //  q - fromRotation quaternion.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::Rotate(const fquat& q) noexcept
    {
        return *this = fromRotation(q) * (*this);
    }

    //------------------------------------------------------------------------------
    // Rotate
    //------------------------------------------------------------------------------
    //
    // Post-multiplies rotation from axis-angle.
    //
    // Params:
    //  axis - Axis (normalized).
    //  angle - Angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::Rotate(const fvec3& axis, radian angle) noexcept
    {
        (*this) *= fromRotation(axis, angle);
        return *this;
    }


    //------------------------------------------------------------------------------
    template <bool V> //---
    inline fmat4_t<V>& fmat4_t<V>::PreRotateX(radian angle) noexcept
    {
        if (angle.m_Value == 0.0f) return *this;

        return (*this) *= fromRotationX(angle);

        /*
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            const floatx4 col1      = this->m_Columns[1]; // [m_01, m_11, m_21, m_31]
            const floatx4 col2      = this->m_Columns[2]; // [m_02, m_12, m_22, m_32]

            const floatx4 cs        = _mm_set1_ps(c);
            const floatx4 sn        = _mm_set1_ps(-s);
            const floatx4 sx        = _mm_set1_ps(s);

            // Update col1 and col2
            const floatx4 temp_col1 = _mm_add_ps(_mm_mul_ps(cs, col1), _mm_mul_ps(sx, col2));
            const floatx4 temp_col2 = _mm_add_ps(_mm_mul_ps(sn, col1), _mm_mul_ps(cs, col2));

            this->m_Columns[1] = temp_col1;
            this->m_Columns[2] = temp_col2;
        }
        else
        {
            float m01 = this->m_01;
            float m11 = this->m_11;
            float m21 = this->m_21;
            float m31 = this->m_31;
            this->m_01 = c * m01 + s * this->m_02; this->m_02 = c * this->m_02 - s * m01;
            this->m_11 = c * m11 + s * this->m_12; this->m_12 = c * this->m_12 - s * m11;
            this->m_21 = c * m21 + s * this->m_22; this->m_22 = c * this->m_22 - s * m21;
            this->m_31 = c * m31 + s * this->m_32; this->m_32 = c * this->m_32 - s * m31;
        }

        return *this;
        */
    }

    //------------------------------------------------------------------------------
    template <bool V> //---
    inline fmat4_t<V>& fmat4_t<V>::PreRotateY(radian angle) noexcept
    {
        if (angle.m_Value == 0.0f) return *this;

        return (*this) *= fromRotationY(angle);

        /*
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            floatx4 col0 = this->m_Columns[0]; // [m_00, m_10, m_20, m_30]
            floatx4 col2 = this->m_Columns[2]; // [m_02, m_12, m_22, m_32]

            floatx4 cs = _mm_set1_ps(c);
            floatx4 sn = _mm_set1_ps(-s);
            floatx4 sx = _mm_set1_ps(s);

            // Update col0 and col2
            floatx4 temp_col0 = _mm_add_ps(_mm_mul_ps(cs, col0), _mm_mul_ps(sn, col2));
            floatx4 temp_col2 = _mm_add_ps(_mm_mul_ps(sx, col0), _mm_mul_ps(cs, col2));

            this->m_Columns[0] = temp_col0;
            this->m_Columns[2] = temp_col2;
        }
        else
        {
            float m00 = this->m_00;
            float m10 = this->m_10;
            float m20 = this->m_20;
            float m30 = this->m_30;
            this->m_00 = c * m00 - s * this->m_02; this->m_02 = s * m00 + c * this->m_02;
            this->m_10 = c * m10 - s * this->m_12; this->m_12 = s * m10 + c * this->m_12;
            this->m_20 = c * m20 - s * this->m_22; this->m_22 = s * m20 + c * this->m_22;
            this->m_30 = c * m30 - s * this->m_32; this->m_32 = s * m30 + c * this->m_32;
        }

        return *this;
        */
    }

    //------------------------------------------------------------------------------
    template <bool V> //----
    inline fmat4_t<V>& fmat4_t<V>::PreRotateZ(radian angle) noexcept
    {
        if (angle.m_Value == 0.0f) return *this;

        return (*this) *= fromRotationZ(angle);

        /*
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            floatx4 col0 = this->m_Columns[0]; // [m_00, m_10, m_20, m_30]
            floatx4 col1 = this->m_Columns[1]; // [m_01, m_11, m_21, m_31]

            floatx4 cs = _mm_set1_ps(c);
            floatx4 sn = _mm_set1_ps(-s);
            floatx4 sx = _mm_set1_ps(s);

            // Update col0 and col1
            floatx4 temp_col0 = _mm_add_ps(_mm_mul_ps(cs, col0), _mm_mul_ps(sx, col1));
            floatx4 temp_col1 = _mm_add_ps(_mm_mul_ps(sn, col0), _mm_mul_ps(cs, col1));

            this->m_Columns[0] = temp_col0;
            this->m_Columns[1] = temp_col1;
        }
        else
        {
            float m00 = this->m_00;
            float m10 = this->m_10;
            float m20 = this->m_20;
            float m30 = this->m_30;
            this->m_00 = c * m00 + s * this->m_01; this->m_01 = c * this->m_01 - s * m00;
            this->m_10 = c * m10 + s * this->m_11; this->m_11 = c * this->m_11 - s * m10;
            this->m_20 = c * m20 + s * this->m_21; this->m_21 = c * this->m_21 - s * m20;
            this->m_30 = c * m30 + s * this->m_31; this->m_31 = c * this->m_31 - s * m30;
        }

        return *this;
        */
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::RotateX(radian angle) noexcept
    {
        if (angle.m_Value == 0.0f) return *this;

        return *this = fmat4_t::fromRotationX(angle) * (*this);
        /*
         float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            for (int i = 0; i < 4; ++i)
            {
                float y_val = this->m_Elements[4 * i + 1];
                float z_val = this->m_Elements[4 * i + 2];
                float new_y = c * y_val - s * z_val;
                float new_z = s * y_val + c * z_val;
                this->m_Elements[4 * i + 1] = new_y;
                this->m_Elements[4 * i + 2] = new_z;
            }
        }
        else
        {
            float y0 = this->m_10; float z0 = this->m_20;
            this->m_10 = c * y0 - s * z0; this->m_20 = s * y0 + c * z0;
            float y1 = this->m_11; float z1 = this->m_21;
            this->m_11 = c * y1 - s * z1; this->m_21 = s * y1 + c * z1;
            float y2 = this->m_12; float z2 = this->m_22;
            this->m_12 = c * y2 - s * z2; this->m_22 = s * y2 + c * z2;
            float y3 = this->m_13; float z3 = this->m_23;
            this->m_13 = c * y3 - s * z3; this->m_23 = s * y3 + c * z3;
        }

        return *this;
        */
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::RotateY(radian angle) noexcept
    {
        if (angle.m_Value == 0.0f) return *this;

        return *this = fmat4_t::fromRotationY(angle) * (*this);
        /*
        float c, s;
        xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            for (int i = 0; i < 4; ++i)
            {
                float x_val = this->m_Elements[4 * i + 0];
                float z_val = this->m_Elements[4 * i + 2];
                float new_x = c * x_val + s * z_val;
                float new_z = -s * x_val + c * z_val;
                this->m_Elements[4 * i + 0] = new_x;
                this->m_Elements[4 * i + 2] = new_z;
            }
        }
        else
        {
            float x0 = this->m_00; float z0 = this->m_20;
            this->m_00 = c * x0 + s * z0; this->m_20 = -s * x0 + c * z0;
            float x1 = this->m_01; float z1 = this->m_21;
            this->m_01 = c * x1 + s * z1; this->m_21 = -s * x1 + c * z1;
            float x2 = this->m_02; float z2 = this->m_22;
            this->m_02 = c * x2 + s * z2; this->m_22 = -s * x2 + c * z2;
            float x3 = this->m_03; float z3 = this->m_23;
            this->m_03 = c * x3 + s * z3; this->m_23 = -s * x3 + c * z3;
        }

        return *this;
        */
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::RotateZ(radian angle) noexcept
    {
        if (angle.m_Value == 0.0f) return *this;

        return *this = fmat4_t::fromRotationZ(angle) * (*this);

        /*
         float c, s;
         xmath::SinCos(angle, s, c);

        if constexpr (V)
        {
            for (int i = 0; i < 4; ++i)
            {
                float x_val = this->m_Elements[4 * i + 0];
                float y_val = this->m_Elements[4 * i + 1];
                float new_x = c * x_val - s * y_val;
                float new_y = s * x_val + c * y_val;
                this->m_Elements[4 * i + 0] = new_x;
                this->m_Elements[4 * i + 1] = new_y;
            }
        }
        else
        {
            float x0 = this->m_00; float y0 = this->m_10;
            this->m_00 = c * x0 - s * y0; this->m_10 = s * x0 + c * y0;
            float x1 = this->m_01; float y1 = this->m_11;
            this->m_01 = c * x1 - s * y1; this->m_11 = s * x1 + c * y1;
            float x2 = this->m_02; float y2 = this->m_12;
            this->m_02 = c * x2 - s * y2; this->m_12 = s * x2 + c * y2;
            float x3 = this->m_03; float y3 = this->m_13;
            this->m_03 = c * x3 - s * y3; this->m_13 = s * x3 + c * y3;
        }

        return *this;
        */
    }

    //------------------------------------------------------------------------------
    // Scale
    //------------------------------------------------------------------------------
    //
    // Post-multiplies scale.
    //
    // Params:
    //  s - Scale vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Asserts finite and non-zero s.
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::Scale(const fvec3& s) noexcept
    {
        assert(s.isFinite() && s.m_X != 0.0f && s.m_Y != 0.0f && s.m_Z != 0.0f);
        return *this = fromScale(s) * (*this);
    }

    //------------------------------------------------------------------------------
    // Scale
    //------------------------------------------------------------------------------
    //
    // Post-multiplies uniform scale.
    //
    // Params:
    //  s - Uniform scale.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::Scale(float s) noexcept
    {
        return Scale(fvec3(s));
    }

    //------------------------------------------------------------------------------
    // Pre-multiply mutable chaining
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // PreTranslate
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies translation.
    //
    // Params:
    //  t - Translation vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Asserts finite t.
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::PreTranslate(const fvec3& t) noexcept
    {
        assert(t.isFinite());
        fmat4_t<V> trans = fromTranslation(t);
        *this = (*this) * trans;
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotate
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies rotation from quaternion.
    //
    // Params:
    //  q - fromRotation quaternion.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::PreRotate(const fquat& q) noexcept
    {
        fmat4_t<V> rot = fromRotation(q);
        *this = (*this) * rot;
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotate
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies rotation from axis-angle.
    //
    // Params:
    //  axis - Axis (normalized).
    //  angle - Angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::PreRotate(const fvec3& axis, radian angle) noexcept
    {
        fmat4_t<V> rot = fromRotation(axis, angle);
        *this = (*this) * rot;
        return *this;
    }



    //------------------------------------------------------------------------------
    // PreScale
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies the matrix with a scale transformation.
    //
    // Params:
    //  s - Scale factors (x, y, z).
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::PreScale(const fvec3& s) noexcept
    {
        if constexpr (V)
        {
            this->m_Columns[0] = _mm_mul_ps(this->m_Columns[0], _mm_set_ps1(s.m_X));
            this->m_Columns[1] = _mm_mul_ps(this->m_Columns[1], _mm_set_ps1(s.m_Y));
            this->m_Columns[2] = _mm_mul_ps(this->m_Columns[2], _mm_set_ps1(s.m_Z));
            floatx4 scale_trans = _mm_set_ps(0.0f, s.m_Z, s.m_Y, s.m_X);
            this->m_Columns[3] = _mm_mul_ps(this->m_Columns[3], _mm_set_ps(1.0f, scale_trans.m128_f32[2], scale_trans.m128_f32[1], scale_trans.m128_f32[0]));
        }
        else
        {
            this->m_00 *= s.m_X; this->m_01 *= s.m_X; this->m_02 *= s.m_X; this->m_03 *= s.m_X;
            this->m_10 *= s.m_Y; this->m_11 *= s.m_Y; this->m_12 *= s.m_Y; this->m_13 *= s.m_Y;
            this->m_20 *= s.m_Z; this->m_21 *= s.m_Z; this->m_22 *= s.m_Z; this->m_23 *= s.m_Z;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreScale
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies uniform scale.
    //
    // Params:
    //  s - Uniform scale.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::PreScale(float s) noexcept
    {
        return PreScale(fvec3(s));
    }

    //------------------------------------------------------------------------------
    // Clear methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // ClearTranslation
    //------------------------------------------------------------------------------
    //
    // Clears translation component.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::ClearTranslation(void) noexcept
    {
        this->m_03 = 0.0f; this->m_13 = 0.0f; this->m_23 = 0.0f;
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClearRotation
    //------------------------------------------------------------------------------
    //
    // Clears rotation (keeps scale/translation).
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::ClearRotation(void) noexcept
    {
        fvec3 scale = ExtractScale();
        setupScale(scale);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClearScale
    //------------------------------------------------------------------------------
    //
    // Clears scale (keeps rotation/translation).
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat4_t<V>& fmat4_t<V>::ClearScale(void) noexcept
    {
        fquat rot = ExtractRotation();
        fvec3 trans = ExtractPosition();
        *this = fromRotation(rot) * fromTranslation(trans);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Immutable copies
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // TranslateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and post-translates.
    //
    // Params:
    //  t - Translation vector.
    //
    // Returns:
    //  Translated copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::TranslateCopy(const fvec3& t) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.Translate(t);
    }

    //------------------------------------------------------------------------------
    // RotateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and post-rotates from quaternion.
    //
    // Params:
    //  q - fromRotation quaternion.
    //
    // Returns:
    //  Rotated copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::RotateCopy(const fquat& q) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.Rotate(q);
    }

    //------------------------------------------------------------------------------
    // RotateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and post-rotates from axis-angle.
    //
    // Params:
    //  axis - Axis.
    //  angle - Angle.
    //
    // Returns:
    //  Rotated copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::RotateCopy(const fvec3& axis, radian angle) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.Rotate(axis, angle);
    }

    //------------------------------------------------------------------------------
    // ScaleCopy
    //------------------------------------------------------------------------------
    //
    // Copies and post-scales.
    //
    // Params:
    //  s - Scale vector.
    //
    // Returns:
    //  Scaled copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::ScaleCopy(const fvec3& s) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.Scale(s);
    }

    //------------------------------------------------------------------------------
    // PreTranslateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and pre-translates.
    //
    // Params:
    //  t - Translation vector.
    //
    // Returns:
    //  Translated copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::PreTranslateCopy(const fvec3& t) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.PreTranslate(t);
    }

    //------------------------------------------------------------------------------
    // PreRotateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and pre-rotates from quaternion.
    //
    // Params:
    //  q - fromRotation quaternion.
    //
    // Returns:
    //  Rotated copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::PreRotateCopy(const fquat& q) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.PreRotate(q);
    }

    //------------------------------------------------------------------------------
    // PreRotateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and pre-rotates from axis-angle.
    //
    // Params:
    //  axis - Axis.
    //  angle - Angle.
    //
    // Returns:
    //  Rotated copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::PreRotateCopy(const fvec3& axis, radian angle) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.PreRotate(axis, angle);
    }

    //------------------------------------------------------------------------------
    // PreScaleCopy
    //------------------------------------------------------------------------------
    //
    // Copies and pre-scales.
    //
    // Params:
    //  s - Scale vector.
    //
    // Returns:
    //  Scaled copy.
    //
    template <bool V>
    inline fmat4_t<V> fmat4_t<V>::PreScaleCopy(const fvec3& s) const noexcept
    {
        fmat4_t<V> m = *this;
        return m.PreScale(s);
    }

    //------------------------------------------------------------------------------
    // Safety and validation
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all elements are finite.
    //
    // Returns:
    //  True if finite.
    //
    template <bool V>
    inline bool fmat4_t<V>::isFinite(void) const noexcept
    {
        for (float e : this->m_Elements) if (!xmath::isFinite(e)) return false;
        return true;
    }

    //------------------------------------------------------------------------------
    // isIdentity
    //------------------------------------------------------------------------------
    //
    // Checks if approximately identity.
    //
    // Returns:
    //  True if identity (epsilon 1e-6).
    //
    template <bool V>
    inline bool fmat4_t<V>::isIdentity(void) const noexcept
    {
        return Equals(fromIdentity(), 1e-6f);
    }

    //------------------------------------------------------------------------------
    // SanityCheck
    //------------------------------------------------------------------------------
    //
    // Performs sanity asserts.
    //
    // Notes:
    //  Asserts finite, non-singular, normalized basis, orthogonal.
    //
    template <bool V>
    inline void fmat4_t<V>::SanityCheck(void) const noexcept
    {
        assert(isFinite());

        fvec3 Scale = ExtractScale();
        assert(Scale.LengthSq() > 0.01f);

        fquat Rotation(*this);
        assert(Rotation.isFinite());

        fvec3 Translation = ExtractPosition();
        assert(Translation.isFinite());

        fmat4 Test(Translation, Rotation, Scale);

        for (std::int32_t i = 0; i < 4 * 4; i++)
        {
            const float a = reinterpret_cast<const float*>(this)[i];
            const float b = reinterpret_cast<const float*>(&Test)[i];
            assert(xmath::Abs(a - b) < 0.1f);
        }

        // Additional checks
        assert(std::abs(Determinant()) > 1e-6f); // Not singular
        assert(Forward().isNormalized());
        assert(Up().isNormalized());
        assert(Right().isNormalized());
        assert(std::abs(Forward().Dot(Up())) < 1e-6f); // Ortho
    }
}
