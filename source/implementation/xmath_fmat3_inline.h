
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fmat3_t
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
    constexpr fmat3_t<V>::fmat3_t(float diagonal) noexcept
    {
        assert(xmath::isFinite(diagonal));
        if constexpr (V)
        {
            floatx4 diag = _mm_set1_ps(diagonal);
            floatx4 zero = _mm_setzero_ps();
            this->m_Columns[0] = _mm_blend_ps(zero, diag, 0x1); // 1 0 0 0
            this->m_Columns[1] = _mm_blend_ps(zero, diag, 0x2); // 0 1 0 0
            this->m_Columns[2] = _mm_blend_ps(zero, diag, 0x4); // 0 0 1 0
        }
        else
        {
            this->m_00 = diagonal;  this->m_01 = 0.0f;      this->m_02 = 0.0f;
            this->m_10 = 0.0f;      this->m_11 = diagonal;  this->m_12 = 0.0f;
            this->m_20 = 0.0f;      this->m_21 = 0.0f;      this->m_22 = diagonal;
        }
    }

    //------------------------------------------------------------------------------
    // fmat3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from array.
    //
    // Params:
    //  arr - Array of 9 floats in column-major order.
    //
    // Notes:
    //  No asserts; assumes valid input.
    //
    template <bool V>
    constexpr fmat3_t<V>::fmat3_t(const std::array<float, 9>& arr) noexcept
    {
        if constexpr (V)
        {
            this->m_Columns[0] = _mm_set_ps(0.0f, arr[2], arr[1], arr[0]);
            this->m_Columns[1] = _mm_set_ps(0.0f, arr[5], arr[4], arr[3]);
            this->m_Columns[2] = _mm_set_ps(0.0f, arr[8], arr[7], arr[6]);
        }
        else
        {
            this->m_Elements = arr;
        }
    }

    //------------------------------------------------------------------------------
    // fmat3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from span.
    //
    // Params:
    //  span - Span of 9 floats in column-major order.
    //
    // Notes:
    //  Assigns in column-major.
    //
    template <bool V>
    constexpr fmat3_t<V>::fmat3_t(std::span<const float, 9> span) noexcept
    {
        if constexpr (V)
        {
            this->m_Columns[0] = _mm_set_ps(0.0f, span[2], span[1], span[0]);
            this->m_Columns[1] = _mm_set_ps(0.0f, span[5], span[4], span[3]);
            this->m_Columns[2] = _mm_set_ps(0.0f, span[8], span[7], span[6]);
        }
        else
        {
            for (size_t i = 0; i < 9; ++i) this->m_Elements[i] = span[i];
        }
    }

    //------------------------------------------------------------------------------
    // fmat3_t
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
    inline fmat3_t<V>::fmat3_t(const fquat& q) noexcept
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
        }
        else
        {
            this->m_00 = 1 - 2 * (yy + zz); this->m_01 = 2 * (xy - zw); this->m_02 = 2 * (xz + yw);
            this->m_10 = 2 * (xy + zw); this->m_11 = 1 - 2 * (xx + zz); this->m_12 = 2 * (yz - xw);
            this->m_20 = 2 * (xz - yw); this->m_21 = 2 * (yz + xw); this->m_22 = 1 - 2 * (xx + yy);
        }
    }

    //------------------------------------------------------------------------------
    // fmat3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from Euler angles.
    //
    // Params:
    //  Euler - Euler angles in radians (pitch, yaw, roll).
    //
    // Notes:
    //  Delegates to quaternion constructor.
    //
    template <bool V>
    inline fmat3_t<V>::fmat3_t(const radian3& Euler) noexcept
        : fmat3_t(fquat(Euler))
    {
    }

    //------------------------------------------------------------------------------
    // fmat3_t
    //------------------------------------------------------------------------------
    //
    // Constructor from rotation and scale (RS).
    //
    // Params:
    //  rotation - Rotation quaternion.
    //  scale - Scale vector.
    //
    // Notes:
    //  Asserts finite inputs, normalized rotation, non-zero scale.
    //  Order: scale * rotation (post-multiply).
    //
    template <bool V>
    inline fmat3_t<V>::fmat3_t(const fquat& rotation, const fvec3& scale) noexcept
    {
        setup(rotation, scale);
    }

    //------------------------------------------------------------------------------
    // fmat3_t
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
    constexpr fmat3_t<V>::fmat3_t(const fmat3_t<!V>& other) noexcept
    {
        if constexpr (V)
        {
            this->m_Columns[0] = _mm_set_ps(0.0f, other.m_20, other.m_10, other.m_00);
            this->m_Columns[1] = _mm_set_ps(0.0f, other.m_21, other.m_11, other.m_01);
            this->m_Columns[2] = _mm_set_ps(0.0f, other.m_22, other.m_12, other.m_02);
        }
        else
        {
            this->m_00 = other.m_00; this->m_01 = other.m_01; this->m_02 = other.m_02;
            this->m_10 = other.m_10; this->m_11 = other.m_11; this->m_12 = other.m_12;
            this->m_20 = other.m_20; this->m_21 = other.m_21; this->m_22 = other.m_22;
        }
    }

    //------------------------------------------------------------------------------
    // fromIdentity
    //------------------------------------------------------------------------------
    //
    // Static identity matrix.
    //
    // Returns:
    //  Identity matrix.
    //
    template <bool V>
    constexpr fmat3_t<V> fmat3_t<V>::fromIdentity(void) noexcept
    {
        fmat3_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = 1.0f;  m.m_01 = 0.0f;  m.m_02 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = 1.0f;  m.m_12 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = 1.0f;
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
    //  Zero matrix.
    //
    template <bool V>
    constexpr fmat3_t<V> fmat3_t<V>::fromZero(void) noexcept
    {
        fmat3_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_setzero_ps();
            m.m_Columns[1] = _mm_setzero_ps();
            m.m_Columns[2] = _mm_setzero_ps();
        }
        else
        {
            m.m_00 = 0.0f;  m.m_01 = 0.0f;  m.m_02 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = 0.0f;  m.m_12 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = 0.0f;
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
    //  q - Rotation quaternion.
    //
    // Returns:
    //  Rotation matrix.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::fromRotation(const fquat& q) noexcept
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

        fmat3_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 2 * (xz - yw), 2 * (xy + zw), 1 - 2 * (yy + zz));
            m.m_Columns[1] = _mm_set_ps(0.0f, 2 * (yz + xw), 1 - 2 * (xx + zz), 2 * (xy - zw));
            m.m_Columns[2] = _mm_set_ps(0.0f, 1 - 2 * (xx + yy), 2 * (yz - xw), 2 * (xz + yw));
        }
        else
        {
            m.m_00 = 1 - 2 * (yy + zz); m.m_01 = 2 * (xy - zw); m.m_02 = 2 * (xz + yw);
            m.m_10 = 2 * (xy + zw); m.m_11 = 1 - 2 * (xx + zz); m.m_12 = 2 * (yz - xw);
            m.m_20 = 2 * (xz - yw); m.m_21 = 2 * (yz + xw); m.m_22 = 1 - 2 * (xx + yy);
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
    //  axis - Rotation axis (normalized).
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Rotation matrix.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::fromRotation(const fvec3& axis, radian angle) noexcept
    {
        assert(axis.isNormalized());
        assert(xmath::isFinite(angle.m_Value));

        return fmat3_t<V>(fquat(axis, angle));
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
    //  Sets diagonal to scale components.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::fromScale(const fvec3& s) noexcept
    {
        assert(s.isFinite() && s.m_X != 0.0f && s.m_Y != 0.0f && s.m_Z != 0.0f);
        fmat3_t<V> m;
        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, s.m_X);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, s.m_Y, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, s.m_Z, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = s.m_X; m.m_01 = 0.0f;  m.m_02 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = s.m_Y; m.m_12 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = s.m_Z;
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // fromRotationX
    //------------------------------------------------------------------------------
    //
    // Creates a rotation matrix around the X-axis.
    //
    // Params:
    //  Angle - Rotation angle in radians.
    //
    // Returns:
    //  Rotation matrix.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::fromRotationX(radian Angle) noexcept
    {
        fmat3_t<V> m;
        float c, s;
        xmath::SinCos(Angle, s, c);

        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
            m.m_Columns[1] = _mm_set_ps(0.0f, s, c, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, c, -s, 0.0f);
        }
        else
        {
            m.m_00 = 1.0f;  m.m_01 = 0.0f;  m.m_02 = 0.0f;
            m.m_10 = 0.0f;  m.m_11 = c;     m.m_12 = -s;
            m.m_20 = 0.0f;  m.m_21 = s;     m.m_22 = c;
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // fromRotationY
    //------------------------------------------------------------------------------
    //
    // Creates a rotation matrix around the Y-axis.
    //
    // Params:
    //  Angle - Rotation angle in radians.
    //
    // Returns:
    //  Rotation matrix.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::fromRotationY(radian Angle) noexcept
    {
        fmat3_t<V> m;
        float c, s;
        xmath::SinCos(Angle, s, c);

        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, s, 0.0f, c);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
            m.m_Columns[2] = _mm_set_ps(0.0f, c, 0.0f, -s);
        }
        else
        {
            m.m_00 = c;     m.m_01 = 0.0f;  m.m_02 = -s;
            m.m_10 = 0.0f;  m.m_11 = 1.0f;  m.m_12 = 0.0f;
            m.m_20 = s;     m.m_21 = 0.0f;  m.m_22 = c;
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // fromRotationZ
    //------------------------------------------------------------------------------
    //
    // Creates a rotation matrix around the Z-axis.
    //
    // Params:
    //  Angle - Rotation angle in radians.
    //
    // Returns:
    //  Rotation matrix.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::fromRotationZ(radian Angle) noexcept
    {
        fmat3_t<V> m;
        float c, s;
        xmath::SinCos(Angle, s, c);

        if constexpr (V)
        {
            m.m_Columns[0] = _mm_set_ps(0.0f, 0.0f, s, c);
            m.m_Columns[1] = _mm_set_ps(0.0f, 0.0f, c, -s);
            m.m_Columns[2] = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
        }
        else
        {
            m.m_00 = c;     m.m_01 = -s;    m.m_02 = 0.0f;
            m.m_10 = s;     m.m_11 = c;     m.m_12 = 0.0f;
            m.m_20 = 0.0f;  m.m_21 = 0.0f;  m.m_22 = 1.0f;
        }
        return m;
    }

    //------------------------------------------------------------------------------
    // Setup methods (mutable, replace)
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
    inline fmat3_t<V>& fmat3_t<V>::setupIdentity(void) noexcept
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
    inline fmat3_t<V>& fmat3_t<V>::setupZero(void) noexcept
    {
        *this = fromZero();
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupRotation
    //------------------------------------------------------------------------------
    //
    // Sets to rotation matrix from quaternion.
    //
    // Params:
    //  q - Rotation quaternion.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::setupRotation(const fquat& q) noexcept
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
    inline fmat3_t<V>& fmat3_t<V>::setupRotation(const radian3& euler) noexcept
    {
        *this = fmat3_t<V>(euler);
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
    inline fmat3_t<V>& fmat3_t<V>::setupScale(const fvec3& s) noexcept
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
    inline fmat3_t<V>& fmat3_t<V>::setupScale(float s) noexcept
    {
        *this = fromScale(fvec3(s));
        return *this;
    }

    //------------------------------------------------------------------------------
    // setup
    //------------------------------------------------------------------------------
    //
    // Sets to rotation-scale (RS) matrix.
    //
    // Params:
    //  rotation - Rotation quaternion.
    //  scale - Scale vector.
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Asserts finite inputs, normalized rotation, non-zero scale.
    //  Directly sets rotation, scales columns for efficiency.
    //  Order: scale * rotation (post-multiply).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::setup(const fquat& rotation, const fvec3& scale) noexcept
    {
        assert(rotation.isNormalized());
        assert(scale.isFinite() && scale.m_X != 0.0f && scale.m_Y != 0.0f && scale.m_Z != 0.0f);

        // Set rotation
        *this = fromRotation(rotation);

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
            this->m_00 *= scale.m_X; this->m_01 *= scale.m_X; this->m_02 *= scale.m_X;
            this->m_10 *= scale.m_Y; this->m_11 *= scale.m_Y; this->m_12 *= scale.m_Y;
            this->m_20 *= scale.m_Z; this->m_21 *= scale.m_Z; this->m_22 *= scale.m_Z;
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
    //  row - Row index (0-2).
    //
    // Returns:
    //  Row vector.
    //
    // Notes:
    //  Asserts valid row.
    //  Gathers from columns for row-major view.
    //
    template <bool V>
    constexpr fvec3 fmat3_t<V>::operator[](size_t row) const noexcept
    {
        assert(row < 3);
        if constexpr (V)
        {
            alignas(16) float col0[4], col1[4], col2[4];
            _mm_store_ps(col0, this->m_Columns[0]);
            _mm_store_ps(col1, this->m_Columns[1]);
            _mm_store_ps(col2, this->m_Columns[2]);
            return fvec3(col0[row], col1[row], col2[row]);
        }
        else
        {
            return fvec3(this->m_Cells[0][row], this->m_Cells[1][row], this->m_Cells[2][row]);
        }
    }

    //------------------------------------------------------------------------------
    // operator()
    //------------------------------------------------------------------------------
    //
    // Access element (row-major API).
    //
    // Params:
    //  row - Row index (0-2).
    //  col - Column index (0-2).
    //
    // Returns:
    //  Reference to element.
    //
    // Notes:
    //  Asserts valid indices; transposes internally.
    //
    template <bool V>
    constexpr float& fmat3_t<V>::operator()(size_t row, size_t col) noexcept
    {
        assert(row < 3 && col < 3);
        return this->m_Cells[col][row]; // Transposed access
    }

    //------------------------------------------------------------------------------
    // operator()
    //------------------------------------------------------------------------------
    //
    // Const access element (row-major API).
    //
    // Params:
    //  row - Row index (0-2).
    //  col - Column index (0-2).
    //
    // Returns:
    //  Const reference to element.
    //
    // Notes:
    //  Asserts valid indices; transposes internally.
    //
    template <bool V>
    constexpr const float& fmat3_t<V>::operator()(size_t row, size_t col) const noexcept
    {
        assert(row < 3 && col < 3);
        return this->m_Cells[col][row]; // Transposed access
    }

    //------------------------------------------------------------------------------
    // operator std::span<const float,9>
    //------------------------------------------------------------------------------
    //
    // Conversion to span of 9 floats.
    //
    // Returns:
    //  Span of the matrix elements in column-major order.
    //
    template <bool V>
    constexpr fmat3_t<V>::operator std::span<const float, 9>() const noexcept
    {
        return std::span<const float, 9>(this->m_Elements.data(), 9);
    }

    //------------------------------------------------------------------------------
    // operator fquat
    //------------------------------------------------------------------------------
    //
    // Conversion to quaternion.
    //
    // Returns:
    //  Extracted rotation quaternion.
    //
    template <bool V>
    inline fmat3_t<V>::operator xmath::fquat() const noexcept
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
    inline fmat3_t<V> fmat3_t<V>::operator+(const fmat3_t<V>& other) const noexcept
    {
        assert(this->isFinite());
        assert(other.isFinite());
        fmat3_t<V> result;
        if constexpr (V)
        {
            result.m_Columns[0] = _mm_add_ps(this->m_Columns[0], other.m_Columns[0]);
            result.m_Columns[1] = _mm_add_ps(this->m_Columns[1], other.m_Columns[1]);
            result.m_Columns[2] = _mm_add_ps(this->m_Columns[2], other.m_Columns[2]);
        }
        else
        {
            for (size_t i = 0; i < 9; ++i) result.m_Elements[i] = this->m_Elements[i] + other.m_Elements[i];
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
    inline fmat3_t<V> fmat3_t<V>::operator-(const fmat3_t<V>& other) const noexcept
    {
        assert(this->isFinite());
        assert(other.isFinite());
        fmat3_t<V> result;
        if constexpr (V)
        {
            result.m_Columns[0] = _mm_sub_ps(this->m_Columns[0], other.m_Columns[0]);
            result.m_Columns[1] = _mm_sub_ps(this->m_Columns[1], other.m_Columns[1]);
            result.m_Columns[2] = _mm_sub_ps(this->m_Columns[2], other.m_Columns[2]);
        }
        else
        {
            for (size_t i = 0; i < 9; ++i) result.m_Elements[i] = this->m_Elements[i] - other.m_Elements[i];
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
    inline fmat3_t<V> fmat3_t<V>::operator*(const fmat3_t<V>& other) const noexcept
    {
        assert(this->isFinite());
        assert(other.isFinite());

        fmat3_t<V> result;
        if constexpr (V)
        {
            // Optimized SSE version (3 columns)
            floatx4 col0 = other.m_Columns[0];
            floatx4 bc0 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(0, 0, 0, 0));
            floatx4 bc1 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(1, 1, 1, 1));
            floatx4 bc2 = _mm_shuffle_ps(col0, col0, _MM_SHUFFLE(2, 2, 2, 2));
            result.m_Columns[0] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_mul_ps(this->m_Columns[2], bc2));

            floatx4 col1 = other.m_Columns[1];
            bc0 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(0, 0, 0, 0));
            bc1 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(1, 1, 1, 1));
            bc2 = _mm_shuffle_ps(col1, col1, _MM_SHUFFLE(2, 2, 2, 2));
            result.m_Columns[1] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_mul_ps(this->m_Columns[2], bc2));

            floatx4 col2 = other.m_Columns[2];
            bc0 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(0, 0, 0, 0));
            bc1 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(1, 1, 1, 1));
            bc2 = _mm_shuffle_ps(col2, col2, _MM_SHUFFLE(2, 2, 2, 2));
            result.m_Columns[2] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(this->m_Columns[0], bc0), _mm_mul_ps(this->m_Columns[1], bc1)),
                _mm_mul_ps(this->m_Columns[2], bc2));
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    result(i, j) = 0.0f;
                    for (int k = 0; k < 3; ++k)
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
    inline fmat3_t<V>& fmat3_t<V>::operator+=(const fmat3_t<V>& other) noexcept
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
    inline fmat3_t<V>& fmat3_t<V>::operator-=(const fmat3_t<V>& other) noexcept
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
    inline fmat3_t<V>& fmat3_t<V>::operator*=(const fmat3_t<V>& other) noexcept
    {
        *this = *this * other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Matrix-vector multiplication (vec3).
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
    inline fvec3 fmat3_t<V>::operator*(const fvec3& v) const noexcept
    {
        assert(this->isFinite());
        assert(v.isFinite());

        if constexpr (V)
        {
            floatx4 vec = _mm_set_ps(0.0f, v.m_Z, v.m_Y, v.m_X);
            floatx4 result = _mm_setzero_ps();
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[0], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0, 0, 0, 0))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[1], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(1, 1, 1, 1))));
            result = _mm_add_ps(result, _mm_mul_ps(this->m_Columns[2], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(2, 2, 2, 2))));
            return fvec3(result.m128_f32[0], result.m128_f32[1], result.m128_f32[2]);
        }
        else
        {
            return fvec3
            (this->m_00 * v.m_X + this->m_01 * v.m_Y + this->m_02 * v.m_Z
                , this->m_10 * v.m_X + this->m_11 * v.m_Y + this->m_12 * v.m_Z
                , this->m_20 * v.m_X + this->m_21 * v.m_Y + this->m_22 * v.m_Z
            );
        }
    }

    //------------------------------------------------------------------------------
    // Equals
    //------------------------------------------------------------------------------
    //
    // Checks if this matrix is approximately equal to another within a tolerance.
    //
    // Params:
    //  other - The matrix to compare.
    //  tolerance - The epsilon tolerance.
    //
    // Returns:
    //  True if equal within tolerance, false otherwise.
    //
    template <bool V>
    inline bool fmat3_t<V>::Equals(const fmat3_t& other, float tolerance) const noexcept
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (std::abs((*this)(i, j) - other(i, j)) > tolerance) return false;
            }
        }
        return true;
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
    //  Optimized for column-major storage.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::Transpose(void) const noexcept
    {
        assert(this->isFinite());
        fmat3_t<V> result;
        if constexpr (V)
        {
            alignas(16) float col0[4], col1[4], col2[4];
            _mm_store_ps(col0, this->m_Columns[0]);
            _mm_store_ps(col1, this->m_Columns[1]);
            _mm_store_ps(col2, this->m_Columns[2]);

            float c00 = col0[0]; // m_00
            float c01 = col1[0]; // m_01
            float c02 = col2[0]; // m_02

            float c10 = col0[1]; // m_10
            float c11 = col1[1]; // m_11
            float c12 = col2[1]; // m_12

            float c20 = col0[2]; // m_20
            float c21 = col1[2]; // m_21
            float c22 = col2[2]; // m_22

            // Set transposed columns
            result.m_Columns[0] = _mm_set_ps(0.0f, c02, c01, c00);
            result.m_Columns[1] = _mm_set_ps(0.0f, c12, c11, c10);
            result.m_Columns[2] = _mm_set_ps(0.0f, c22, c21, c20);
        }
        else
        {
            for (size_t i = 0; i < 3; ++i)
                for (size_t j = 0; j < 3; ++j)
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
    inline float fmat3_t<V>::Determinant(void) const noexcept
    {
        assert(this->isFinite());

        // Scalar implementation for both
        return this->m_00 * (this->m_11 * this->m_22 - this->m_12 * this->m_21) -
            this->m_01 * (this->m_10 * this->m_22 - this->m_12 * this->m_20) +
            this->m_02 * (this->m_10 * this->m_21 - this->m_11 * this->m_20);
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
    //  Uses adjoint / det for 3x3.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::Inverse(void) const noexcept
    {
        assert(this->isFinite());

        float det = Determinant();
        assert(std::abs(det) >= 0.00001f);
        float inv_det = 1.0f / det;

        fmat3_t<V> result;
        if constexpr (V)
        {
            // Cofactors (adjoint / det)
            result.m_Columns[0] = _mm_set_ps(0.0f,
                inv_det * (this->m_10 * this->m_21 - this->m_11 * this->m_20),  // C20
                -inv_det * (this->m_10 * this->m_22 - this->m_12 * this->m_20), // C10
                inv_det * (this->m_11 * this->m_22 - this->m_12 * this->m_21)); // C00

            result.m_Columns[1] = _mm_set_ps(0.0f,
                -inv_det * (this->m_00 * this->m_21 - this->m_01 * this->m_20),  // C21
                inv_det * (this->m_00 * this->m_22 - this->m_02 * this->m_20),   // C11
                -inv_det * (this->m_01 * this->m_22 - this->m_02 * this->m_21)); // C01

            result.m_Columns[2] = _mm_set_ps(0.0f,
                inv_det * (this->m_00 * this->m_11 - this->m_01 * this->m_10),   // C22
                -inv_det * (this->m_00 * this->m_12 - this->m_02 * this->m_10),  // C12
                inv_det * (this->m_01 * this->m_12 - this->m_02 * this->m_11));  // C02
        }
        else
        {
            // Compute cofactors
            result.m_00 = inv_det * (this->m_11 * this->m_22 - this->m_12 * this->m_21);
            result.m_01 = -inv_det * (this->m_01 * this->m_22 - this->m_02 * this->m_21);
            result.m_02 = inv_det * (this->m_01 * this->m_12 - this->m_02 * this->m_11);

            result.m_10 = -inv_det * (this->m_10 * this->m_22 - this->m_12 * this->m_20);
            result.m_11 = inv_det * (this->m_00 * this->m_22 - this->m_02 * this->m_20);
            result.m_12 = -inv_det * (this->m_00 * this->m_12 - this->m_02 * this->m_10);

            result.m_20 = inv_det * (this->m_10 * this->m_21 - this->m_11 * this->m_20);
            result.m_21 = -inv_det * (this->m_00 * this->m_21 - this->m_01 * this->m_20);
            result.m_22 = inv_det * (this->m_00 * this->m_11 - this->m_01 * this->m_10);
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
    inline fmat3_t<V>& fmat3_t<V>::Orthogonalize(void) noexcept
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
    //  Assumes orthogonal basis (RS matrix); for stability, normalizes columns.
    //  Preserves mirroring if present.
    //  If input from quaternion, extracts back the same (up to sign).
    //
    template <bool V>
    inline fquat fmat3_t<V>::ExtractRotation(void) const noexcept
    {
        fvec3 scale = ExtractScale();

        fmat3_t<V> normalized = *this;
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
            normalized.m_00 /= scale.m_X; normalized.m_01 /= scale.m_X; normalized.m_02 /= scale.m_X;
            normalized.m_10 /= scale.m_Y; normalized.m_11 /= scale.m_Y; normalized.m_12 /= scale.m_Y;
            normalized.m_20 /= scale.m_Z; normalized.m_21 /= scale.m_Z; normalized.m_22 /= scale.m_Z;
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
    inline fvec3 fmat3_t<V>::ExtractScale(void) const noexcept
    {
        return fvec3(fvec3(this->m_00, this->m_10, this->m_20).Length(),
            fvec3(this->m_01, this->m_11, this->m_21).Length(),
            fvec3(this->m_02, this->m_12, this->m_22).Length());
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
    inline fvec3 fmat3_t<V>::Forward(void) const noexcept
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
    inline fvec3 fmat3_t<V>::Back(void) const noexcept
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
    inline fvec3 fmat3_t<V>::Up(void) const noexcept
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
    inline fvec3 fmat3_t<V>::Down(void) const noexcept
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
    inline fvec3 fmat3_t<V>::Left(void) const noexcept
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
    inline fvec3 fmat3_t<V>::Right(void) const noexcept
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
    //  Transformed vector.
    //
    // Notes:
    //  Uses the 3x3 submatrix.
    //  Optimized for SIMD and scalar.
    //
    template <bool V>
    inline fvec3 fmat3_t<V>::RotateVector(const fvec3& v) const noexcept
    {
        if constexpr (V)
        {
            floatx4 vec = _mm_set_ps(0.0f, v.m_Z, v.m_Y, v.m_X);
            floatx4 result = _mm_setzero_ps();
            result = _mm_add_ps(_mm_add_ps(
                _mm_mul_ps(this->m_Columns[0], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0, 0, 0, 0))),
                _mm_mul_ps(this->m_Columns[1], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(1, 1, 1, 1)))),
                _mm_mul_ps(this->m_Columns[2], _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(2, 2, 2, 2))));
            return fvec3(result.m128_f32[0], result.m128_f32[1], result.m128_f32[2]);
        }
        else
        {
            float x = this->m_00 * v.m_X + this->m_01 * v.m_Y + this->m_02 * v.m_Z;
            float y = this->m_10 * v.m_X + this->m_11 * v.m_Y + this->m_12 * v.m_Z;
            float z = this->m_20 * v.m_X + this->m_21 * v.m_Y + this->m_22 * v.m_Z;
            return fvec3(x, y, z);
        }
    }

    //------------------------------------------------------------------------------
    // InvRotateVector
    //------------------------------------------------------------------------------
    //
    // Transforms a vector by the inverse of the 3x3 rotation/scale part.
    //
    // Params:
    //  v - Vector to transform.
    //
    // Returns:
    //  Transformed vector.
    //
    // Notes:
    //  Uses inverse of 3x3 submatrix.
    //  Optimized for SIMD and scalar.
    //
    template <bool V>
    inline fvec3 fmat3_t<V>::InvRotateVector(const fvec3& v) const noexcept
    {
        float det = Determinant();
        assert(std::abs(det) >= 0.00001f);
        float inv_det = 1.0f / det;

        if constexpr (V)
        {
            // Cofactors
            floatx4 inv_row0 = _mm_set_ps(0.0f,
                inv_det * (this->m_10 * this->m_21 - this->m_11 * this->m_20),  // C20
                -inv_det * (this->m_10 * this->m_22 - this->m_12 * this->m_20), // C10
                inv_det * (this->m_11 * this->m_22 - this->m_12 * this->m_21)); // C00

            floatx4 inv_row1 = _mm_set_ps(0.0f,
                -inv_det * (this->m_00 * this->m_21 - this->m_01 * this->m_20),  // C21
                inv_det * (this->m_00 * this->m_22 - this->m_02 * this->m_20),   // C11
                -inv_det * (this->m_01 * this->m_22 - this->m_02 * this->m_21)); // C01

            floatx4 inv_row2 = _mm_set_ps(0.0f,
                inv_det * (this->m_00 * this->m_11 - this->m_01 * this->m_10),   // C22
                -inv_det * (this->m_00 * this->m_12 - this->m_02 * this->m_10),  // C12
                inv_det * (this->m_01 * this->m_12 - this->m_02 * this->m_11));  // C02

            // Apply to vector
            floatx4 vec = _mm_set_ps(0.0f, v.m_Z, v.m_Y, v.m_X);
            floatx4 result = _mm_setzero_ps();
            result = _mm_add_ps(_mm_add_ps(
                _mm_mul_ps(inv_row0, _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0, 0, 0, 0))),
                _mm_mul_ps(inv_row1, _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(1, 1, 1, 1)))),
                _mm_mul_ps(inv_row2, _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(2, 2, 2, 2))));
            return fvec3(result.m128_f32[0], result.m128_f32[1], result.m128_f32[2]);
        }
        else
        {
            return fvec3(
                inv_det * (this->m_11 * this->m_22 - this->m_12 * this->m_21) * v.m_X +
                inv_det * (this->m_02 * this->m_21 - this->m_01 * this->m_22) * v.m_Y +
                inv_det * (this->m_01 * this->m_12 - this->m_02 * this->m_11) * v.m_Z,
                inv_det * (this->m_12 * this->m_20 - this->m_10 * this->m_22) * v.m_X +
                inv_det * (this->m_00 * this->m_22 - this->m_02 * this->m_20) * v.m_Y +
                inv_det * (this->m_02 * this->m_10 - this->m_00 * this->m_12) * v.m_Z,
                inv_det * (this->m_10 * this->m_21 - this->m_11 * this->m_20) * v.m_X +
                inv_det * (this->m_01 * this->m_20 - this->m_00 * this->m_21) * v.m_Y +
                inv_det * (this->m_00 * this->m_11 - this->m_01 * this->m_10) * v.m_Z
            );
        }
    }

    //------------------------------------------------------------------------------
    // TransformDirection
    //------------------------------------------------------------------------------
    //
    // Transforms direction (rotation/scale).
    //
    // Params:
    //  d - Direction to transform.
    //
    // Returns:
    //  Transformed direction.
    //
    template <bool V>
    inline fvec3 fmat3_t<V>::TransformDirection(const fvec3& d) const noexcept
    {
        return RotateVector(d);
    }

    //------------------------------------------------------------------------------
    // Mutable chaining methods (post-multiply)
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Rotate
    //------------------------------------------------------------------------------
    //
    // Post-multiplies rotation from quaternion.
    //
    // Params:
    //  q - Rotation quaternion.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::Rotate(const fquat& q) noexcept
    {
        *this = fromRotation(q) * (*this);
        return *this;
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
    inline fmat3_t<V>& fmat3_t<V>::Rotate(const fvec3& axis, radian angle) noexcept
    {
        *this *= fromRotation(axis, angle);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateX
    //------------------------------------------------------------------------------
    //
    // Post-multiplies rotation around the X-axis.
    //
    // Params:
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::RotateX(radian angle) noexcept
    {
        *this = fromRotationX(angle) * (*this);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateY
    //------------------------------------------------------------------------------
    //
    // Post-multiplies rotation around the Y-axis.
    //
    // Params:
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::RotateY(radian angle) noexcept
    {
        *this = fromRotationY(angle) * (*this);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateZ
    //------------------------------------------------------------------------------
    //
    // Post-multiplies rotation around the Z-axis.
    //
    // Params:
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::RotateZ(radian angle) noexcept
    {
        *this = fromRotationZ(angle) * (*this);
        return *this;
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
    inline fmat3_t<V>& fmat3_t<V>::Scale(const fvec3& s) noexcept
    {
        assert(s.isFinite() && s.m_X != 0.0f && s.m_Y != 0.0f && s.m_Z != 0.0f);
        *this = fromScale(s) * (*this);
        return *this;
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
    inline fmat3_t<V>& fmat3_t<V>::Scale(float s) noexcept
    {
        return Scale(fvec3(s));
    }

    //------------------------------------------------------------------------------
    // Pre-multiply mutable chaining
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // PreRotate
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies rotation from quaternion.
    //
    // Params:
    //  q - Rotation quaternion.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::PreRotate(const fquat& q) noexcept
    {
        *this *= fromRotation(q);
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
    inline fmat3_t<V>& fmat3_t<V>::PreRotate(const fvec3& axis, radian angle) noexcept
    {
        *this *= fromRotation(axis, angle);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateX
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies rotation around the X-axis.
    //
    // Params:
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::PreRotateX(radian angle) noexcept
    {
        *this *= fromRotationX(angle);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateY
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies rotation around the Y-axis.
    //
    // Params:
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::PreRotateY(radian angle) noexcept
    {
        *this *= fromRotationY(angle);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateZ
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies rotation around the Z-axis.
    //
    // Params:
    //  angle - Rotation angle in radians.
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::PreRotateZ(radian angle) noexcept
    {
        *this *= fromRotationZ(angle);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreScale
    //------------------------------------------------------------------------------
    //
    // Pre-multiplies scale.
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
    inline fmat3_t<V>& fmat3_t<V>::PreScale(const fvec3& s) noexcept
    {
        assert(s.isFinite() && s.m_X != 0.0f && s.m_Y != 0.0f && s.m_Z != 0.0f);
        if constexpr (V)
        {
            this->m_Columns[0] = _mm_mul_ps(this->m_Columns[0], _mm_set1_ps(s.m_X));
            this->m_Columns[1] = _mm_mul_ps(this->m_Columns[1], _mm_set1_ps(s.m_Y));
            this->m_Columns[2] = _mm_mul_ps(this->m_Columns[2], _mm_set1_ps(s.m_Z));
        }
        else
        {
            this->m_00 *= s.m_X; this->m_01 *= s.m_X; this->m_02 *= s.m_X;
            this->m_10 *= s.m_Y; this->m_11 *= s.m_Y; this->m_12 *= s.m_Y;
            this->m_20 *= s.m_Z; this->m_21 *= s.m_Z; this->m_22 *= s.m_Z;
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
    inline fmat3_t<V>& fmat3_t<V>::PreScale(float s) noexcept
    {
        return PreScale(fvec3(s));
    }

    //------------------------------------------------------------------------------
    // Clear methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // ClearRotation
    //------------------------------------------------------------------------------
    //
    // Clears rotation (keeps scale).
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::ClearRotation(void) noexcept
    {
        fvec3 scale = ExtractScale();
        setupScale(scale);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClearScale
    //------------------------------------------------------------------------------
    //
    // Clears scale (keeps rotation).
    //
    // Returns:
    //  Reference (chainable).
    //
    template <bool V>
    inline fmat3_t<V>& fmat3_t<V>::ClearScale(void) noexcept
    {
        fquat rot = ExtractRotation();
        setupRotation(rot);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Immutable versions (Copy suffix)
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // RotateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and post-rotates from quaternion.
    //
    // Params:
    //  q - Rotation quaternion.
    //
    // Returns:
    //  Rotated copy.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::RotateCopy(const fquat& q) const noexcept
    {
        fmat3_t<V> m = *this;
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
    inline fmat3_t<V> fmat3_t<V>::RotateCopy(const fvec3& axis, radian angle) const noexcept
    {
        fmat3_t<V> m = *this;
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
    inline fmat3_t<V> fmat3_t<V>::ScaleCopy(const fvec3& s) const noexcept
    {
        fmat3_t<V> m = *this;
        return m.Scale(s);
    }

    //------------------------------------------------------------------------------
    // PreRotateCopy
    //------------------------------------------------------------------------------
    //
    // Copies and pre-rotates from quaternion.
    //
    // Params:
    //  q - Rotation quaternion.
    //
    // Returns:
    //  Rotated copy.
    //
    template <bool V>
    inline fmat3_t<V> fmat3_t<V>::PreRotateCopy(const fquat& q) const noexcept
    {
        fmat3_t<V> m = *this;
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
    inline fmat3_t<V> fmat3_t<V>::PreRotateCopy(const fvec3& axis, radian angle) const noexcept
    {
        fmat3_t<V> m = *this;
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
    inline fmat3_t<V> fmat3_t<V>::PreScaleCopy(const fvec3& s) const noexcept
    {
        fmat3_t<V> m = *this;
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
    inline bool fmat3_t<V>::isFinite(void) const noexcept
    {
        for (size_t i = 0; i < (V ? 12 : 9); ++i) if (!xmath::isFinite(this->m_Elements[i])) return false;
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
    inline bool fmat3_t<V>::isIdentity(void) const noexcept
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
    inline void fmat3_t<V>::SanityCheck(void) const noexcept
    {
        assert(isFinite());

        fvec3 scale = ExtractScale();
        assert(scale.LengthSq() > 0.01f);

        fquat rotation = ExtractRotation();
        assert(rotation.isFinite());

        fmat3_t test;
        test.setup(rotation, scale);

        for (size_t i = 0; i < 3; ++i)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                assert(std::abs((*this)(i, j) - test(i, j)) < 0.1f);
            }
        }

        // Additional checks
        assert(std::abs(Determinant()) > 1e-6f); // Not singular
        assert(Forward().isNormalized());
        assert(Up().isNormalized());
        assert(Right().isNormalized());
        assert(std::abs(Forward().Dot(Up())) < 1e-6f); // Ortho
        assert(std::abs(Forward().Dot(Right())) < 1e-6f);
        assert(std::abs(Up().Dot(Right())) < 1e-6f);
    }
}

