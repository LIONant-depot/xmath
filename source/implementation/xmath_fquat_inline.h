#pragma once
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from individual components.
    //
    // Params:
    //  x - X component (imaginary i).
    //  y - Y component (imaginary j).
    //  z - Z component (imaginary k).
    //  w - W component (real part).
    //
    // Notes:
    //  Uses _mm_set_ps for SIMD; assumes right-hand rule.
    //  Asserts finite components in debug.
    //
    template <bool V >
    constexpr fquat_t<V>::fquat_t(float x, float y, float z, float w) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(w, z, y, x);
        }
        else
        {
            this->m_X = x;
            this->m_Y = y;
            this->m_Z = z;
            this->m_W = w;
        }
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from axis and angle.
    //
    // Params:
    //  axis - Normalized rotation axis.
    //  angle - Rotation angle in radians.
    //
    // Notes:
    //  Uses half-angle trig: w = cos(angle/2), xyz = axis * sin(angle/2).
    //  Assumes normalized axis; asserts finite and unit length.
    //
    template <bool V >
    constexpr fquat_t<V>::fquat_t(const fvec3& axis, radian angle) noexcept
    {
        *this = fromAxisAngle(axis, angle);
    }

    //------------------------------------------------------------------------------
    // SetRotationX
    //------------------------------------------------------------------------------
    //
    // Sets the quaternion to a rotation around the X-axis.
    //
    // Params:
    //  rx - Rotation angle around X.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses half-angle: x = sin(rx/2), y=0, z=0, w=cos(rx/2).
    //  Assumes right-hand rule; asserts finite rx.
    //
    template <bool V >
    inline fquat_t<V>& fquat_t<V>::setupRotationX(radian rx) noexcept
    {
        radian half = rx * 0.5f;
        float s, c;
        xmath::SinCos(half, s, c);
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(c, 0.0f, 0.0f, s);
        }
        else
        {
            this->m_X = s;
            this->m_Y = 0.0f;
            this->m_Z = 0.0f;
            this->m_W = c;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // SetRotationY
    //------------------------------------------------------------------------------
    //
    // Sets the quaternion to a rotation around the Y-axis.
    //
    // Params:
    //  ry - Rotation angle around Y.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses half-angle: x=0, y=sin(ry/2), z=0, w=cos(ry/2).
    //  Assumes right-hand rule; asserts finite ry.
    //
    template <bool V >
    inline fquat_t<V>& fquat_t<V>::setupRotationY(radian ry) noexcept
    {
        radian half = ry * 0.5f;
        float s, c;
        xmath::SinCos(half, s, c);
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(c, 0.0f, s, 0.0f);
        }
        else
        {
            this->m_X = 0.0f;
            this->m_Y = s;
            this->m_Z = 0.0f;
            this->m_W = c;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // SetRotationZ
    //------------------------------------------------------------------------------
    //
    // Sets the quaternion to a rotation around the Z-axis.
    //
    // Params:
    //  rz - Rotation angle around Z.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses half-angle: x=0, y=0, z=sin(rz/2), w=cos(rz/2).
    //  Assumes right-hand rule; asserts finite rz.
    //
    template <bool V >
    inline fquat_t<V>& fquat_t<V>::setupRotationZ(radian rz) noexcept
    {
        radian half = rz * 0.5f;
        float s, c;
        xmath::SinCos(half, s, c);
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(c, s, 0.0f, 0.0f);
        }
        else
        {
            this->m_X = 0.0f;
            this->m_Y = 0.0f;
            this->m_Z = s;
            this->m_W = c;
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from Euler angles.
    //
    // Params:
    //  euler - Euler angles in radians (pitch x, yaw y, roll z).
    //
    // Notes:
    //  Uses YXZ order (yaw, pitch, roll) for right-hand system.
    //  Multiplies individual axis quaternions.
    //  Asserts finite angles.
    //
    template <bool V >
    inline fquat_t<V>::fquat_t(const radian3& euler) noexcept
    {
        fquat_t qx = setupRotationX(euler.m_Pitch);
        fquat_t qy = setupRotationY(euler.m_Yaw);
        fquat_t qz = setupRotationZ(euler.m_Roll);
        *this = qy * qx * qz; // YXZ order
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from from-to vectors.
    //
    // Params:
    //  from - Starting direction (normalized).
    //  to - Target direction (normalized).
    //
    // Notes:
    //  Computes rotation aligning from to to using cross and dot.
    //  Handles 180 deg case with perpendicular axis.
    //  Assumes normalized vectors; asserts finite and unit.
    //
    template <bool V >
    inline fquat_t<V>::fquat_t(const fvec3& from, const fvec3& to, const fvec3& up ) noexcept
    {
        *this = FromToRotation(from, to, up);
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor for look rotation.
    //
    // Params:
    //  forward - Forward direction.
    //  up - Up direction (default Up()).
    //
    // Notes:
    //  Normalizes forward/up; computes right and adjusted up.
    //  Builds matrix then to quaternion (or direct quat).
    //  Asserts finite and non-zero.
    //
    template <bool V >
    inline fquat_t<V>::fquat_t(const fvec3& forward, const fvec3& up) noexcept
    {
        *this = LookRotation(forward, up);
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from SIMD register.
    //
    // Params:
    //  reg - __m128 register with x,y,z,w.
    //
    // Notes:
    //  Direct assignment.
    //
    template <bool V >
    constexpr fquat_t<V>::fquat_t(const floatx4& reg) noexcept requires V
        : details::f_quat::simd_data{ .m_XYZW = reg }
    {}

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from other SIMD variant.
    //
    // Params:
    //  other - Other fquat_t.
    //
    // Notes:
    //  Copies components.
    //
    template <bool V >
    constexpr fquat_t<V>::fquat_t(const fquat_t<!V>& other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZW = _mm_set_ps(other.m_W, other.m_Z, other.m_Y, other.m_X);
        }
        else
        {
            this->m_X = other.m_X;
            this->m_Y = other.m_Y;
            this->m_Z = other.m_Z;
            this->m_W = other.m_W;
        }
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // Constructor from double array.
    //
    // Params:
    //  conversion - Array of 4 doubles.
    //
    // Notes:
    //  Static cast to float.
    //
    template <bool V >
    constexpr fquat_t<V>::fquat_t(const std::array<double, 4>& conversion) noexcept
    {
        this->m_X = static_cast<float>(conversion[0]);
        this->m_Y = static_cast<float>(conversion[1]);
        this->m_Z = static_cast<float>(conversion[2]);
        this->m_W = static_cast<float>(conversion[3]);
    }

    //------------------------------------------------------------------------------
    // operator std::array<double,4>
    //------------------------------------------------------------------------------
    //
    // Converts to double array.
    //
    // Returns:
    //  {x, y, z, w} as doubles.
    //
    template <bool V >
    constexpr fquat_t<V>::operator std::array<double, 4>() const noexcept
    {
        return { static_cast<double>(this->m_X), static_cast<double>(this->m_Y), static_cast<double>(this->m_Z), static_cast<double>(this->m_W) };
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fromIdentity
    //------------------------------------------------------------------------------
    //
    // Returns identity quaternion (0,0,0,1).
    //
    template <bool V >
    constexpr fquat_t<V> fquat_t<V>::fromIdentity(void) noexcept
    {
        if constexpr (V) return fquat_t{ floatx4{.m128_f32{ 0,0,0,1 }} };
        else             return{ 0.0f, 0.0f, 0.0f, 1.0f };
    }

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Returns zero quaternion (0,0,0,0).
    //
    template <bool V >
    constexpr fquat_t<V> fquat_t<V>::fromZero(void) noexcept
    {
        if constexpr (V) return fquat_t{ floatx4{.m128_f32{ 0,0,0,0 }} };
        else             return{ 0.0f, 0.0f, 0.0f, 0.0f };
    }

    //------------------------------------------------------------------------------

    template <bool V >
    inline fquat_t<V> fquat_t<V>::fromAxisAngle(const fvec3& axis, radian angle) noexcept
    {
        radian half_angle = angle * 0.5f;
        float s, c;
        xmath::SinCos(half_angle, s, c);
        if constexpr (V)
        {
            __m128 sin_vec = _mm_set1_ps(s);
            __m128 axis_vec = _mm_mul_ps(axis.m_XYZW, sin_vec);
            return fquat_t{ _mm_insert_ps(axis_vec, _mm_set_ss(c), 0x30) }; // Insert c into w
        }
        else
        {
            return
            { axis.m_X * s
            , axis.m_Y* s
            , axis.m_Z* s
            , c
            };
        }

    }

    //------------------------------------------------------------------------------
    // Static methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes dot product of two quaternions.
    //
    // Params:
    //  a - First quaternion.
    //  b - Second quaternion.
    //
    // Returns:
    //  x1*x2 + y1*y2 + z1*z2 + w1*w2.
    //
    // Notes:
    //  Uses _mm_dp_ps for SIMD.
    //
    template <bool V >
    inline float fquat_t<V>::Dot(const fquat_t& a, const fquat_t& b) noexcept
    {
        if constexpr (V)
        {
            return _mm_cvtss_f32(_mm_dp_ps(a.m_XYZW, b.m_XYZW, 0xFF));
        }
        else
        {
            return a.m_X * b.m_X + a.m_Y * b.m_Y + a.m_Z * b.m_Z + a.m_W * b.m_W;
        }
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Linear interpolation between two quaternions.
    //
    // Params:
    //  a - Start quaternion.
    //  b - End quaternion.
    //  t - Interpolation factor (clamped 0-1).
    //
    // Returns:
    //  Normalized lerp result.
    //
    // Notes:
    //  Uses dot to choose shortest path; normalizes result.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::Lerp(const fquat_t& a, const fquat_t& b, float t) noexcept
    {
        t = xmath::Clamp(t, 0.0f, 1.0f);
        float dot = Dot(a, b);
        fquat_t<V> end = (dot >= 0.0f) ? b : -b;
        fquat_t<V> result;
        if constexpr (V)
        {
            __m128 t_vec = _mm_set1_ps(t);
            __m128 one_minus_t = _mm_sub_ps(_mm_set1_ps(1.0f), t_vec);
            result.m_XYZW = _mm_add_ps(_mm_mul_ps(a.m_XYZW, one_minus_t), _mm_mul_ps(end.m_XYZW, t_vec));
        }
        else
        {
            float omt = 1.0f - t;
            result.m_X = a.m_X * omt + end.m_X * t;
            result.m_Y = a.m_Y * omt + end.m_Y * t;
            result.m_Z = a.m_Z * omt + end.m_Z * t;
            result.m_W = a.m_W * omt + end.m_W * t;
        }
        return result.Normalize();
    }

    //------------------------------------------------------------------------------
    // LerpFast
    //------------------------------------------------------------------------------
    //
    // Linear interpolation between two quaternions.
    //
    // Params:
    //  a - Start quaternion.
    //  b - End quaternion.
    //  t - Interpolation factor (clamped 0-1).
    //
    // Returns:
    //  Normalized lerp result.
    //
    // Notes:
    //  Uses dot to choose shortest path; normalizes result.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::LerpFast(const fquat_t& Start, const fquat_t& End, float T ) noexcept
    {
        float LenSquared;
        float OneOverL;
        fquat_t Temp;
        float x0, y0, z0, w0;

        // Determine if quats are further than 90 degrees
        const float dot = Start.Dot(End); 

        // If dot is negative flip one of the quaterions
        if (dot < 0.0f) Temp = -Start;
        else            Temp = Start;

        // Compute interpolated values
        Temp += T * (End - Temp);

        // get squared length of new quaternion
        LenSquared = Temp.Dot(Temp);

        // Use home-baked polynomial to compute 1/sqrt(LenSquared)
        // Input range is 0.5 <-> 1.0
        // Ouput range is 1.414213 <-> 1.0

        if (LenSquared < 0.857f)
            OneOverL = (((0.699368f) * LenSquared) + -1.819985f) * LenSquared + 2.126369f;    //0.0000792
        else
            OneOverL = (((0.454012f) * LenSquared) + -1.403517f) * LenSquared + 1.949542f;    //0.0000373

        // Renormalize and return quaternion
        return Temp * OneOverL;
    }

    //------------------------------------------------------------------------------
    // Slerp
    //------------------------------------------------------------------------------
    //
    // Spherical linear interpolation.
    //
    // Params:
    //  a - Start.
    //  b - End.
    //  t - Factor (clamped 0-1).
    //
    // Returns:
    //  Slerp result.
    //
    // Notes:
    //  Handles shortest path; uses acos, sin for angle.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::Slerp(const fquat_t& a, const fquat_t& b, float t) noexcept
    {
        t = xmath::Clamp(t, 0.0f, 1.0f);
        float cos_omega = Dot(a, b);
        fquat_t<V> end = (cos_omega >= 0.0f) ? b : -b;
        cos_omega = xmath::Abs(cos_omega);
        if (cos_omega > 0.999f) return Lerp(a, end, t);
        radian omega    = xmath::Acos(cos_omega);
        float sin_omega = xmath::Sin(omega);
        float s1        = xmath::Sin((1.0f - t) * omega) / sin_omega;
        float s2        = xmath::Sin(t * omega) / sin_omega;
        if constexpr (V)
        {
            __m128 s1_vec = _mm_set1_ps(s1);
            __m128 s2_vec = _mm_set1_ps(s2);
            return fquat_t<V>(_mm_add_ps(_mm_mul_ps(a.m_XYZW, s1_vec), _mm_mul_ps(end.m_XYZW, s2_vec)));
        }
        else
        {
            return fquat_t<V>(a.m_X * s1 + end.m_X * s2, a.m_Y * s1 + end.m_Y * s2, a.m_Z * s1 + end.m_Z * s2, a.m_W * s1 + end.m_W * s2);
        }
    }

    //------------------------------------------------------------------------------
    // BlendAccurate (adapted as SlerpAccurate for precision)
    //------------------------------------------------------------------------------
    //
    // Computes accurate spherical linear interpolation between two quaternions.
    //
    // Params:
    //  a - Starting quaternion (normalized).
    //  b - Ending quaternion (normalized).
    //  t - Interpolation factor (0 to 1).
    //
    // Returns:
    //  Interpolated quaternion on the shortest arc.
    //
    // Notes:
    //  Adapted from old BlendAccurate; uses acos/sin for theta.
    //  Falls back to linear interp for small angles to avoid instability.
    //  Handles shortest path via sign flip if dot < 0.
    //  Assumes normalized inputs; asserts finite in debug.
    //  Right-hand rule; SIMD optimized where possible.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::SlerpAccurate(const fquat_t<V>& a, const fquat_t<V>& b, float t) noexcept
    {
        t = xmath::Clamp(t, 0.0f, 1.0f);
        float cs = Dot(a, b);
        bool flip = false;
        if (cs < 0.0f)
        {
            cs = -cs;
            flip = true;
        }
        float inv_t;
        if ((1.0f - cs) < 0.000001f)
        {
            inv_t = 1.0f - t;
        }
        else
        {
            radian theta = xmath::Acos(cs);
            float s = xmath::Sin(theta);
            inv_t = xmath::Sin((1.0f - t) * theta) / s;
            t = xmath::Sin(t * theta) / s;
        }

        if (flip)
        {
            t = -t;
        }
        fquat_t<V> result;
        if constexpr (V)
        {
            __m128 inv_t_vec = _mm_set1_ps(inv_t);
            __m128 t_vec = _mm_set1_ps(t);
            result.m_XYZW = _mm_add_ps(_mm_mul_ps(a.m_XYZW, inv_t_vec), _mm_mul_ps(b.m_XYZW, t_vec));
        }
        else
        {
            result.m_X = a.m_X * inv_t + b.m_X * t;
            result.m_Y = a.m_Y * inv_t + b.m_Y * t;
            result.m_Z = a.m_Z * inv_t + b.m_Z * t;
            result.m_W = a.m_W * inv_t + b.m_W * t;
        }
        return result;
    }


    //------------------------------------------------------------------------------
    // AngleBetween
    //------------------------------------------------------------------------------
    //
    // Angle between two quaternions in radians.
    //
    // Params:
    //  a - First.
    //  b - Second.
    //
    // Returns:
    //  2 * acos(dot) , between 0 and pi.
    //
    template <bool V >
    inline radian fquat_t<V>::AngleBetween(const fquat_t& a, const fquat_t& b) noexcept
    {
        float dot = xmath::Clamp(Dot(a, b), -1.0f, 1.0f);
        return 2.0f * xmath::Acos(dot);
    }

    //------------------------------------------------------------------------------
    // FromToRotation
    //------------------------------------------------------------------------------
    //
    // Rotation from from to to.
    //
    // Params:
    //  from - Normalized from vector.
    //  to - Normalized to vector.
    //
    // Returns:
    //  Quaternion aligning from to to.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::FromToRotation(const fvec3& from, const fvec3& to, const fvec3& up) noexcept
    {
        assert(from.isFinite() && to.isFinite() && up.isFinite());
        assert(from.isNormalized(1e-6f) && to.isNormalized(1e-6f) && up.isNormalized(1e-6f));
        const float d = from.Dot(to);
        if (d >= 1.0f - 1e-6f) 
        {
            return fromIdentity();
        }
        if (d <= -1.0f + 1e-6f) 
        {
            fvec3 axis = from.Cross(up).NormalizeSafeCopy();
            if (axis.isNearlyZero(1e-6f)) 
            {
                axis = from.Perpendicular(fvec3::fromRight()).NormalizeSafeCopy();
            }
            return fquat_t<V>(axis, radian(xmath::pi_v));
        }
        const fvec3 cross = from.Cross(to);
        const float s = xmath::Sqrt((1.0f + d) * 2.0f);
        const float inv_s = 1.0f / s;
        return fquat_t<V>(cross.m_X * inv_s, cross.m_Y * inv_s, cross.m_Z * inv_s, s * 0.5f);
    }


    //------------------------------------------------------------------------------
    // LookRotation
    //------------------------------------------------------------------------------
    //
    // Rotation looking in forward direction with up.
    //
    // Params:
    //  forward - Forward.
    //  up - Up.
    //
    // Returns:
    //  Quaternion for look.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::LookRotation(const fvec3& forward, const fvec3& up) noexcept
    {
        assert(forward.isFinite() && up.isFinite());
        assert(!forward.isNearlyZero(1e-6f) && !up.isNearlyZero(1e-6f));
        fvec3 f = forward.NormalizeSafeCopy();
        fvec3 u = up.NormalizeSafeCopy();
        fvec3 r = u.Cross(f).NormalizeSafeCopy();
        u = f.Cross(r);
        float trace = r.m_X + u.m_Y + f.m_Z;
        if (trace > 0.0f) 
        {
            float s = 0.5f / xmath::Sqrt(trace + 1.0f);
            return fquat_t<V>( (u.m_Z - f.m_Y) * s, (f.m_X - r.m_Z) * s, (r.m_Y - u.m_X) * s, 0.25f / s );
        }
        else if (r.m_X > u.m_Y && r.m_X > f.m_Z) 
        {
            float s = 2.0f * xmath::Sqrt(1.0f + r.m_X - u.m_Y - f.m_Z);
            return fquat_t<V>(0.25f * s, (r.m_Y + u.m_X) / s, (r.m_Z + f.m_X) / s, (u.m_Z - f.m_Y) / s);
        }
        else if (u.m_Y > f.m_Z) 
        {
            float s = 2.0f * xmath::Sqrt(1.0f + u.m_Y - r.m_X - f.m_Z);
            return fquat_t<V>((r.m_Y + u.m_X) / s, 0.25f * s, (u.m_Z + f.m_Y) / s, (f.m_X - r.m_Z) / s);
        }
        else 
        {
            float s = 2.0f * xmath::Sqrt(1.0f + f.m_Z - r.m_X - u.m_Y);
            return fquat_t<V>((r.m_Z + f.m_X) / s, (u.m_Z + f.m_Y) / s, 0.25f * s, (r.m_Y - u.m_X) / s);
        }
    }

    //------------------------------------------------------------------------------
    // RandomUnitQuaternion
    //------------------------------------------------------------------------------
    //
    // Generates a uniform random unit quaternion.
    //
    // Returns:
    //  Random rotation quaternion.
    //
    // Notes:
    //  Uses std::random_device and uniform distribution.
    //  Formula for uniform on S3.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::RandomUnitQuaternion(void) noexcept
    {
        thread_local std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        float u1 = dist(gen);
        float u2 = dist(gen);
        float u3 = dist(gen);
        float s1 = xmath::Sqrt(1.0f - u1);
        float s2 = xmath::Sqrt(u1);
        radian theta1 = 2.0f * xmath::pi_v * u2;
        radian theta2 = 2.0f * xmath::pi_v * u3;
        float x = s1 * xmath::Sin(theta1);
        float y = s1 * xmath::Cos(theta1);
        float z = s2 * xmath::Sin(theta2);
        float w = s2 * xmath::Cos(theta2);
        return fquat_t<V>(x, y, z, w);
    }

    //------------------------------------------------------------------------------
    // operator std::string
    //------------------------------------------------------------------------------
    //
    // Converts quaternion to string representation.
    //
    // Returns:
    //  String in format "(x, y, z, w)".
    //
    // Notes:
    //  Delegates to ToString for consistency.
    //
    template <bool V>
    inline fquat_t<V>::operator std::string() const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // ToString
    //------------------------------------------------------------------------------
    //
    // Returns string representation of the quaternion.
    //
    // Returns:
    //  String in format "(x, y, z, w)".
    //
    // Notes:
    //  Uses std::to_string for each component; assumes finite.
    //
    template <bool V>
    std::string fquat_t<V>::ToString() const noexcept
    {
        return std::format("( {}, {}, {}, {} )", this->m_X, this->m_Y, this->m_Z, this->m_W);
    }

    //------------------------------------------------------------------------------
    // operator<<
    //------------------------------------------------------------------------------
    //
    // Streams quaternion to output stream.
    //
    // Params:
    //  os - Output stream.
    //  quat - Quaternion to stream.
    //
    // Returns:
    //  Reference to the output stream.
    //
    // Notes:
    //  Uses ToString for formatting.
    //
    template <bool V>
    inline std::ostream& operator<<(std::ostream& os, const fquat_t<V>& quat) noexcept
    {
        return os << quat.ToString();
    }

    //------------------------------------------------------------------------------
    // Dot (instance)
    //------------------------------------------------------------------------------
    //
    // Computes dot product with another quaternion (static method as member).
    //
    // Params:
    //  other - The second quaternion.
    //
    // Returns:
    //  x1*x2 + y1*y2 + z1*z2 + w1*w2.
    //
    // Notes:
    //  Delegates to static Dot; assumes finite inputs.
    //
    template <bool V>
    inline float fquat_t<V>::Dot(const fquat_t& other) const noexcept
    {
        return Dot(*this, other);
    }

    //------------------------------------------------------------------------------
    // Lerp (instance)
    //------------------------------------------------------------------------------
    //
    // Linear interpolation with another quaternion (static method as member).
    //
    // Params:
    //  other - End quaternion.
    //  t - Interpolation factor (clamped 0-1).
    //
    // Returns:
    //  Normalized lerp result.
    //
    // Notes:
    //  Delegates to static Lerp; ensures shortest path.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::Lerp(const fquat_t& other, float t) const noexcept
    {
        return Lerp(*this, other, t);
    }

    //------------------------------------------------------------------------------

    template <bool V>
    inline fquat_t<V> fquat_t<V>::LerpFast(const fquat_t& other, float t) const noexcept
    {
        return LerpFast(*this, other, t);
    }

    // In xmath_fquat_inline.h

    //------------------------------------------------------------------------------
    // LerpUnclamped
    //------------------------------------------------------------------------------
    //
    // Linear interpolation between two quaternions without clamping t.
    //
    // Params:
    //  a - Start quaternion.
    //  b - End quaternion.
    //  t - Interpolation factor (can be outside [0,1] for extrapolation).
    //
    // Returns:
    //  Normalized lerp result, shortest path.
    //
    // Notes:
    //  Similar to Lerp but no t clamp; flips b if dot < 0.
    //  Asserts finite inputs in debug.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::LerpUnclamped(const fquat_t<V>& a, const fquat_t<V>& b, float t) noexcept
    {
        assert(a.isFinite() && b.isFinite() && xmath::isFinite(t));
        float dot = Dot(a, b);
        fquat_t<V> end = (dot >= 0.0f) ? b : -b;
        fquat_t<V> result;
        float omt = 1.0f - t;
        if constexpr (V)
        {
            __m128 omt_vec = _mm_set1_ps(omt);
            __m128 t_vec = _mm_set1_ps(t);
            result.m_XYZW = _mm_add_ps(_mm_mul_ps(a.m_XYZW, omt_vec), _mm_mul_ps(end.m_XYZW, t_vec));
        }
        else
        {
            result.m_X = a.m_X * omt + end.m_X * t;
            result.m_Y = a.m_Y * omt + end.m_Y * t;
            result.m_Z = a.m_Z * omt + end.m_Z * t;
            result.m_W = a.m_W * omt + end.m_W * t;
        }
        result.Normalize();
        assert(result.isNormalized(1e-6f));
        return result;
    }

    //------------------------------------------------------------------------------
    // SlerpUnclamped
    //------------------------------------------------------------------------------
    //
    // Spherical linear interpolation without clamping t.
    //
    // Params:
    //  a - Start.
    //  b - End.
    //  t - Factor (can be outside [0,1] for extrapolation).
    //
    // Returns:
    //  Slerp result, shortest path.
    //
    // Notes:
    //  Similar to Slerp but no t clamp; fallback to lerp if cos > 0.999.
    //  Flips b if dot < 0; asserts finite in debug.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::SlerpUnclamped(const fquat_t<V>& a, const fquat_t<V>& b, float t) noexcept
    {
        assert(a.isFinite() && b.isFinite() && xmath::isFinite(t));
        float cos_omega = Dot(a, b);
        fquat_t<V> end = (cos_omega >= 0.0f) ? b : -b;
        cos_omega = xmath::Abs(cos_omega);
        if (cos_omega > 0.999f)
        {
            return LerpUnclamped(a, end, t);
        }
        float omega = std::acos(cos_omega);
        float sin_omega = std::sin(omega);
        float s1 = std::sin((1.0f - t) * omega) / sin_omega;
        float s2 = std::sin(t * omega) / sin_omega;
        fquat_t<V> result;
        if constexpr (V)
        {
            __m128 s1_vec = _mm_set1_ps(s1);
            __m128 s2_vec = _mm_set1_ps(s2);
            result.m_XYZW = _mm_add_ps(_mm_mul_ps(a.m_XYZW, s1_vec), _mm_mul_ps(end.m_XYZW, s2_vec));
        }
        else
        {
            result.m_X = a.m_X * s1 + end.m_X * s2;
            result.m_Y = a.m_Y * s1 + end.m_Y * s2;
            result.m_Z = a.m_Z * s1 + end.m_Z * s2;
            result.m_W = a.m_W * s1 + end.m_W * s2;
        }
        assert(result.isFinite());
        return result;
    }

    //------------------------------------------------------------------------------
    // LerpUnclamped (instance)
    //------------------------------------------------------------------------------
    //
    // Instance LerpUnclamped with other.
    //
    // Params:
    //  other - End.
    //  t - Factor.
    //
    // Returns:
    //  LerpUnclamped result.
    //
    // Notes:
    //  Delegates to static LerpUnclamped.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::LerpUnclamped(const fquat_t<V>& other, float t) const noexcept
    {
        return LerpUnclamped(*this, other, t);
    }

    //------------------------------------------------------------------------------
    // SlerpUnclamped (instance)
    //------------------------------------------------------------------------------
    //
    // Instance SlerpUnclamped with other.
    //
    // Params:
    //  other - End.
    //  t - Factor.
    //
    // Returns:
    //  SlerpUnclamped result.
    //
    // Notes:
    //  Delegates to static SlerpUnclamped.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::SlerpUnclamped(const fquat_t<V>& other, float t) const noexcept
    {
        return SlerpUnclamped(*this, other, t);
    }

    //------------------------------------------------------------------------------
    // Squad
    //------------------------------------------------------------------------------
    //
    // Spherical cubic interpolation between two quaternions with tangents.
    //
    // Params:
    //  a - Start quaternion.
    //  a_tangent - Tangent at start.
    //  b - End quaternion.
    //  b_tangent - Tangent at end.
    //  t - Interpolation factor (0 to 1).
    //
    // Returns:
    //  Smoothly interpolated quaternion.
    //
    // Notes:
    //  Uses Slerp for interpolation; assumes normalized inputs.
    //  Asserts finite inputs in debug; right-hand rule.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::Squad(const fquat_t<V>& a, const fquat_t<V>& a_tangent, const fquat_t<V>& b, const fquat_t<V>& b_tangent, float t) noexcept
    {
        assert(a.isFinite() && a_tangent.isFinite() && b.isFinite() && b_tangent.isFinite() && xmath::isFinite(t));
        t = xmath::Clamp(t, 0.0f, 1.0f);
        fquat_t<V> q1 = Slerp(a, a_tangent, t);
        fquat_t<V> q2 = Slerp(b, b_tangent, t);
        fquat_t<V> result = Slerp(q1, q2, 2.0f * t * (1.0f - t));
        assert(result.isFinite());
        return result;
    }

    //------------------------------------------------------------------------------
    // Squad (instance)
    //------------------------------------------------------------------------------
    //
    // Instance Squad with other quaternions and tangents.
    //
    // Params:
    //  a_tangent - Tangent at this quaternion.
    //  b - End quaternion.
    //  b_tangent - Tangent at end.
    //  t - Interpolation factor.
    //
    // Returns:
    //  Squad result.
    //
    // Notes:
    //  Delegates to static Squad.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::Squad(const fquat_t<V>& a_tangent, const fquat_t<V>& b, const fquat_t<V>& b_tangent, float t) const noexcept
    {
        return Squad(*this, a_tangent, b, b_tangent, t);
    }

    //------------------------------------------------------------------------------
    // setupFromToRotation
    //------------------------------------------------------------------------------
    //
    // Sets quaternion to rotation from one vector to another.
    //
    // Params:
    //  from - Start direction (normalized).
    //  to - Target direction (normalized).
    //  up - Up direction for degenerate cases (default Up).
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Delegates to FromToRotation; asserts finite and normalized inputs in debug.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::setupFromToRotation(const fvec3& from, const fvec3& to, const fvec3& up) noexcept
    {
        assert(from.isFinite() && to.isFinite() && up.isFinite());
        assert(from.isNormalized(1e-6f) && to.isNormalized(1e-6f) && up.isNormalized(1e-6f));
        *this = FromToRotation(from, to, up);
        assert(isFinite());
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupLookRotation
    //------------------------------------------------------------------------------
    //
    // Sets quaternion to look in forward direction with up.
    //
    // Params:
    //  forward - Forward direction.
    //  up - Up direction.
    //
    // Returns:
    //  Reference (chainable).
    //
    // Notes:
    //  Delegates to LookRotation; asserts finite and non-zero inputs in debug.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::setupLookRotation(const fvec3& forward, const fvec3& up) noexcept
    {
        assert(forward.isFinite() && up.isFinite());
        assert(!forward.isNearlyZero(1e-6f) && !up.isNearlyZero(1e-6f));
        *this = LookRotation(forward, up);
        assert(isFinite());
        return *this;
    }

    //------------------------------------------------------------------------------
    // isNearlyZero
    //------------------------------------------------------------------------------
    //
    // Checks if quaternion length is near zero.
    //
    // Params:
    //  tolerance - Allowed deviation (default 1e-6f).
    //
    // Returns:
    //  True if length^2 < tolerance.
    //
    // Notes:
    //  Uses LengthSq; asserts finite in debug.
    //
    template <bool V>
    inline bool fquat_t<V>::isNearlyZero(float tolerance) const noexcept
    {
        assert(isFinite());
        return LengthSq() < tolerance;
    }

    //------------------------------------------------------------------------------
    // Slerp (instance)
    //------------------------------------------------------------------------------
    //
    // Spherical linear interpolation with another quaternion (static method as member).
    //
    // Params:
    //  other - End quaternion.
    //  t - Interpolation factor (clamped 0-1).
    //
    // Returns:
    //  Slerp result.
    //
    // Notes:
    //  Delegates to static Slerp; handles shortest path.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::Slerp(const fquat_t& other, float t) const noexcept
    {
        return Slerp(*this, other, t);
    }

    //------------------------------------------------------------------------------
    // SlerpAccurate (instance)
    //------------------------------------------------------------------------------
    //
    // Computes accurate spherical linear interpolation between two quaternions.
    //
    // Params:
    //  other - Ending quaternion (normalized).
    //  t     - Interpolation factor (0 to 1).
    //
    // Returns:
    //  Interpolated quaternion on the shortest arc.
    //
    // Notes:
    //  Adapted from old BlendAccurate; uses acos/sin for theta.
    //  Falls back to linear interp for small angles to avoid instability.
    //  Handles shortest path via sign flip if dot < 0.
    //  Assumes normalized inputs; asserts finite in debug.
    //  Right-hand rule; SIMD optimized where possible.
    //
    template <bool V >
    inline fquat_t<V> fquat_t<V>::SlerpAccurate(const fquat_t& other, float t) const noexcept
    {
        return SlerpAccurate(*this, other, t);
    }

    //------------------------------------------------------------------------------
    // AngleBetween (instance)
    //------------------------------------------------------------------------------
    //
    // Angle between this and another quaternion (static method as member).
    //
    // Params:
    //  other - The other quaternion.
    //
    // Returns:
    //  2 * acos(dot), between 0 and pi.
    //
    // Notes:
    //  Delegates to static AngleBetween.
    //
    template <bool V>
    inline radian fquat_t<V>::AngleBetween(const fquat_t& other) const noexcept
    {
        return AngleBetween(*this, other);
    }

    //------------------------------------------------------------------------------
    // Length
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean length (magnitude) of the quaternion.
    //
    // Returns:
    //  sqrt(x^2 + y^2 + z^2 + w^2).
    //
    // Notes:
    //  Uses Dot and xmath::Sqrt; asserts finite in debug.
    //
    template <bool V>
    inline float fquat_t<V>::Length() const noexcept
    {
        return xmath::Sqrt(LengthSq());
    }

    //------------------------------------------------------------------------------
    // LengthSq
    //------------------------------------------------------------------------------
    //
    // Computes the squared Euclidean length of the quaternion.
    //
    // Returns:
    //  x^2 + y^2 + z^2 + w^2.
    //
    // Notes:
    //  Uses Dot; faster than Length for comparisons.
    //
    template <bool V>
    inline float fquat_t<V>::LengthSq() const noexcept
    {
        return Dot(*this, *this);
    }

    //------------------------------------------------------------------------------
    // NormalizeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy of the quaternion.
    //
    // Returns:
    //  Quaternion / Length, or fromIdentity if length < 1e-6.
    //
    // Notes:
    //  Uses Length; asserts finite in debug.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::NormalizeCopy() const noexcept
    {
        float len = Length();
        if (len < 1e-6f) return fromIdentity();
        if constexpr (V)
        {
            return fquat_t<V>(_mm_div_ps(this->m_XYZW, _mm_set1_ps(len)));
        }
        else
        {
            return fquat_t<V>(this->m_X / len, this->m_Y / len, this->m_Z / len, this->m_W / len);
        }
    }

    //------------------------------------------------------------------------------
    // Normalize
    //------------------------------------------------------------------------------
    //
    // Normalizes the quaternion in-place to unit length.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses NormalizeCopy; sets to fromIdentity if length < 1e-6.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::Normalize() noexcept
    {
        *this = NormalizeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // NormalizeSafeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy or fromIdentity if length is near zero.
    //
    // Returns:
    //  Normalized quaternion or fromIdentity if length < 1e-6.
    //
    // Notes:
    //  Avoids division by zero; uses Length.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::NormalizeSafeCopy() const noexcept
    {
        float len = Length();
        if (len < 1e-6f) return fromIdentity();
        if constexpr (V)
        {
            return fquat_t<V>(_mm_div_ps(this->m_XYZW, _mm_set1_ps(len)));
        }
        else
        {
            return fquat_t<V>(this->m_X / len, this->m_Y / len, this->m_Z / len, this->m_W / len);
        }
    }

    //------------------------------------------------------------------------------
    // NormalizeSafe
    //------------------------------------------------------------------------------
    //
    // Normalizes the quaternion in-place or sets to fromIdentity if length is near zero.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses NormalizeSafeCopy; avoids division by zero.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::NormalizeSafe() noexcept
    {
        *this = NormalizeSafeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all quaternion components are finite (no NaN or infinity).
    //
    // Returns:
    //  True if all components are finite.
    //
    // Notes:
    //  Uses xmath::isFinite on x, y, z, w components.
    //
    template <bool V>
    inline bool fquat_t<V>::isFinite() const noexcept
    {
        return xmath::isFinite(this->m_X) && xmath::isFinite(this->m_Y) &&
               xmath::isFinite(this->m_Z) && xmath::isFinite(this->m_W);
    }

    //------------------------------------------------------------------------------
    // isNormalized
    //------------------------------------------------------------------------------
    //
    // Checks if the quaternion is normalized (length ~ 1).
    //
    // Params:
    //  tolerance - Allowed deviation from 1 (default 1e-6).
    //
    // Returns:
    //  True if |length^2 - 1| < tolerance.
    //
    template <bool V>
    inline bool fquat_t<V>::isNormalized(float tolerance) const noexcept
    {
        float len_sq = LengthSq();
        return xmath::Abs(len_sq - 1.0f) < tolerance;
    }

    //------------------------------------------------------------------------------
    // isNearlyIdentity
    //------------------------------------------------------------------------------
    //
    // Checks if the quaternion is nearly identity (0,0,0,1).
    //
    // Params:
    //  tolerance - Allowed deviation (default 1e-6).
    //
    // Returns:
    //  True if x,y,z < tolerance and |w - 1| < tolerance.
    //
    template <bool V>
    inline bool fquat_t<V>::isNearlyIdentity(float tolerance) const noexcept
    {
        return xmath::Abs(this->m_X) < tolerance && xmath::Abs(this->m_Y) < tolerance &&
               xmath::Abs(this->m_Z) < tolerance && xmath::Abs(this->m_W - 1.0f) < tolerance;
    }

    //------------------------------------------------------------------------------
    // Equals
    //------------------------------------------------------------------------------
    //
    // Checks if two quaternions are equal within a tolerance.
    //
    // Params:
    //  other - Quaternion to compare.
    //  tolerance - Allowed deviation (default 1e-6).
    //
    // Returns:
    //  True if all components differ by < tolerance.
    //
    template <bool V>
    inline bool fquat_t<V>::Equals(const fquat_t& other, float tolerance) const noexcept
    {
        return xmath::Abs(this->m_X - other.m_X) < tolerance &&
               xmath::Abs(this->m_Y - other.m_Y) < tolerance &&
               xmath::Abs(this->m_Z - other.m_Z) < tolerance &&
               xmath::Abs(this->m_W - other.m_W) < tolerance;
    }

    //------------------------------------------------------------------------------
    // ConjugateCopy
    //------------------------------------------------------------------------------
    //
    // Returns the conjugate of the quaternion.
    //
    // Returns:
    //  (-x, -y, -z, w).
    //
    // Notes:
    //  Negates imaginary components; assumes finite.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::ConjugateCopy() const noexcept
    {
        if constexpr (V)
        {
            __m128 neg = _mm_set_ps(1.0f, -1.0f, -1.0f, -1.0f);
            return fquat_t<V>(_mm_mul_ps(this->m_XYZW, neg));
        }
        else
        {
            return fquat_t<V>(-this->m_X, -this->m_Y, -this->m_Z, this->m_W);
        }
    }

    //------------------------------------------------------------------------------
    // Conjugate
    //------------------------------------------------------------------------------
    //
    // Conjugates the quaternion in-place.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses ConjugateCopy; negates x,y,z.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::Conjugate() noexcept
    {
        *this = ConjugateCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // InverseCopy
    //------------------------------------------------------------------------------
    //
    // Returns the inverse of the quaternion.
    //
    // Returns:
    //  Conjugate / Length^2, or fromIdentity if length < 1e-6.
    //
    // Notes:
    //  Uses ConjugateCopy and LengthSq; assumes finite.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::InverseCopy() const noexcept
    {
        float len_sq = LengthSq();
        if (len_sq < 1e-6f) return fromIdentity();
        fquat_t<V> conj = ConjugateCopy();
        if constexpr (V)
        {
            return fquat_t<V>(_mm_div_ps(conj.m_XYZW, _mm_set1_ps(len_sq)));
        }
        else
        {
            return fquat_t<V>(conj.m_X / len_sq, conj.m_Y / len_sq, conj.m_Z / len_sq, conj.m_W / len_sq);
        }
    }

    //------------------------------------------------------------------------------
    // Inverse
    //------------------------------------------------------------------------------
    //
    // Inverts the quaternion in-place.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses InverseCopy; sets to fromIdentity if length < 1e-6.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::Inverse() noexcept
    {
        *this = InverseCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Axis
    //------------------------------------------------------------------------------
    //
    // Returns the rotation axis of the quaternion.
    //
    // Returns:
    //  Normalized axis (x,y,z)/sin(theta/2) or Right if angle is small.
    //
    // Notes:
    //  Assumes normalized quaternion; uses Right for near-zero angle.
    //
    template <bool V>
    inline fvec3 fquat_t<V>::Axis() const noexcept
    {
        float sin_theta_half = xmath::Sqrt(this->m_X * this->m_X + this->m_Y * this->m_Y + this->m_Z * this->m_Z);
        if (sin_theta_half < 1e-6f) return fvec3::fromRight();
        float inv_sin = 1.0f / sin_theta_half;
        return fvec3(this->m_X * inv_sin, this->m_Y * inv_sin, this->m_Z * inv_sin);
    }

    //------------------------------------------------------------------------------
    // Angle
    //------------------------------------------------------------------------------
    //
    // Returns the rotation angle of the quaternion in radians.
    //
    // Returns:
    //  2 * acos(w).
    //
    // Notes:
    //  Assumes normalized quaternion; clamps w to [-1, 1].
    //
    template <bool V>
    inline radian fquat_t<V>::Angle() const noexcept
    {
        return radian(2.0f * std::acos(xmath::Clamp(this->m_W, -1.0f, 1.0f)));
    }

    //------------------------------------------------------------------------------
    // ToEuler
    //------------------------------------------------------------------------------
    //
    // Converts quaternion to Euler angles (YXZ order).
    //
    // Returns:
    //  radian3(pitch, yaw, roll).
    //
    // Notes:
    //  Uses YXZ order for right-hand system; handles singularities with atan2.
    //
    template <bool V>
    inline radian3 fquat_t<V>::ToEuler() const noexcept
    {
        assert(isFinite());
        float x = this->m_X, y = this->m_Y, z = this->m_Z, w = this->m_W;
        float siny_cosp = 2 * (w * y + z * x);
        float cosy_cosp = 1 - 2 * (y * y + x * x);
        radian yaw = radian(xmath::Atan2(siny_cosp, cosy_cosp));

        float sinp = 2 * (w * x - y * z);
        radian pitch;
        if (std::abs(sinp) >= 1.0f)
        {
            pitch = radian(xmath::CopySign(xmath::pi_over2_v.m_Value, sinp)); // Handle singularity
        }
        else
        {
            pitch = radian(std::asin(sinp));
        }

        float sinr_cosp = 2 * (w * z + y * x);
        float cosr_cosp = 1 - 2 * (x * x + z * z);
        radian roll = radian(xmath::Atan2(sinr_cosp, cosr_cosp));

        return radian3(pitch, yaw, roll);
    }
/*
    template <bool V>
    inline radian3 fquat_t<V>::ToEuler() const noexcept
    {
        float x = this->m_X, y = this->m_Y, z = this->m_Z, w = this->m_W;
        float sinr_cosp = 2.0f * (w * x + y * z);
        float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
        radian roll = xmath::Atan2(sinr_cosp, cosr_cosp);

        float sinp = 2.0f * (w * y - z * x);
        radian pitch;
        if (std::abs(sinp) >= 1.0f) 
        {
            pitch = xmath::radian( xmath::CopySign(xmath::pi_over2_v.m_Value, sinp)); // Handle singularity
        }
        else 
        {
            pitch = xmath::Asin(sinp);
        }

        float siny_cosp = 2.0f * (w * z + x * y);
        float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
        radian yaw = xmath::Atan2(siny_cosp, cosy_cosp);

        return radian3(pitch, yaw, roll);
    }
    */

    //------------------------------------------------------------------------------
    // Forward
    //------------------------------------------------------------------------------
    //
    // Returns the forward direction after quaternion rotation.
    //
    // Returns:
    //  Rotated (0,0,1) vector.
    //
    // Notes:
    //  Uses quaternion rotation formula; assumes right-hand rule.
    //
    template <bool V>
    inline fvec3 fquat_t<V>::Forward() const noexcept
    {
        return RotateVector(fvec3::fromForward());
    }

    //------------------------------------------------------------------------------
    // Up
    //------------------------------------------------------------------------------
    //
    // Returns the up direction after quaternion rotation.
    //
    // Returns:
    //  Rotated (0,1,0) vector.
    //
    // Notes:
    //  Uses RotateVector; assumes right-hand rule.
    //
    template <bool V>
    inline fvec3 fquat_t<V>::Up() const noexcept
    {
        return RotateVector(fvec3::fromUp());
    }

    //------------------------------------------------------------------------------
    // Right
    //------------------------------------------------------------------------------
    //
    // Returns the right direction after quaternion rotation.
    //
    // Returns:
    //  Rotated (1,0,0) vector.
    //
    // Notes:
    //  Uses RotateVector; assumes right-hand rule.
    //
    template <bool V>
    inline fvec3 fquat_t<V>::Right() const noexcept
    {
        return RotateVector(fvec3::fromRight());
    }

    //------------------------------------------------------------------------------
    // ToAxisAngle
    //------------------------------------------------------------------------------
    //
    // Converts quaternion to axis-angle representation.
    //
    // Returns:
    //  Pair of axis (normalized) and angle in radians.
    //
    // Notes:
    //  Uses Axis and Angle; returns Right and 0 if angle is small.
    //
    template <bool V>
    inline std::pair<fvec3, radian> fquat_t<V>::ToAxisAngle() const noexcept
    {
        float sin_theta_half = xmath::Sqrt(this->m_X * this->m_X + this->m_Y * this->m_Y + this->m_Z * this->m_Z);
        if (sin_theta_half < 1e-6f) return { fvec3::fromRight(), radian(0.0f) };
        return { Axis(), Angle() };
    }

    //------------------------------------------------------------------------------
    // Delta
    //------------------------------------------------------------------------------
    //
    // Returns the quaternion representing the rotation from this to other.
    //
    // Params:
    //  other - Target quaternion.
    //
    // Returns:
    //  other * inverse(this).
    //
    // Notes:
    //  Uses InverseCopy and quaternion multiplication.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::Delta(const fquat_t& other) const noexcept
    {
        return other * InverseCopy();
    }

    //------------------------------------------------------------------------------
    // LogCopy
    //------------------------------------------------------------------------------
    //
    // Returns the logarithm of the quaternion.
    //
    // Returns:
    //  log(w + sqrt(w^2 + x^2 + y^2 + z^2)) + (x,y,z)/|v| * acos(w/|q|).
    //
    // Notes:
    //  For unit quaternions, simplifies to (0, axis * angle/2).
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::LogCopy() const noexcept
    {
        float len = Length();
        if (len < 1e-6f) return fromZero();
        float w = this->m_W / len;
        radian angle = xmath::Acos(xmath::Clamp(w, -1.0f, 1.0f));
        float sin_angle = xmath::Sin(angle);
        if (sin_angle < 1e-6f) return fquat_t<V>(0.0f, 0.0f, 0.0f, std::log(len));
        radian coeff = angle / sin_angle;
        return fquat_t<V>(this->m_X * coeff.m_Value / len, this->m_Y * coeff.m_Value / len, this->m_Z * coeff.m_Value / len, xmath::Log(len));
    }

    //------------------------------------------------------------------------------
    // Log
    //------------------------------------------------------------------------------
    //
    // Computes the logarithm of the quaternion in-place.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses LogCopy; assumes finite.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::Log() noexcept
    {
        *this = LogCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // ExpCopy
    //------------------------------------------------------------------------------
    //
    // Returns the exponential of the quaternion.
    //
    // Returns:
    //  exp(w) * (cos(|v|) + v/|v| * sin(|v|)).
    //
    // Notes:
    //  For pure quaternions (w=0), gives rotation quaternion.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::ExpCopy() const noexcept
    {
        radian v_len = radian( xmath::Sqrt(this->m_X * this->m_X + this->m_Y * this->m_Y + this->m_Z * this->m_Z) );
        float exp_w = xmath::Exp(this->m_W);
        if (v_len.m_Value < 1e-6f) return fquat_t<V>(0.0f, 0.0f, 0.0f, exp_w);
        float s = xmath::Sin(v_len) / v_len.m_Value;
        return fquat_t<V>(exp_w * s * this->m_X, exp_w * s * this->m_Y, exp_w * s * this->m_Z, exp_w * xmath::Cos(v_len));
    }

    //------------------------------------------------------------------------------
    // Exp
    //------------------------------------------------------------------------------
    //
    // Computes the exponential of the quaternion in-place.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses ExpCopy; assumes finite.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::Exp() noexcept
    {
        *this = ExpCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateXCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated by rx around the X-axis (local).
    //
    // Params:
    //  rx - Rotation angle in radians.
    //
    // Returns:
    //  this * quaternion(rx, X-axis).
    //
    // Notes:
    //  Post-multiplies; right-hand rule; asserts finite rx.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::RotateXCopy(radian rx) const noexcept
    {
        assert(xmath::isFinite(rx.m_Value));
        fquat_t<V> qx;
        qx.setupRotationX(rx);
        return *this * qx;
    }

    //------------------------------------------------------------------------------
    // RotateYCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated by ry around the Y-axis (local).
    //
    // Params:
    //  ry - Rotation angle in radians.
    //
    // Returns:
    //  this * quaternion(ry, Y-axis).
    //
    // Notes:
    //  Post-multiplies; right-hand rule; asserts finite ry.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::RotateYCopy(radian ry) const noexcept
    {
        assert(xmath::isFinite(ry.m_Value));
        fquat_t<V> qy;
        qy.setupRotationY(ry);
        return *this * qy;
    }

    //------------------------------------------------------------------------------
    // RotateZCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated by rz around the Z-axis (local).
    //
    // Params:
    //  rz - Rotation angle in radians.
    //
    // Returns:
    //  this * quaternion(rz, Z-axis).
    //
    // Notes:
    //  Post-multiplies; right-hand rule; asserts finite rz.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::RotateZCopy(radian rz) const noexcept
    {
        assert(xmath::isFinite(rz.m_Value));
        fquat_t<V> qz;
        qz.setupRotationZ(rz);
        return *this * qz;
    }

    //------------------------------------------------------------------------------
    // RotateX
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion by rx around the X-axis (local).
    //
    // Params:
    //  rx - Rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses RotateXCopy; post-multiplies.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::RotateX(radian rx) noexcept
    {
        *this = RotateXCopy(rx);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateY
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion by ry around the Y-axis (local).
    //
    // Params:
    //  ry - Rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses RotateYCopy; post-multiplies.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::RotateY(radian ry) noexcept
    {
        *this = RotateYCopy(ry);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateZ
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion by rz around the Z-axis (local).
    //
    // Params:
    //  rz - Rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses RotateZCopy; post-multiplies.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::RotateZ(radian rz) noexcept
    {
        *this = RotateZCopy(rz);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateX
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion by rx around the X-axis (world).
    //
    // Params:
    //  rx - Rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Pre-multiplies: quaternion(rx, X-axis) * this.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::PreRotateX(radian rx) noexcept
    {
        *this = PreRotateXCopy(rx);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateY
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion by ry around the Y-axis (world).
    //
    // Params:
    //  ry - Rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Pre-multiplies: quaternion(ry, Y-axis) * this.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::PreRotateY(radian ry) noexcept
    {
        *this = PreRotateYCopy(ry);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateZ
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion by rz around the Z-axis (world).
    //
    // Params:
    //  rz - Rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Pre-multiplies: quaternion(rz, Z-axis) * this.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::PreRotateZ(radian rz) noexcept
    {
        *this = PreRotateZCopy(rz);
        return *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateXCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated by rx around the X-axis (world).
    //
    // Params:
    //  rx - Rotation angle in radians.
    //
    // Returns:
    //  quaternion(rx, X-axis) * this.
    //
    // Notes:
    //  Pre-multiplies; right-hand rule; asserts finite rx.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::PreRotateXCopy(radian rx) const noexcept
    {
        assert(xmath::isFinite(rx.m_Value));
        fquat_t<V> qx;
        qx.setupRotationX(rx);
        return qx * *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateYCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated by ry around the Y-axis (world).
    //
    // Params:
    //  ry - Rotation angle in radians.
    //
    // Returns:
    //  quaternion(ry, Y-axis) * this.
    //
    // Notes:
    //  Pre-multiplies; right-hand rule; asserts finite ry.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::PreRotateYCopy(radian ry) const noexcept
    {
        assert(xmath::isFinite(ry.m_Value));
        fquat_t<V> qy;
        qy.setupRotationY(ry);
        return qy * *this;
    }

    //------------------------------------------------------------------------------
    // PreRotateZCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated by rz around the Z-axis (world).
    //
    // Params:
    //  rz - Rotation angle in radians.
    //
    // Returns:
    //  quaternion(rz, Z-axis) * this.
    //
    // Notes:
    //  Pre-multiplies; right-hand rule; asserts finite rz.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::PreRotateZCopy(radian rz) const noexcept
    {
        assert(xmath::isFinite(rz.m_Value));
        fquat_t<V> qz;
        qz.setupRotationZ(rz);
        return qz * *this;
    }

    //------------------------------------------------------------------------------
    // RotateTowardsCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy rotated towards the target quaternion by maxDelta.
    //
    // Params:
    //  target - Target quaternion.
    //  maxDelta - Maximum rotation angle in radians.
    //
    // Returns:
    //  Slerp result limited to maxDelta.
    //
    // Notes:
    //  Uses Slerp; clamps angle to maxDelta; assumes normalized inputs.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::RotateTowardsCopy(const fquat_t& target, radian maxDelta) const noexcept
    {
        assert(xmath::isFinite(maxDelta.m_Value) && maxDelta.m_Value >= 0.0f);
        float angle = AngleBetween(target).m_Value;
        if (angle < 1e-6f) return *this;
        float t = xmath::Min(maxDelta.m_Value / angle, 1.0f);
        return Slerp(target, t);
    }

    //------------------------------------------------------------------------------
    // RotateTowards
    //------------------------------------------------------------------------------
    //
    // Rotates the quaternion towards the target by maxDelta in-place.
    //
    // Params:
    //  target - Target quaternion.
    //  maxDelta - Maximum rotation angle in radians.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses RotateTowardsCopy.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::RotateTowards(const fquat_t& target, radian maxDelta) noexcept
    {
        *this = RotateTowardsCopy(target, maxDelta);
        return *this;
    }

    //------------------------------------------------------------------------------
    // RotateVector
    //------------------------------------------------------------------------------
    //
    // Rotates a vector by the quaternion.
    //
    // Params:
    //  v - Vector to rotate.
    //
    // Returns:
    //  Rotated vector q * v * q^-1.
    //
    // Notes:
    //  Uses optimized quaternion-vector rotation; assumes normalized quaternion.
    //
    template <bool V>
    inline fvec3 fquat_t<V>::RotateVector(const fvec3& v) const noexcept
    {
        fvec3 t = fvec3(this->m_X, this->m_Y, this->m_Z).Cross(v) * 2.0f;
        return v + t * this->m_W + fvec3(this->m_X, this->m_Y, this->m_Z).Cross(t);
    }

    //------------------------------------------------------------------------------
    // operator+
    //------------------------------------------------------------------------------
    //
    // Adds two quaternions component-wise.
    //
    // Params:
    //  other - Quaternion to add.
    //
    // Returns:
    //  (x1+x2, y1+y2, z1+z2, w1+w2).
    //
    // Notes:
    //  Not typically used for rotations; assumes finite.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::operator+(const fquat_t& other) const noexcept
    {
        if constexpr (V)
        {
            return fquat_t<V>(_mm_add_ps(this->m_XYZW, other.m_XYZW));
        }
        else
        {
            return fquat_t<V>(this->m_X + other.m_X, this->m_Y + other.m_Y, this->m_Z + other.m_Z, this->m_W + other.m_W);
        }
    }

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Subtracts two quaternions component-wise.
    //
    // Params:
    //  other - Quaternion to subtract.
    //
    // Returns:
    //  (x1-x2, y1-y2, z1-z2, w1-w2).
    //
    // Notes:
    //  Not typically used for rotations; assumes finite.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::operator-(const fquat_t& other) const noexcept
    {
        if constexpr (V)
        {
            return fquat_t<V>(_mm_sub_ps(this->m_XYZW, other.m_XYZW));
        }
        else
        {
            return fquat_t<V>(this->m_X - other.m_X, this->m_Y - other.m_Y, this->m_Z - other.m_Z, this->m_W - other.m_W);
        }
    }

    //------------------------------------------------------------------------------
    // operator* (scalar)
    //------------------------------------------------------------------------------
    //
    // Scales the quaternion by a scalar.
    //
    // Params:
    //  scalar - Scale factor.
    //
    // Returns:
    //  (x*scalar, y*scalar, z*scalar, w*scalar).
    //
    // Notes:
    //  Assumes finite scalar; asserts in debug.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::operator*(float scalar) const noexcept
    {
        assert(xmath::isFinite(scalar));
        if constexpr (V)
        {
            return fquat_t<V>(_mm_mul_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else
        {
            return fquat_t<V>(this->m_X * scalar, this->m_Y * scalar, this->m_Z * scalar, this->m_W * scalar);
        }
    }

    //------------------------------------------------------------------------------
    // operator/ (scalar)
    //------------------------------------------------------------------------------
    //
    // Divides the quaternion by a scalar.
    //
    // Params:
    //  scalar - Division factor.
    //
    // Returns:
    //  (x/scalar, y/scalar, z/scalar, w/scalar).
    //
    // Notes:
    //  Assumes non-zero scalar; asserts in debug.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::operator/(float scalar) const noexcept
    {
        assert(xmath::isFinite(scalar) && scalar != 0.0f);
        if constexpr (V)
        {
            return fquat_t<V>(_mm_div_ps(this->m_XYZW, _mm_set1_ps(scalar)));
        }
        else
        {
            float inv = 1.0f / scalar;
            return fquat_t<V>(this->m_X * inv, this->m_Y * inv, this->m_Z * inv, this->m_W * inv);
        }
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Multiplies two quaternions.
    //
    // Params:
    //  other - Second quaternion.
    //
    // Returns:
    //  Combined rotation (this * other).
    //
    // Notes:
    //  Uses SIMD optimized multiplication for V=true; right-hand rule.
    //  Asserts finite inputs in debug.
    //  Result may need normalization if inputs are not unit quaternions.
    //
    template <bool V>
    inline fquat_t<V> fquat_t<V>::operator*(const fquat_t& other) const noexcept
    {
        assert(isFinite() && other.isFinite());
        if constexpr (V)
        {
            const floatx4 w = _mm_mul_ps(this->m_XYZW, other.m_XYZW);
            fquat_t<V> result(_mm_sub_ps(
                _mm_add_ps(
                    _mm_add_ps(
                        _mm_mul_ps(other.m_XYZW, _mm_set1_ps(this->m_W)),
                        _mm_mul_ps(this->m_XYZW, _mm_set1_ps(other.m_W))
                    ),
                    _mm_mul_ps(
                        _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(3, 0, 2, 1)),
                        _mm_shuffle_ps(other.m_XYZW, other.m_XYZW, _MM_SHUFFLE(3, 1, 0, 2))
                    )
                ),
                _mm_mul_ps(
                    _mm_shuffle_ps(this->m_XYZW, this->m_XYZW, _MM_SHUFFLE(3, 1, 0, 2)),
                    _mm_shuffle_ps(other.m_XYZW, other.m_XYZW, _MM_SHUFFLE(3, 0, 2, 1))
                )
            ));
            result.m_W = ((float*)&w)[3] - ((float*)&w)[2] - ((float*)&w)[1] - ((float*)&w)[0];
            assert(result.isFinite());
            return result;
        }
        else
        {
            fquat_t<V> result(
                this->m_W * other.m_X + this->m_X * other.m_W + this->m_Y * other.m_Z - this->m_Z * other.m_Y,
                this->m_W * other.m_Y - this->m_X * other.m_Z + this->m_Y * other.m_W + this->m_Z * other.m_X,
                this->m_W * other.m_Z + this->m_X * other.m_Y - this->m_Y * other.m_X + this->m_Z * other.m_W,
                this->m_W * other.m_W - this->m_X * other.m_X - this->m_Y * other.m_Y - this->m_Z * other.m_Z
            );
            assert(result.isFinite());
            return result;
        }
    }
    //------------------------------------------------------------------------------
    // operator+=
    //------------------------------------------------------------------------------
    //
    // Adds another quaternion in-place.
    //
    // Params:
    //  other - Quaternion to add.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses operator+; not typical for rotations.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::operator+=(const fquat_t& other) noexcept
    {
        *this = *this + other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator-=
    //------------------------------------------------------------------------------
    //
    // Subtracts another quaternion in-place.
    //
    // Params:
    //  other - Quaternion to subtract.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses operator-; not typical for rotations.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::operator-=(const fquat_t& other) noexcept
    {
        *this = *this - other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*=(scalar)
    //------------------------------------------------------------------------------
    //
    // Scales the quaternion in-place by a scalar.
    //
    // Params:
    //  scalar - Scale factor.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses operator*; asserts finite scalar.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::operator*=(float scalar) noexcept
    {
        *this = *this * scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator/=(scalar)
    //------------------------------------------------------------------------------
    //
    // Divides the quaternion in-place by a scalar.
    //
    // Params:
    //  scalar - Division factor.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses operator/; asserts non-zero scalar.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::operator/=(float scalar) noexcept
    {
        *this = *this / scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*=(quaternion)
    //------------------------------------------------------------------------------
    //
    // Multiplies the quaternion in-place by another quaternion.
    //
    // Params:
    //  other - Quaternion to multiply.
    //
    // Returns:
    //  Reference to this quaternion (chainable).
    //
    // Notes:
    //  Uses operator*; applies rotation composition.
    //
    template <bool V>
    inline fquat_t<V>& fquat_t<V>::operator*=(const fquat_t& other) noexcept
    {
        *this = *this * other;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator==
    //------------------------------------------------------------------------------
    //
    // Checks if two quaternions are exactly equal.
    //
    // Params:
    //  other - Quaternion to compare.
    //
    // Returns:
    //  True if all components are identical.
    //
    // Notes:
    //  Use Equals for tolerance-based comparison.
    //
    template <bool V>
    inline bool fquat_t<V>::operator==(const fquat_t& other) const noexcept
    {
        if constexpr (V)
        {
            return _mm_movemask_ps(_mm_cmpeq_ps(this->m_XYZW, other.m_XYZW)) == 0xF;
        }
        else
        {
            return this->m_X == other.m_X && this->m_Y == other.m_Y &&
                   this->m_Z == other.m_Z && this->m_W == other.m_W;
        }
    }

    //------------------------------------------------------------------------------
    // operator!=
    //------------------------------------------------------------------------------
    //
    // Checks if two quaternions are not equal.
    //
    // Params:
    //  other - Quaternion to compare.
    //
    // Returns:
    //  True if any components differ.
    //
    // Notes:
    //  Uses operator==.
    //
    template <bool V>
    inline bool fquat_t<V>::operator!=(const fquat_t& other) const noexcept
    {
        return !(*this == other);
    }

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses a quaternion component by index.
    //
    // Params:
    //  index - Component index (0=x, 1=y, 2=z, 3=w).
    //
    // Returns:
    //  Component value.
    //
    // Notes:
    //  Asserts valid index; assumes finite.
    //
    template <bool V>
    inline float fquat_t<V>::operator[](std::int32_t index) const noexcept
    {
        assert(index >= 0 && index <= 3);
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // operator[] (mutable)
    //------------------------------------------------------------------------------
    //
    // Accesses a quaternion component by index for modification.
    //
    // Params:
    //  index - Component index (0=x, 1=y, 2=z, 3=w).
    //
    // Returns:
    //  Reference to component.
    //
    // Notes:
    //  Asserts valid index.
    //
    template <bool V>
    inline float& fquat_t<V>::operator[](std::int32_t index) noexcept
    {
        assert(index >= 0 && index <= 3);
        return this->m_Elements[index];
    }

    //------------------------------------------------------------------------------
    // operator* (scalar, quaternion)
    //------------------------------------------------------------------------------
    //
    // Scales a quaternion by a scalar (friend).
    //
    // Params:
    //  scalar - Scale factor.
    //  q - Quaternion to scale.
    //
    // Returns:
    //  (x*scalar, y*scalar, z*scalar, w*scalar).
    //
    // Notes:
    //  Commutative; uses operator*.
    //
    template <bool V>
    inline fquat_t<V> operator*(float scalar, const fquat_t<V>& q) noexcept
    {
        return q * scalar;
    }

    //------------------------------------------------------------------------------
    // operator- (unary)
    //------------------------------------------------------------------------------
    //
    // Returns the negation of the quaternion.
    //
    // Params:
    //  q - Quaternion to negate.
    //
    // Returns:
    //  (-x, -y, -z, -w).
    //
    // Notes:
    //  Component-wise negation; not typically used for rotations.
    //
    template <bool V>
    inline fquat_t<V> operator-(const fquat_t<V>& q) noexcept
    {
        if constexpr (V)
        {
            return fquat_t<V>(_mm_sub_ps(_mm_setzero_ps(), q.m_XYZW));
        }
        else
        {
            return fquat_t<V>(-q.m_X, -q.m_Y, -q.m_Z, -q.m_W);
        }
    }

    //------------------------------------------------------------------------------
    // operator* (quaternion, vector)
    //------------------------------------------------------------------------------
    //
    // Rotates a vector by the quaternion.
    //
    // Params:
    //  q - Quaternion.
    //  v - Vector to rotate.
    //
    // Returns:
    //  Rotated vector q * v * q^-1.
    //
    // Notes:
    //  Delegates to RotateVector; assumes normalized quaternion.
    //
    template <bool V>
    inline fvec3 operator*(const fquat_t<V>& q, const fvec3& v) noexcept
    {
        return q.RotateVector(v);
    }
}
