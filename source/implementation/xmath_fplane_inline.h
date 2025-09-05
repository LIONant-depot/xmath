#pragma once
namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from individual components.
    //
    // Params:
    //  X - A component (normal X).
    //  Y - B component (normal Y).
    //  Z - C component (normal Z).
    //  D - D component (offset).
    //
    // Notes:
    //  Uses _mm_set_ps for SIMD; assumes normalized normal.
    //  Asserts finite components in debug.
    //
    template <bool V >
    constexpr fplane_t<V>::fplane_t(float X, float Y, float Z, float D) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZD = _mm_set_ps(-D, Z, Y, X);
        }
        else
        {
            this->m_X = X;
            this->m_Y = Y;
            this->m_Z = Z;
            this->m_D = -D;
        }
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from normal and distance.
    //
    // Params:
    //  Normal - Normalized plane normal.
    //  Distance - Distance from origin along normal.
    //
    // Notes:
    //  Sets D = -Normal.Dot(point on plane); here point = Normal * Distance.
    //  Assumes normalized Normal; asserts finite.
    //
    template <bool V >
    constexpr fplane_t<V>::fplane_t(const fvec3& Normal, float Distance) noexcept
    {
        *this = fromNormalDistance(Normal, Distance);
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from normal and point on plane.
    //
    // Params:
    //  Normal - Normalized plane normal.
    //  Point - Point on the plane.
    //
    // Notes:
    //  Sets D = -Normal.Dot(Point).
    //  Assumes normalized Normal; asserts finite.
    //
    template <bool V >
    inline fplane_t<V>::fplane_t(const fvec3& Normal, const fvec3& Point) noexcept
    {
        *this = fromNormalPoint(Normal, Point);
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from three points on plane.
    //
    // Params:
    //  P1 - First point.
    //  P2 - Second point.
    //  P3 - Third point.
    //
    // Notes:
    //  Computes Normal = (P2-P1).Cross(P3-P1).Normalize; D = -Normal.Dot(P1).
    //  Assumes non-collinear points; asserts finite and valid.
    //
    template <bool V >
    inline fplane_t<V>::fplane_t(const fvec3& P1, const fvec3& P2, const fvec3& P3) noexcept
    {
        *this = fromThreePoints(P1, P2, P3);
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from SIMD register.
    //
    // Params:
    //  Reg - __m128 register with A,B,C,D.
    //
    // Notes:
    //  Direct assignment.
    //
    template <bool V >
    constexpr fplane_t<V>::fplane_t(const floatx4& Reg) noexcept requires V
        : details::f_plane::simd_data{ .m_XYZD = Reg }
    {}

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from other SIMD variant.
    //
    // Params:
    //  Other - Other fplane_t.
    //
    // Notes:
    //  Copies components.
    //
    template <bool V >
    constexpr fplane_t<V>::fplane_t(const fplane_t<!V>& Other) noexcept
    {
        if constexpr (V)
        {
            this->m_XYZD = _mm_set_ps(Other.m_D, Other.m_Z, Other.m_Y, Other.m_X);
        }
        else
        {
            this->m_X = Other.m_X;
            this->m_Y = Other.m_Y;
            this->m_Z = Other.m_Z;
            this->m_D = Other.m_D;
        }
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from double array.
    //
    // Params:
    //  Conversion - Array of 4 doubles.
    //
    // Notes:
    //  Static cast to float.
    //
    template <bool V >
    constexpr fplane_t<V>::fplane_t(const std::array<double, 4>& Conversion) noexcept
    {
        this->m_X = static_cast<float>(Conversion[0]);
        this->m_Y = static_cast<float>(Conversion[1]);
        this->m_Z = static_cast<float>(Conversion[2]);
        this->m_D = static_cast<float>(Conversion[3]);
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // Constructor from float array.
    //
    // Params:
    //  Conversion - Array of 4 floats.
    //
    template <bool V >
    constexpr fplane_t<V>::fplane_t(const std::array<float, 4>& Conversion) noexcept
    {
        this->m_X = Conversion[0];
        this->m_Y = Conversion[1];
        this->m_Z = Conversion[2];
        this->m_D = Conversion[3];
    }

    //------------------------------------------------------------------------------
    // Assignment and conversion operators
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator std::array<double,4>
    //------------------------------------------------------------------------------
    //
    // Converts to double array.
    //
    // Returns:
    //  {X, Y, Z, D} as doubles.
    //
    template <bool V >
    constexpr fplane_t<V>::operator std::array<double, 4>(void) const noexcept
    {
        return { static_cast<double>(this->m_X), static_cast<double>(this->m_Y), static_cast<double>(this->m_Z), static_cast<double>(this->m_D) };
    }

    //------------------------------------------------------------------------------
    // operator std::array<float,4>
    //------------------------------------------------------------------------------
    //
    // Converts to float array.
    //
    // Returns:
    //  {X, Y, Z, D} as floats.
    //
    template <bool V >
    constexpr fplane_t<V>::operator std::array<float, 4>(void) const noexcept
    {
        return { this->m_X, this->m_Y, this->m_Z, this->m_D };
    }

    //------------------------------------------------------------------------------
    // operator std::string
    //------------------------------------------------------------------------------
    //
    // Converts to string representation.
    //
    // Returns:
    //  "fplane_t<X, Y, Z, D>"
    //
    template <bool V >
    inline fplane_t<V>::operator std::string(void) const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // ToString
    //------------------------------------------------------------------------------
    //
    // Returns string representation.
    //
    // Returns:
    //  "fplane_t<X, Y, Z, D>"
    //
    template <bool V >
    inline std::string fplane_t<V>::ToString(void) const noexcept
    {
        return std::format("plane({}, {}, {}, {})", this->m_X, this->m_Y, this->m_Z, this->m_D);
    }

    //------------------------------------------------------------------------------
    // operator<<
    //------------------------------------------------------------------------------
    //
    // Streams the plane to ostream.
    //
    // Params:
    //  os - Output stream.
    //  plane - Plane to output.
    //
    // Returns:
    //  Reference to os.
    //
    template <bool V >
    inline std::ostream& operator<< (std::ostream& os, const fplane_t<V>& plane) noexcept
    {
        os << plane.ToString();
        return os;
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Returns zero plane (0,0,0,0).
    //
    template <bool V >
    constexpr fplane_t<V> fplane_t<V>::fromZero(void) noexcept
    {
        if constexpr (V) return fplane_t{ floatx4{.m128_f32{0,0,0,0}} };
        else             return { 0.0f, 0.0f, 0.0f, 0.0f };
    }

    //------------------------------------------------------------------------------
    // fromIdentity
    //------------------------------------------------------------------------------
    //
    // Returns identity plane (0,0,1,0) - XY plane at origin.
    //
    template <bool V >
    constexpr fplane_t<V> fplane_t<V>::fromIdentity(void) noexcept
    {
        if constexpr (V) return fplane_t{ floatx4{.m128_f32{0,0,1,0}} };
        else             return { 0.0f, 0.0f, 1.0f, 0.0f };
    }

    //------------------------------------------------------------------------------
    // Static methods
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // fromNormalDistance
    //------------------------------------------------------------------------------
    //
    // Creates plane from normal and distance.
    //
    // Params:
    //  Normal - Normalized normal.
    //  Distance - Distance from origin.
    //
    // Returns:
    //  Plane with D = -Distance (assuming normalized normal).
    //
    // Notes:
    //  Assumes normalized Normal; D = -Normal.Dot(Normal * Distance).
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::fromNormalDistance(const fvec3& Normal, float Distance) noexcept
    {
        assert(Normal.isFinite() && Normal.isNormalized(1e-6f) && xmath::isFinite(Distance));
        fplane_t<V> result;
        result.m_X = Normal.m_X;
        result.m_Y = Normal.m_Y;
        result.m_Z = Normal.m_Z;
        result.m_D = -Distance;
        return result;
    }

    //------------------------------------------------------------------------------
    // fromNormalPoint
    //------------------------------------------------------------------------------
    //
    // Creates plane from normal and point.
    //
    // Params:
    //  Normal - Normalized normal.
    //  Point - Point on plane.
    //
    // Returns:
    //  Plane with D = -Normal.Dot(Point).
    //
    // Notes:
    //  Assumes normalized Normal.
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::fromNormalPoint(const fvec3& Normal, const fvec3& Point) noexcept
    {
        assert(Normal.isFinite() && Point.isFinite() && Normal.isNormalized(1e-6f));
        fplane_t<V> result;
        result.m_X = Normal.m_X;
        result.m_Y = Normal.m_Y;
        result.m_Z = Normal.m_Z;
        result.m_D = -Normal.Dot(Point);
        return result;
    }

    //------------------------------------------------------------------------------
    // fromThreePoints
    //------------------------------------------------------------------------------
    //
    // Creates plane from three points.
    //
    // Params:
    //  P1 - First point.
    //  P2 - Second point.
    //  P3 - Third point.
    //
    // Returns:
    //  Plane with Normal = (P2-P1).Cross(P3-P1).Normalize, D = -Normal.Dot(P1).
    //
    // Notes:
    //  Assumes non-collinear points; normalizes.
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::fromThreePoints(const fvec3& P1, const fvec3& P2, const fvec3& P3) noexcept
    {
        assert(P1.isFinite() && P2.isFinite() && P3.isFinite());
        fvec3 normal = (P3 - P1).Cross(P2 - P1).NormalizeSafeCopy();
        assert(!normal.isNearlyZero(1e-6f));
        return fromNormalPoint(normal, P1);
    }

    //------------------------------------------------------------------------------
    // IntersectThreePlanes
    //------------------------------------------------------------------------------
    //
    // Finds intersection point of three planes.
    //
    // Params:
    //  P1 - First plane.
    //  P2 - Second plane.
    //  P3 - Third plane.
    //  Result - Output intersection point.
    //
    // Returns:
    //  True if intersection exists (non-parallel).
    //
    // Notes:
    //  Solves linear system; uses determinant for validity.
    //
    template <bool V >
    inline bool fplane_t<V>::IntersectThreePlanes(const fplane_t& P1, const fplane_t& P2, const fplane_t& P3, fvec3& Result) noexcept
    {
        assert(P1.isFinite() && P2.isFinite() && P3.isFinite());
        fvec3 n1(P1.m_X, P1.m_Y, P1.m_Z);
        fvec3 n2(P2.m_X, P2.m_Y, P2.m_Z);
        fvec3 n3(P3.m_X, P3.m_Y, P3.m_Z);
        fvec3 cross12 = n1.Cross(n2);
        float det = cross12.Dot(n3);
        if (xmath::Abs(det) < 1e-6f) return false; // Parallel or degenerate
        fvec3 cross23 = n2.Cross(n3);
        fvec3 cross31 = n3.Cross(n1);
        Result = (cross23 * (-P1.m_D) + cross31 * (-P2.m_D) + cross12 * (-P3.m_D)) / det;
        return true;
    }

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes plane-point dot product.
    //
    // Params:
    //  Plane - Plane.
    //  Point - Point.
    //
    // Returns:
    //  X*Px + Y*Py + Z*Pz + D.
    //
    // Notes:
    //  Uses _mm_dp_ps for SIMD.
    //
    template <bool V >
    inline float fplane_t<V>::Dot(const fplane_t& Plane, const fvec3& Point) noexcept
    {
        if constexpr (V)
        {
            floatx4 point_vec = _mm_set_ps(0.0f, Point.m_Z, Point.m_Y, Point.m_X);
            return _mm_cvtss_f32(_mm_dp_ps(Plane.m_XYZD, point_vec, 0x71)) + Plane.m_D;
        }
        else
        {
            return Plane.m_X * Point.m_X + Plane.m_Y * Point.m_Y + Plane.m_Z * Point.m_Z + Plane.m_D;
        }
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Signed distance from point to plane.
    //
    // Params:
    //  Plane - Plane.
    //  Point - Point.
    //
    // Returns:
    //  Dot / normal length; signed.
    //
    // Notes:
    //  Assumes normalized plane normal.
    //
    template <bool V >
    inline float fplane_t<V>::Distance(const fplane_t& Plane, const fvec3& Point) noexcept
    {
        assert(Plane.isNormalized(1e-6f));
        return Dot(Plane, Point) / Plane.NormalLength();
    }

    //------------------------------------------------------------------------------
    // Side
    //------------------------------------------------------------------------------
    //
    // Determines which side of the plane the point is on.
    //
    // Params:
    //  Plane - Plane.
    //  Point - Point.
    //
    // Returns:
    //  1 (front), -1 (back), 0 (on plane).
    //
    // Notes:
    //  Uses sign of Dot; tolerance for zero.
    //
    template <bool V >
    inline std::int32_t fplane_t<V>::Side(const fplane_t& Plane, const fvec3& Point) noexcept
    {
        float dot = Dot(Plane, Point);
        if (xmath::Abs(dot) < 1e-6f) return 0;
        return dot > 0.0f ? 1 : -1;
    }

    //------------------------------------------------------------------------------
    // IntersectLine
    //------------------------------------------------------------------------------
    //
    // Intersects plane with infinite line.
    //
    // Params:
    //  Plane - Plane.
    //  T - Output parameter.
    //  Start - Line start.
    //  Direction - Line direction (normalized).
    //
    // Returns:
    //  True if intersects; T is parameter.
    //
    // Notes:
    //  T = -Dot(Plane, Start) / Dot(Plane.Normal, Direction).
    //
    template <bool V >
    inline bool fplane_t<V>::IntersectLine(const fplane_t& Plane, float& T, const fvec3& Start, const fvec3& Direction) noexcept
    {
        assert(Direction.isNormalized(1e-6f));
        float denom = Plane.Normal().Dot(Direction);
        if (xmath::Abs(denom) < 1e-6f) return false; // Parallel
        T = -Dot(Plane, Start) / denom;
        return true;
    }

    //------------------------------------------------------------------------------
    // IntersectLineSegment
    //------------------------------------------------------------------------------
    //
    // Intersects plane with line segment.
    //
    // Params:
    //  Plane - Plane.
    //  T - Output parameter (0-1).
    //  P0 - Segment start.
    //  P1 - Segment end.
    //
    // Returns:
    //  True if intersects within segment; T is fraction.
    //
    // Notes:
    //  Uses IntersectLine; checks T in [0,1].
    //
    template <bool V >
    inline bool fplane_t<V>::IntersectLineSegment(const fplane_t& Plane, float& T, const fvec3& P0, const fvec3& P1) noexcept
    {
        fvec3 dir = P1 - P0;
        float len = dir.Length();
        if (len < 1e-6f) return false; // Degenerate segment
        if (!IntersectLine(Plane, T, P0, dir / len)) return false;
        T /= len;
        return T >= 0.0f && T <= 1.0f;
    }

    //------------------------------------------------------------------------------
    // Static methods as members
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Instance version of static Dot.
    //
    template <bool V >
    inline float fplane_t<V>::Dot(const fvec3& Point) const noexcept
    {
        return Dot(*this, Point);
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Instance version of static Distance.
    //
    template <bool V >
    inline float fplane_t<V>::Distance(const fvec3& Point) const noexcept
    {
        return Distance(*this, Point);
    }

    //------------------------------------------------------------------------------
    // Side
    //------------------------------------------------------------------------------
    //
    // Instance version of static Side.
    //
    template <bool V >
    inline std::int32_t fplane_t<V>::Side(const fvec3& Point) const noexcept
    {
        return Side(*this, Point);
    }

    //------------------------------------------------------------------------------
    // IntersectLine
    //------------------------------------------------------------------------------
    //
    // Instance version of static IntersectLine.
    //
    template <bool V >
    inline bool fplane_t<V>::IntersectLine(float& T, const fvec3& Start, const fvec3& Direction) const noexcept
    {
        return IntersectLine(*this, T, Start, Direction);
    }

    //------------------------------------------------------------------------------
    // IntersectLineSegment
    //------------------------------------------------------------------------------
    //
    // Instance version of static IntersectLineSegment.
    //
    template <bool V >
    inline bool fplane_t<V>::IntersectLineSegment(float& T, const fvec3& P0, const fvec3& P1) const noexcept
    {
        return IntersectLineSegment(*this, T, P0, P1);
    }

    //------------------------------------------------------------------------------
    // Instance methods - Basic operations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // NormalLength
    //------------------------------------------------------------------------------
    //
    // Length of the normal vector.
    //
    // Returns:
    //  Sqrt(X^2 + Y^2 + Z^2).
    //
    template <bool V >
    inline float fplane_t<V>::NormalLength(void) const noexcept
    {
        return xmath::Sqrt(NormalLengthSq());
    }

    //------------------------------------------------------------------------------
    // NormalLengthSq
    //------------------------------------------------------------------------------
    //
    // Squared length of the normal vector.
    //
    // Returns:
    //  X^2 + Y^2 + Z^2.
    //
    template <bool V >
    inline float fplane_t<V>::NormalLengthSq(void) const noexcept
    {
        if constexpr (V)
        {
            return _mm_cvtss_f32(_mm_dp_ps(this->m_XYZD, this->m_XYZD, 0x71));
        }
        else
        {
            return this->m_X * this->m_X + this->m_Y * this->m_Y + this->m_Z * this->m_Z;
        }
    }

    //------------------------------------------------------------------------------
    // NormalizeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy of the plane.
    //
    // Returns:
    //  Plane with normal normalized, D scaled accordingly.
    //
    // Notes:
    //  Divides components by normal length; assumes non-zero.
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::NormalizeCopy(void) const noexcept
    {
        assert(!isNearlyZero(1e-6f));
        float len = NormalLength();
        assert(len > 0.0f);
        float inv_len = 1.0f / len;
        if constexpr (V)
        {
            return fplane_t<V>(_mm_mul_ps(this->m_XYZD, _mm_set1_ps(inv_len)));
        }
        else
        {
            fplane_t<V> result;
            result.m_X = this->m_X * inv_len;
            result.m_Y = this->m_Y * inv_len;
            result.m_Z = this->m_Z * inv_len;
            result.m_D = this->m_D * inv_len;
            return result;
        }
    }

    //------------------------------------------------------------------------------
    // Normalize
    //------------------------------------------------------------------------------
    //
    // Normalizes the plane in-place.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    // Notes:
    //  Uses NormalizeCopy.
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::Normalize(void) noexcept
    {
        *this = NormalizeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // NormalizeSafeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a safe normalized copy; falls back to identity if zero.
    //
    // Returns:
    //  Normalized or identity plane.
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::NormalizeSafeCopy(void) const noexcept
    {
        float len_sq = NormalLengthSq();
        if (len_sq < 1e-6f) return fromIdentity();
        float inv_len = 1.0f / xmath::Sqrt(len_sq);
        if constexpr (V)
        {
            return fplane_t<V>(_mm_mul_ps(this->m_XYZD, _mm_set1_ps(inv_len)));
        }
        else
        {
            return fplane_t<V>(this->m_X * inv_len, this->m_Y * inv_len, this->m_Z * inv_len, this->m_D * inv_len);
        }
    }

    //------------------------------------------------------------------------------
    // NormalizeSafe
    //------------------------------------------------------------------------------
    //
    // Safe normalizes in-place; falls back to identity if zero.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::NormalizeSafe(void) noexcept
    {
        *this = NormalizeSafeCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all components are finite.
    //
    // Returns:
    //  True if no NaN or Inf.
    //
    template <bool V >
    inline bool fplane_t<V>::isFinite(void) const noexcept
    {
        return xmath::isFinite(this->m_X) && xmath::isFinite(this->m_Y) && xmath::isFinite(this->m_Z) && xmath::isFinite(this->m_D);
    }

    //------------------------------------------------------------------------------
    // isNormalized
    //------------------------------------------------------------------------------
    //
    // Checks if normal is unit length within tolerance.
    //
    // Params:
    //  Tolerance - Epsilon for comparison.
    //
    // Returns:
    //  True if |length^2 - 1| < tolerance.
    //
    template <bool V >
    inline bool fplane_t<V>::isNormalized(float Tolerance) const noexcept
    {
        return xmath::Abs(NormalLengthSq() - 1.0f) < Tolerance;
    }

    //------------------------------------------------------------------------------
    // isNearlyZero
    //------------------------------------------------------------------------------
    //
    // Checks if plane is nearly zero within tolerance.
    //
    // Params:
    //  Tolerance - Epsilon for comparison.
    //
    // Returns:
    //  True if all components < tolerance.
    //
    template <bool V >
    inline bool fplane_t<V>::isNearlyZero(float Tolerance) const noexcept
    {
        return xmath::Abs(this->m_X) < Tolerance && xmath::Abs(this->m_Y) < Tolerance &&
            xmath::Abs(this->m_Z) < Tolerance && xmath::Abs(this->m_D) < Tolerance;
    }

    //------------------------------------------------------------------------------
    // Equals
    //------------------------------------------------------------------------------
    //
    // Checks if two planes are equal within tolerance.
    //
    // Params:
    //  Other - Other plane.
    //  Tolerance - Epsilon for comparison.
    //
    // Returns:
    //  True if all components differ by < tolerance.
    //
    template <bool V >
    inline bool fplane_t<V>::Equals(const fplane_t& Other, float Tolerance) const noexcept
    {
        return xmath::Abs(this->m_X - Other.m_X) < Tolerance
            && xmath::Abs(this->m_Y - Other.m_Y) < Tolerance
            && xmath::Abs(this->m_Z - Other.m_Z) < Tolerance
            && xmath::Abs(this->m_D - Other.m_D) < Tolerance;
    }

    //------------------------------------------------------------------------------
    // Instance methods - Plane specifics
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Normal
    //------------------------------------------------------------------------------
    //
    // Returns the plane normal.
    //
    // Returns:
    //  fvec3(X, Y, Z).
    //
    template <bool V >
    inline fvec3 fplane_t<V>::Normal(void) const noexcept
    {
        return fvec3(this->m_X, this->m_Y, this->m_Z);
    }

    //------------------------------------------------------------------------------
    // Offset
    //------------------------------------------------------------------------------
    //
    // Returns the plane offset D.
    //
    // Returns:
    //  D.
    //
    template <bool V >
    inline float fplane_t<V>::Offset(void) const noexcept
    {
        return this->m_D;
    }

    //------------------------------------------------------------------------------
    // GetOrigin
    //------------------------------------------------------------------------------
    //
    // Returns a point on the plane (projection of origin).
    //
    // Returns:
    //  -Normal * (D / Normal.LengthSq()) if not zero, else (0,0,0).
    //
    template <bool V >
    inline fvec3 fplane_t<V>::GetOrigin(void) const noexcept
    {
        float len_sq = NormalLengthSq();
        if (len_sq < 1e-6f) return fvec3::fromZero();
        return -Normal() * (this->m_D / len_sq);
    }

    //------------------------------------------------------------------------------
    // ComputeOffset
    //------------------------------------------------------------------------------
    //
    // Computes D from a point on the plane.
    //
    // Params:
    //  Point - Point on plane.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    // Notes:
    //  D = -Normal.Dot(Point).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::ComputeOffset(const fvec3& Point) noexcept
    {
        this->m_D = -Normal().Dot(Point);
        return *this;
    }

    //------------------------------------------------------------------------------
    // DecomposeVector
    //------------------------------------------------------------------------------
    //
    // Decomposes vector into parallel and perpendicular to normal.
    //
    // Params:
    //  V - Vector to decompose.
    //  Parallel - Output parallel component.
    //  Perpendicular - Output perpendicular component.
    //
    // Notes:
    //  Parallel = (V.Dot(Normal) / Normal.LengthSq()) * Normal.
    //  Perpendicular = V - Parallel.
    //
    template <bool V >
    inline void fplane_t<V>::DecomposeVector(const fvec3& V, fvec3& Parallel, fvec3& Perpendicular) const noexcept
    {
        fvec3 n = Normal();
        float dot = V.Dot(n);
        float len_sq = NormalLengthSq();
        Parallel = (dot / len_sq) * n;
        Perpendicular = V - Parallel;
    }

    //------------------------------------------------------------------------------
    // OrthoVectors
    //------------------------------------------------------------------------------
    //
    // Computes two orthogonal vectors in the plane.
    //
    // Params:
    //  AxisA - Output first ortho vector.
    //  AxisB - Output second ortho vector.
    //
    // Notes:
    //  Finds arbitrary perpendicular to normal, then cross for second.
    //
    template <bool V >
    inline void fplane_t<V>::OrthoVectors(fvec3& AxisA, fvec3& AxisB) const noexcept
    {
        fvec3 n = Normal();
        AxisA = n.Perpendicular(fvec3(0.0f, 1.0f, 0.0f)).NormalizeSafeCopy();
        if (AxisA.isNearlyZero(1e-6f))
        {
            AxisA = n.Perpendicular(fvec3(0.0f, 0.0f, 1.0f)).NormalizeSafeCopy();
        }
        AxisB = n.Cross(AxisA).NormalizeSafeCopy();
    }

    //------------------------------------------------------------------------------
    // ReflectVector
    //------------------------------------------------------------------------------
    //
    // Reflects a vector over the plane.
    //
    // Params:
    //  V - Vector to reflect.
    //
    // Returns:
    //  V - 2 * proj_normal(V).
    //
    // Notes:
    //  Uses projection; assumes normalized normal.
    //
    template <bool V >
    inline fvec3 fplane_t<V>::ReflectVector(const fvec3& V) const noexcept
    {
        fvec3 n = Normal();
        return V - 2.0f * (V.Dot(n)) * n;
    }

    //------------------------------------------------------------------------------
    // Project
    //------------------------------------------------------------------------------
    //
    // Projects a point onto the plane.
    //
    // Params:
    //  Point - Point to project.
    //
    // Returns:
    //  Point - Distance * Normal.
    //
    // Notes:
    //  Uses signed distance.
    //
    template <bool V >
    inline fvec3 fplane_t<V>::Project(const fvec3& Point) const noexcept
    {
        return Point - Distance(Point) * Normal();
    }

    //------------------------------------------------------------------------------
    // IsPointOnPlane
    //------------------------------------------------------------------------------
    //
    // Checks if point is on the plane within epsilon.
    //
    // Params:
    //  Point - Point to check.
    //  Epsilon - Tolerance.
    //
    // Returns:
    //  True if Abs(Dot) < epsilon.
    //
    template <bool V >
    inline bool fplane_t<V>::IsPointOnPlane(const fvec3& Point, float Epsilon) const noexcept
    {
        return xmath::Abs(Dot(Point)) < Epsilon;
    }

    //------------------------------------------------------------------------------
    // SameSide
    //------------------------------------------------------------------------------
    //
    // Checks if two points are on the same side of the plane.
    //
    // Params:
    //  P0 - First point.
    //  P1 - Second point.
    //
    // Returns:
    //  True if Side(P0) == Side(P1).
    //
    template <bool V >
    inline bool fplane_t<V>::SameSide(const fvec3& P0, const fvec3& P1) const noexcept
    {
        return Side(P0) == Side(P1);
    }

    //------------------------------------------------------------------------------
    // ClipNGon
    //------------------------------------------------------------------------------
    //
    // Clips an n-gon polygon against the plane.
    //
    // Params:
    //  Dst - Output vertices.
    //  DstCount - Output count.
    //  Src - Input vertices.
    //  SrcCount - Input count.
    //
    // Returns:
    //  True if any output; Dst filled with clipped polygon.
    //
    // Notes:
    //  Uses Sutherland-Hodgman; assumes convex input.
    //
    template <bool V >
    inline bool fplane_t<V>::ClipNGon(fvec3* Dst, std::int32_t& DstCount, const fvec3* Src, std::int32_t SrcCount) const noexcept
    {
        DstCount = 0;
        if (SrcCount < 3) return false;
        fvec3 prev = Src[SrcCount - 1];
        float prev_dot = Dot(prev);
        bool prev_inside = prev_dot >= 0.0f;

        for (int i = 0; i < SrcCount; ++i)
        {
            fvec3 curr = Src[i];
            float curr_dot = Dot(curr);
            bool curr_inside = curr_dot >= 0.0f;

            if (curr_inside != prev_inside)
            {
                float t = prev_dot / (prev_dot - curr_dot);
                Dst[DstCount++] = prev + t * (curr - prev);
            }

            if (curr_inside)
            {
                Dst[DstCount++] = curr;
            }

            prev = curr;
            prev_dot = curr_dot;
            prev_inside = curr_inside;
        }

        return DstCount >= 3;
    }

    //------------------------------------------------------------------------------
    // Instance methods - Setup operations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // setupFromNormalDistance
    //------------------------------------------------------------------------------
    //
    // Sets plane from normal and distance.
    //
    // Params:
    //  Normal - Normalized normal.
    //  Distance - Distance.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::setupFromNormalDistance(const fvec3& Normal, float Distance) noexcept
    {
        *this = fromNormalDistance(Normal, Distance);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupFromNormalPoint
    //------------------------------------------------------------------------------
    //
    // Sets plane from normal and point.
    //
    // Params:
    //  Normal - Normalized normal.
    //  Point - Point.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::setupFromNormalPoint(const fvec3& Normal, const fvec3& Point) noexcept
    {
        *this = fromNormalPoint(Normal, Point);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupFromThreePoints
    //------------------------------------------------------------------------------
    //
    // Sets plane from three points.
    //
    // Params:
    //  P1 - First point.
    //  P2 - Second point.
    //  P3 - Third point.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::setupFromThreePoints(const fvec3& P1, const fvec3& P2, const fvec3& P3) noexcept
    {
        *this = fromThreePoints(P1, P2, P3);
        return *this;
    }

    //------------------------------------------------------------------------------
    // TranslateCopy
    //------------------------------------------------------------------------------
    //
    // Returns a translated copy of the plane.
    //
    // Params:
    //  Translation - Translation vector.
    //
    // Returns:
    //  Plane with D adjusted by -Normal.Dot(Translation).
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::TranslateCopy(const fvec3& Translation) const noexcept
    {
        fplane_t<V> result = *this;
        result.m_D -= Normal().Dot(Translation);
        return result;
    }

    //------------------------------------------------------------------------------
    // Translate
    //------------------------------------------------------------------------------
    //
    // Translates the plane in-place.
    //
    // Params:
    //  Translation - Translation vector.
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::Translate(const fvec3& Translation) noexcept
    {
        *this = TranslateCopy(Translation);
        return *this;
    }

    //------------------------------------------------------------------------------
    // FlipCopy
    //------------------------------------------------------------------------------
    //
    // Returns a flipped copy of the plane (reversed normal).
    //
    // Returns:
    //  Plane with -Normal, -D.
    //
    template <bool V >
    inline fplane_t<V> fplane_t<V>::FlipCopy(void) const noexcept
    {
        return -(*this);
    }

    //------------------------------------------------------------------------------
    // Flip
    //------------------------------------------------------------------------------
    //
    // Flips the plane in-place (reverses normal).
    //
    // Returns:
    //  Reference to this plane (chainable).
    //
    template <bool V >
    inline fplane_t<V>& fplane_t<V>::Flip(void) noexcept
    {
        *this = FlipCopy();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Operator overloads
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Returns the negated plane.
    //
    // Returns:
    //  (-X, -Y, -Z, -D).
    //
    // Notes:
    //  Flips normal and offset.
    //
    template <bool V>
    inline fplane_t<V> fplane_t<V>::operator-() const noexcept
    {
        if constexpr (V)
        {
            return fplane_t<V>(_mm_sub_ps(_mm_setzero_ps(), this->m_XYZD));
        }
        else
        {
            fplane_t<V> result;
            result.m_X = -this->m_X;
            result.m_Y = -this->m_Y;
            result.m_Z = -this->m_Z;
            result.m_D = -this->m_D;
            return result;
        }
    }
    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses a plane component by index.
    //
    // Params:
    //  Index - Component index (0=X, 1=Y, 2=Z, 3=D).
    //
    // Returns:
    //  Component value.
    //
    // Notes:
    //  Asserts valid index.
    //
    template <bool V >
    inline float fplane_t<V>::operator[] (std::int32_t Index) const noexcept
    {
        assert(Index >= 0 && Index <= 3);
        return this->m_Elements[Index];
    }

    //------------------------------------------------------------------------------
    // operator[] (mutable)
    //------------------------------------------------------------------------------
    //
    // Accesses a plane component by index for modification.
    //
    // Params:
    //  Index - Component index (0=X, 1=Y, 2=Z, 3=D).
    //
    // Returns:
    //  Reference to component.
    //
    // Notes:
    //  Asserts valid index.
    //
    template <bool V >
    inline float& fplane_t<V>::operator[] (std::int32_t Index) noexcept
    {
        assert(Index >= 0 && Index <= 3);
        return this->m_Elements[Index];
    }

    //------------------------------------------------------------------------------
    // operator==
    //------------------------------------------------------------------------------
    //
    // Checks if two planes are exactly equal.
    //
    // Params:
    //  Other - Other plane.
    //
    // Returns:
    //  True if all components identical.
    //
    // Notes:
    //  Use Equals for tolerance.
    //
    template <bool V >
    inline bool fplane_t<V>::operator== (const fplane_t& Other) const noexcept
    {
        if constexpr (V)
        {
            return _mm_movemask_ps(_mm_cmpeq_ps(this->m_XYZD, Other.m_XYZD)) == 0xF;
        }
        else
        {
            return this->m_X == Other.m_X && this->m_Y == Other.m_Y &&
                this->m_Z == Other.m_Z && this->m_D == Other.m_D;
        }
    }

    //------------------------------------------------------------------------------
    // operator!=
    //------------------------------------------------------------------------------
    //
    // Checks if two planes are not equal.
    //
    // Params:
    //  Other - Other plane.
    //
    // Returns:
    //  True if any components differ.
    //
    // Notes:
    //  Uses operator==.
    //
    template <bool V >
    inline bool fplane_t<V>::operator!= (const fplane_t& Other) const noexcept
    {
        return !(*this == Other);
    }

    //------------------------------------------------------------------------------
    // Friend operators
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // operator* (matrix, plane)
    //------------------------------------------------------------------------------
    //
    // Transforms a plane by a matrix.
    //
    // Params:
    //  M - Transformation matrix.
    //  P - Plane to transform.
    //
    // Returns:
    //  Transformed plane.
    //
    // Notes:
    //  Transforms normal by inverse transpose; adjusts D.
    //  Assumes affine matrix; uses fmat4 methods.
    //
    template<bool V>
    fplane_t<V> operator* (const fmat4& M, const fplane_t<V>& P) noexcept
    {
        fvec3 normal(P.m_X, P.m_Y, P.m_Z);
        fvec3 point = P.GetOrigin();
        fvec3 new_point = M * point;
        fvec3 new_normal = M.InverseTranspose().TransformNormal(normal).NormalizeSafeCopy();
        return fplane_t<V>(new_normal, new_point);
    }
}