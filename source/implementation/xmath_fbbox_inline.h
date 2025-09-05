namespace xmath
{
    //------------------------------------------------------------------------------
    // Constructors
    //------------------------------------------------------------------------------
    //
    // Constructs a bounding box from a single point, setting both min and max to the point.
    //
    // Parameters:
    // Point - The point defining both min and max corners.
    //
    template <bool V>
    constexpr fbbox_t<V>::fbbox_t(const fvec3& Point) noexcept
    {
        this->m_Min = Point;
        this->m_Max = Point;
    }

    //------------------------------------------------------------------------------
    // Constructs a bounding box from minimum and maximum points.
    //
    // Parameters:
    // Min - The minimum corner of the box.
    // Max - The maximum corner of the box.
    //
    template <bool V>
    constexpr fbbox_t<V>::fbbox_t(const fvec3& Min, const fvec3& Max) noexcept
    {
        this->m_Min = Min;
        this->m_Max = Max;
    }

    //------------------------------------------------------------------------------
    // Constructs a bounding box from an array of vertices by computing their bounds.
    //
    // Parameters:
    // Verts - Array of vertices to enclose.
    // NVerts - Number of vertices in the array.
    //
    // Notes:
    // Asserts NVerts > 0 and Verts is not null.
    //
    template <bool V>
    inline fbbox_t<V>::fbbox_t(const fvec3* Verts, std::int32_t NVerts) noexcept
    {
        assert(NVerts > 0 && Verts);
        setupFromVerts(Verts, NVerts);
    }

    //------------------------------------------------------------------------------
    // Constructs a bounding box from a center point and radius, approximating a sphere.
    //
    // Parameters:
    // Center - The center point of the box.
    // Radius - The radius to extend from the center in all directions.
    //
    template <bool V>
    constexpr fbbox_t<V>::fbbox_t(const fvec3& Center, float Radius) noexcept
    {
        this->m_Min = Center - fvec3(Radius);
        this->m_Max = Center + fvec3(Radius);
    }

    //------------------------------------------------------------------------------
    // Constructs a bounding box from SIMD min/max registers (SIMD only).
    //
    // Parameters:
    // MinMax - Array of two __m128 registers for min and max corners.
    //
    // Notes:
    // Only available when T_USE_SIMD_V is true.
    //
    template <bool V>
    constexpr fbbox_t<V>::fbbox_t(const floatx4 MinMax[2]) noexcept requires V
    {
        this->m_MinMax[0] = MinMax[0];
        this->m_MinMax[1] = MinMax[1];
    }

    //------------------------------------------------------------------------------
    // Constructs a bounding box from a box with opposite SIMD setting.
    //
    // Parameters:
    // Other - The source bounding box to convert.
    //
    template <bool V>
    constexpr fbbox_t<V>::fbbox_t(const fbbox_t<!V>& Other) noexcept
    {
        this->m_Min = Other.m_Min;
        this->m_Max = Other.m_Max;
    }

    //------------------------------------------------------------------------------
    // Constructs a bounding box from an array of 6 floats (min x,y,z, max x,y,z).
    //
    // Parameters:
    // Conversion - Array containing min (x,y,z) followed by max (x,y,z).
    //
    template <bool V>
    constexpr fbbox_t<V>::fbbox_t(const std::array<float, 6>& Conversion) noexcept
    {
        this->m_Min = fvec3(Conversion[0], Conversion[1], Conversion[2]);
        this->m_Max = fvec3(Conversion[3], Conversion[4], Conversion[5]);
    }

    //------------------------------------------------------------------------------
    // Assignment and conversion operators
    //------------------------------------------------------------------------------
    //
    // Converts the bounding box to an array of 6 floats (min x,y,z, max x,y,z).
    //
    // Returns:
    // Array containing min (x,y,z) followed by max (x,y,z).
    //
    template <bool V>
    constexpr fbbox_t<V>::operator std::array<float, 6>(void) const noexcept
    {
        return {
            this->m_Min.m_X, this->m_Min.m_Y, this->m_Min.m_Z,
            this->m_Max.m_X, this->m_Max.m_Y, this->m_Max.m_Z
        };
    }

    //------------------------------------------------------------------------------
    // Converts the bounding box to a string representation.
    //
    // Returns:
    // A string in the format "Min: (x,y,z), Max: (x,y,z)".
    //
    template <bool V>
    inline fbbox_t<V>::operator std::string(void) const noexcept
    {
        return toString();
    }

    //------------------------------------------------------------------------------
    // Returns a string representation of the bounding box.
    //
    // Returns:
    // A string in the format "Min: (x,y,z), Max: (x,y,z)".
    //
    template <bool V>
    inline std::string fbbox_t<V>::toString(void) const noexcept
    {
        return std::format("Min: {}, Max: {}", this->m_Min.ToString(), this->m_Max.ToString());
    }

    //------------------------------------------------------------------------------
    // Overloads the stream output operator to print the bounding box.
    //
    // Parameters:
    // Os - The output stream.
    // BBox - The bounding box to print.
    //
    // Returns:
    // Reference to the output stream.
    //
    template <bool V>
    inline std::ostream& operator<<(std::ostream& Os, const fbbox_t<V>& BBox) noexcept
    {
        return Os << BBox.toString();
    }

    //------------------------------------------------------------------------------
    // Static properties
    //------------------------------------------------------------------------------
    //
    // Returns a zero bounding box with min and max at (0,0,0).
    //
    // Returns:
    // A zeroed bounding box.
    //
    template <bool V>
    constexpr fbbox_t<V> fbbox_t<V>::fromZero(void) noexcept
    {
        return fbbox_t<V>(fvec3::fromZero(), fvec3::fromZero());
    }

    //------------------------------------------------------------------------------
    // Returns an identity bounding box with min at (-1,-1,-1) and max at (1,1,1).
    //
    // Returns:
    // An identity bounding box.
    //
    template <bool V>
    constexpr fbbox_t<V> fbbox_t<V>::fromIdentity(void) noexcept
    {
        return fbbox_t<V>(fvec3(-1.f), fvec3(1.f));
    }

    //------------------------------------------------------------------------------
    // Creates a bounding box from a single point.
    //
    // Parameters:
    // Point - The point defining both min and max corners.
    //
    // Returns:
    // A bounding box enclosing the point.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::fromPoint(const fvec3& Point) noexcept
    {
        return fbbox_t<V>(Point);
    }

    //------------------------------------------------------------------------------
    // Creates a bounding box from minimum and maximum points.
    //
    // Parameters:
    // Min - The minimum corner of the box.
    // Max - The maximum corner of the box.
    //
    // Returns:
    // A bounding box with specified min and max.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::fromMinMax(const fvec3& Min, const fvec3& Max) noexcept
    {
        return fbbox_t<V>(Min, Max);
    }

    //------------------------------------------------------------------------------
    // Creates a bounding box from a center point and radius.
    //
    // Parameters:
    // Center - The center point of the box.
    // Radius - The radius to extend from the center in all directions.
    //
    // Returns:
    // A bounding box approximating a sphere.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::fromCenterRadius(const fvec3& Center, float Radius) noexcept
    {
        return fbbox_t<V>(Center, Radius);
    }

    //------------------------------------------------------------------------------
    // Creates a bounding box from an array of vertices.
    //
    // Parameters:
    // Verts - Array of vertices to enclose.
    // NVerts - Number of vertices in the array.
    //
    // Returns:
    // A bounding box enclosing all vertices.
    //
    // Notes:
    // Asserts NVerts > 0 and Verts is not null.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::fromVerts(const fvec3* Verts, std::int32_t NVerts) noexcept
    {
        return fbbox_t<V>(Verts, NVerts);
    }

    //------------------------------------------------------------------------------
    // Static methods
    //------------------------------------------------------------------------------
    //
    // Checks if two bounding boxes intersect.
    //
    // Parameters:
    // BBox1 - The first bounding box.
    // BBox2 - The second bounding box.
    //
    // Returns:
    // True if the boxes overlap, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fbbox_t& BBox1, const fbbox_t& BBox2) noexcept
    {
        return BBox1.m_Max.m_X >= BBox2.m_Min.m_X && BBox1.m_Min.m_X <= BBox2.m_Max.m_X &&
               BBox1.m_Max.m_Y >= BBox2.m_Min.m_Y && BBox1.m_Min.m_Y <= BBox2.m_Max.m_Y &&
               BBox1.m_Max.m_Z >= BBox2.m_Min.m_Z && BBox1.m_Min.m_Z <= BBox2.m_Max.m_Z;
    }

    //------------------------------------------------------------------------------
    // Checks if a bounding box contains a point.
    //
    // Parameters:
    // BBox - The bounding box.
    // Point - The point to check.
    //
    // Returns:
    // True if the point is inside or on the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fbbox_t& BBox, const fvec3& Point) noexcept
    {
        return BBox.ContainsPoint(Point);
    }

    //------------------------------------------------------------------------------
    // Checks if a bounding box intersects a plane.
    //
    // Parameters:
    // BBox - The bounding box.
    // Plane - The plane to check against.
    //
    // Returns:
    // True if the box intersects the plane, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fbbox_t& BBox, const fplane& Plane) noexcept
    {
        fvec3 center = BBox.getCenter();
        fvec3 extent = BBox.getSize() * 0.5f;
        float dist = Plane.Distance(center);
        float maxExtent = fvec3::Dot(extent, Plane.m_Normal.AbsCopy());
        return std::abs(dist) <= maxExtent;
    }

    //------------------------------------------------------------------------------
    // Checks if a bounding box intersects a line segment.
    //
    // Parameters:
    // BBox - The bounding box.
    // T - Output parameter for the intersection parameter (0 to 1).
    // P0 - The start point of the line segment.
    // P1 - The end point of the line segment.
    //
    // Returns:
    // True if the line segment intersects the box, false otherwise.
    //
    // Notes:
    // T is set to the entry point of intersection if true.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fbbox_t& BBox, float& T, const fvec3& P0, const fvec3& P1) noexcept
    {
        fvec3 dir = P1 - P0;
        float tmin = 0.0f, tmax = 1.0f;
        for (int i = 0; i < 3; ++i)
        {
            if (std::abs(dir[i]) < 1e-6f)
            {
                if (P0[i] < BBox.m_Min[i] || P0[i] > BBox.m_Max[i])
                    return false;
            }
            else
            {
                float ood = 1.0f / dir[i];
                float t1 = (BBox.m_Min[i] - P0[i]) * ood;
                float t2 = (BBox.m_Max[i] - P0[i]) * ood;
                if (t1 > t2) std::swap(t1, t2);
                tmin = std::max(tmin, t1);
                tmax = std::min(tmax, t2);
                if (tmin > tmax) return false;
            }
        }
        T = tmin;
        return true;
    }

    //------------------------------------------------------------------------------
    // Checks if a bounding box intersects a set of vertices (convex hull).
    //
    // Parameters:
    // BBox - The bounding box.
    // Verts - Array of vertices.
    // NVerts - Number of vertices in the array.
    //
    // Returns:
    // True if any vertex is inside the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fbbox_t& BBox, const std::span<fvec3> Verts) noexcept
    {
        for (auto& E : Verts)
        {
            if (BBox.containsPoint(E))
                return true;
        }
        return false;
    }

    //------------------------------------------------------------------------------
    // Checks if a bounding box intersects a triangle.
    //
    // Parameters:
    // BBox - The bounding box.
    // P0 - First vertex of the triangle.
    // P1 - Second vertex of the triangle.
    // P2 - Third vertex of the triangle.
    //
    // Returns:
    // True if the triangle intersects the box, false otherwise.
    //
    // Notes:
    // Uses AABB of triangle and plane intersection for efficiency.
    //
    template <bool V>
    inline bool fbbox_t<V>::IntersectTriangle(const fbbox_t& BBox, const fvec3& P0, const fvec3& P1, const fvec3& P2) noexcept
    {
        fbbox_t<V> triBox(P0);
        triBox += P1;
        triBox += P2;
        if (!intersect(BBox, triBox)) return false;

        fvec3 normal = fvec3::Cross(P1 - P0, P2 - P0);
        fplane plane(normal, fvec3::Dot(normal, P0));
        if (!intersect(BBox, plane)) return false;

        for (int i = 0; i < 3; ++i)
        {
            float boxMin = BBox.m_Min[i], boxMax = BBox.m_Max[i];
            float triMin = std::min({ P0[i], P1[i], P2[i] });
            float triMax = std::max({ P0[i], P1[i], P2[i] });
            if (boxMax < triMin || triMax < boxMin) return false;
        }
        return true;
    }

    //------------------------------------------------------------------------------
    // Computes the intersection of two bounding boxes.
    //
    // Parameters:
    // BBox1 - The first bounding box.
    // BBox2 - The second bounding box.
    //
    // Returns:
    // A bounding box representing the intersection volume.
    //
    // Notes:
    // Returns an invalid box if the boxes do not intersect.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::Intersection(const fbbox_t& BBox1, const fbbox_t& BBox2) noexcept
    {
        return fbbox_t<V>(
            fvec3::Max(BBox1.m_Min, BBox2.m_Min),
            fvec3::Min(BBox1.m_Max, BBox2.m_Max)
        );
    }

    //------------------------------------------------------------------------------
    // Static methods as members
    //------------------------------------------------------------------------------
    //
    // Checks if this box intersects another box.
    //
    // Parameters:
    // Other - The other bounding box.
    //
    // Returns:
    // True if the boxes overlap, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fbbox_t& Other) const noexcept
    {
        return Intersect(*this, Other);
    }

    //------------------------------------------------------------------------------
    // Checks if this box contains a point.
    //
    // Parameters:
    // Point - The point to check.
    //
    // Returns:
    // True if the point is inside or on the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fvec3& Point) const noexcept
    {
        return ContainsPoint(Point);
    }

    //------------------------------------------------------------------------------
    // Checks if this box intersects a plane.
    //
    // Parameters:
    // Plane - The plane to check against.
    //
    // Returns:
    // True if the box intersects the plane, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const fplane& Plane) const noexcept
    {
        return Intersect(*this, Plane);
    }

    //------------------------------------------------------------------------------
    // Checks if this box intersects a line segment.
    //
    // Parameters:
    // T - Output parameter for the intersection parameter (0 to 1).
    // P0 - The start point of the line segment.
    // P1 - The end point of the line segment.
    //
    // Returns:
    // True if the line segment intersects the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(float& T, const fvec3& P0, const fvec3& P1) const noexcept
    {
        return Intersect(*this, T, P0, P1);
    }

    //------------------------------------------------------------------------------
    // Checks if this box intersects a set of vertices.
    //
    // Parameters:
    // Verts - Array of vertices.
    // NVerts - Number of vertices in the array.
    //
    // Returns:
    // True if any vertex is inside the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Intersect(const std::span<fvec3> Verts) const noexcept
    {
        return Intersect(*this, Verts);
    }

    //------------------------------------------------------------------------------
    // Checks if this box intersects a triangle.
    //
    // Parameters:
    // P0 - First vertex of the triangle.
    // P1 - Second vertex of the triangle.
    // P2 - Third vertex of the triangle.
    //
    // Returns:
    // True if the triangle intersects the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::IntersectTriangle(const fvec3& P0, const fvec3& P1, const fvec3& P2) const noexcept
    {
        return IntersectTriangle(*this, P0, P1, P2);
    }

    //------------------------------------------------------------------------------
    // Computes the intersection with another box.
    //
    // Parameters:
    // Other - The other bounding box.
    //
    // Returns:
    // A bounding box representing the intersection volume.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::Intersection(const fbbox_t& Other) const noexcept
    {
        return Intersection(*this, Other);
    }

    //------------------------------------------------------------------------------
    // Instance methods - Basic operations
    //------------------------------------------------------------------------------
    //
    // Checks if all components of the box are finite (not NaN or infinite).
    //
    // Returns:
    // True if all components are finite, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::isFinite(void) const noexcept
    {
        return this->m_Min.isFinite() && this->m_Max.isFinite();
    }

    //------------------------------------------------------------------------------
    // Checks if the box is valid (min <= max and finite).
    //
    // Returns:
    // True if min <= max for all axes and components are finite, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::isValid(void) const noexcept
    {
        return this->m_Min.m_X <= this->m_Max.m_X &&
            this->m_Min.m_Y <= this->m_Max.m_Y &&
            this->m_Min.m_Z <= this->m_Max.m_Z &&
            isFinite();
    }

    //------------------------------------------------------------------------------
    // Computes the size (extents) of the box.
    //
    // Returns:
    // A vector representing the size (max - min).
    //
    template <bool V>
    inline fvec3 fbbox_t<V>::getSize(void) const noexcept
    {
        return this->m_Max - this->m_Min;
    }

    //------------------------------------------------------------------------------
    // Computes the center of the box.
    //
    // Returns:
    // A vector representing the center point (average of min and max).
    //
    template <bool V>
    inline fvec3 fbbox_t<V>::getCenter(void) const noexcept
    {
        return (this->m_Min + this->m_Max) * 0.5f;
    }

    //------------------------------------------------------------------------------
    // Computes the radius of the box (half the diagonal length).
    //
    // Returns:
    // The radius as a float.
    //
    template <bool V>
    inline float fbbox_t<V>::getRadius(void) const noexcept
    {
        return getSize().Length() * 0.5f;
    }

    //------------------------------------------------------------------------------
    // Computes the squared radius of the box.
    //
    // Returns:
    // The squared radius as a float.
    //
    // Notes:
    // Faster than getRadius() for comparisons.
    //
    template <bool V>
    inline float fbbox_t<V>::getRadiusSquared(void) const noexcept
    {
        return getSize().LengthSq() * 0.25f;
    }

    //------------------------------------------------------------------------------
    // Computes the surface area of the box.
    //
    // Returns:
    // The surface area as a float.
    //
    template <bool V>
    inline float fbbox_t<V>::getSurfaceArea(void) const noexcept
    {
        fvec3 size = getSize();
        return 2.0f * (size.m_X * size.m_Y + size.m_Y * size.m_Z + size.m_Z * size.m_X);
    }

    //------------------------------------------------------------------------------
    // Checks if all components of min and max are within a given range.
    //
    // Parameters:
    // Min - The minimum value of the range.
    // Max - The maximum value of the range.
    //
    // Returns:
    // True if all components are in [Min, Max], false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::isInRange(float Min, float Max) const noexcept
    {
        return this->m_Min.isInRange(Min, Max) && this->m_Max.isInRange(Min, Max);
    }

    //------------------------------------------------------------------------------
    // Checks if this box fully contains another box.
    //
    // Parameters:
    // Other - The other bounding box to check.
    //
    // Returns:
    // True if this box contains Other, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::Contains(const fbbox_t& Other) const noexcept
    {
        return this->m_Min.m_X <= Other.m_Min.m_X && this->m_Max.m_X >= Other.m_Max.m_X &&
            this->m_Min.m_Y <= Other.m_Min.m_Y && this->m_Max.m_Y >= Other.m_Max.m_Y &&
            this->m_Min.m_Z <= Other.m_Min.m_Z && this->m_Max.m_Z >= Other.m_Max.m_Z;
    }

    //------------------------------------------------------------------------------
    // Checks if this box contains a point.
    //
    // Parameters:
    // Point - The point to check.
    //
    // Returns:
    // True if the point is inside or on the box, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::ContainsPoint(const fvec3& Point) const noexcept
    {
        return Point.m_X >= this->m_Min.m_X && Point.m_X <= this->m_Max.m_X &&
            Point.m_Y >= this->m_Min.m_Y && Point.m_Y <= this->m_Max.m_Y &&
            Point.m_Z >= this->m_Min.m_Z && Point.m_Z <= this->m_Max.m_Z;
    }

    //------------------------------------------------------------------------------
    // Computes the closest point on or in the box to a given point.
    //
    // Parameters:
    // Point - The reference point.
    //
    // Returns:
    // The closest point on or in the box.
    //
    template <bool V>
    inline fvec3 fbbox_t<V>::getClosestPoint(const fvec3& Point) const noexcept
    {
        return fvec3::Clamp(Point, this->m_Min, this->m_Max);
    }

    //------------------------------------------------------------------------------
    // Retrieves the 8 vertices of the box.
    //
    // Parameters:
    // Dst - Array to store the vertices (must have space for at least 8).
    // NVerts - Number of vertices to write (must be >= 8).
    //
    // Notes:
    // Asserts NVerts >= 8.
    //
    template <bool V>
    inline void fbbox_t<V>::getVerts(fvec3* Dst, std::int32_t NVerts) const noexcept
    {
        assert(NVerts >= 8);
        fvec3 min = this->m_Min;
        fvec3 max = this->m_Max;
        Dst[0] = min;
        Dst[1] = fvec3(min.m_X, min.m_Y, max.m_Z);
        Dst[2] = fvec3(min.m_X, max.m_Y, min.m_Z);
        Dst[3] = fvec3(min.m_X, max.m_Y, max.m_Z);
        Dst[4] = fvec3(max.m_X, min.m_Y, min.m_Z);
        Dst[5] = fvec3(max.m_X, min.m_Y, max.m_Z);
        Dst[6] = fvec3(max.m_X, max.m_Y, min.m_Z);
        Dst[7] = max;
    }

    //------------------------------------------------------------------------------
    // Instance methods - Setup operations
    //------------------------------------------------------------------------------
    //
    // Sets up the box from a single point, making min and max equal to the point.
    //
    // Parameters:
    // Point - The point defining both min and max corners.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::setupFromPoint(const fvec3& Point) noexcept
    {
        this->m_Min = Point;
        this->m_Max = Point;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Sets up the box from minimum and maximum points.
    //
    // Parameters:
    // Min - The minimum corner of the box.
    // Max - The maximum corner of the box.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::setupFromMinMax(const fvec3& Min, const fvec3& Max) noexcept
    {
        this->m_Min = Min;
        this->m_Max = Max;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Sets up the box from a center point and radius.
    //
    // Parameters:
    // Center - The center point of the box.
    // Radius - The radius to extend from the center in all directions.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::setupFromCenterRadius(const fvec3& Center, float Radius) noexcept
    {
        this->m_Min = Center - fvec3(Radius);
        this->m_Max = Center + fvec3(Radius);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Sets up the box from an array of vertices by computing their bounds.
    //
    // Parameters:
    // Verts - Array of vertices to enclose.
    // NVerts - Number of vertices in the array.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    // Notes:
    // Asserts NVerts > 0 and Verts is not null.
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::setupFromVerts(const fvec3* Verts, std::int32_t NVerts) noexcept
    {
        assert(NVerts > 0 && Verts);
        this->m_Min = Verts[0];
        this->m_Max = Verts[0];
        for (std::int32_t i = 1; i < NVerts; ++i)
        {
            this->m_Min = fvec3::Min(this->m_Min, Verts[i]);
            this->m_Max = fvec3::Max(this->m_Max, Verts[i]);
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Sets the box to zero (min and max at (0,0,0)).
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::setZero(void) noexcept
    {
        this->m_Min = fvec3::fromZero();
        this->m_Max = fvec3::fromZero();
        return *this;
    }

    //------------------------------------------------------------------------------
    // Sets the box to identity (min at (-1,-1,-1), max at (1,1,1)).
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::setIdentity(void) noexcept
    {
        this->m_Min = fvec3(-1.f);
        this->m_Max = fvec3(1.f);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Expands the box by a delta vector in all directions.
    //
    // Parameters:
    // Delta - The vector to expand by (positive values increase size).
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::Inflate(const fvec3& Delta) noexcept
    {
        this->m_Min -= Delta;
        this->m_Max += Delta;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Shrinks the box by a delta vector in all directions.
    //
    // Parameters:
    // Delta - The vector to shrink by (positive values decrease size).
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::Deflate(const fvec3& Delta) noexcept
    {
        this->m_Min += Delta;
        this->m_Max -= Delta;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Translates the box by a delta vector.
    //
    // Parameters:
    // Delta - The translation vector.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::Translate(const fvec3& Delta) noexcept
    {
        this->m_Min += Delta;
        this->m_Max += Delta;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Transforms the box by a 4x4 matrix in-place.
    //
    // Parameters:
    // Matrix - The 4x4 transformation matrix.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    // Notes:
    // Transforms all 8 vertices and recomputes the bounds.
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::Transform(const fmat4& Matrix) noexcept
    {
        *this = TransformCopy(Matrix);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Computes a transformed copy of the box.
    //
    // Parameters:
    // Matrix - The 4x4 transformation matrix.
    //
    // Returns:
    // A new bounding box enclosing the transformed vertices.
    //
    // Notes:
    // Transforms all 8 vertices and recomputes the bounds.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::TransformCopy(const fmat4& Matrix) const noexcept
    {
        fvec3 verts[8];
        getVerts(verts, 8);
        fbbox_t<V> result(Matrix * verts[0]);
        for (int i = 1; i < 8; ++i)
        {
            result += Matrix * verts[i];
        }
        return result;
    }

    //------------------------------------------------------------------------------
    // Operator overloads
    //------------------------------------------------------------------------------
    //
    // Expands the box to include another box.
    //
    // Parameters:
    // Other - The other bounding box to include.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::operator+=(const fbbox_t& Other) noexcept
    {
        this->m_Min = fvec3::Min(this->m_Min, Other.m_Min);
        this->m_Max = fvec3::Max(this->m_Max, Other.m_Max);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Expands the box to include a point.
    //
    // Parameters:
    // Point - The point to include.
    //
    // Returns:
    // Reference to this box (chainable).
    //
    template <bool V>
    inline fbbox_t<V>& fbbox_t<V>::operator+=(const fvec3& Point) noexcept
    {
        this->m_Min = fvec3::Min(this->m_Min, Point);
        this->m_Max = fvec3::Max(this->m_Max, Point);
        return *this;
    }

    //------------------------------------------------------------------------------
    // Creates a new box expanded to include another box.
    //
    // Parameters:
    // Other - The other bounding box to include.
    //
    // Returns:
    // A new bounding box enclosing both boxes.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::operator+(const fbbox_t& Other) const noexcept
    {
        fbbox_t<V> result = *this;
        return result += Other;
    }

    //------------------------------------------------------------------------------
    // Creates a new box expanded to include a point.
    //
    // Parameters:
    // Point - The point to include.
    //
    // Returns:
    // A new bounding box enclosing the point.
    //
    template <bool V>
    inline fbbox_t<V> fbbox_t<V>::operator+(const fvec3& Point) const noexcept
    {
        fbbox_t<V> result = *this;
        return result += Point;
    }

    //------------------------------------------------------------------------------
    // Checks for equality with another box.
    //
    // Parameters:
    // Other - The other bounding box to compare.
    //
    // Returns:
    // True if min and max corners are equal, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::operator==(const fbbox_t& Other) const noexcept
    {
        return this->m_Min == Other.m_Min && this->m_Max == Other.m_Max;
    }

    //------------------------------------------------------------------------------
    // Checks for inequality with another box.
    //
    // Parameters:
    // Other - The other bounding box to compare.
    //
    // Returns:
    // True if min or max corners differ, false otherwise.
    //
    template <bool V>
    inline bool fbbox_t<V>::operator!=(const fbbox_t& Other) const noexcept
    {
        return !(*this == Other);
    }

    //------------------------------------------------------------------------------
    // Accesses a component by index (const).
    //
    // Parameters:
    // Index - The index (0-2 for min x,y,z; 3-5 for max x,y,z).
    //
    // Returns:
    // The component value.
    //
    // Notes:
    // Asserts Index is in [0,5].
    //
    template <bool V>
    inline float fbbox_t<V>::operator[](std::int32_t Index) const noexcept
    {
        assert(Index >= 0 && Index < 6);
        return this->m_Elements[Index];
    }

    //------------------------------------------------------------------------------
    // Accesses a component by index (mutable).
    //
    // Parameters:
    // Index - The index (0-2 for min x,y,z; 3-5 for max x,y,z).
    //
    // Returns:
    // Reference to the component.
    //
    // Notes:
    // Asserts Index is in [0,5].
    //
    template <bool V>
    inline float& fbbox_t<V>::operator[](std::int32_t Index) noexcept
    {
        assert(Index >= 0 && Index < 6);
        return this->m_Elements[Index];
    }

    //------------------------------------------------------------------------------
    // Friend operators
    //------------------------------------------------------------------------------
    //
    // Creates a new box expanded to include a point (point + box).
    //
    // Parameters:
    // Point - The point to include.
    // BBox - The bounding box to expand.
    //
    // Returns:
    // A new bounding box enclosing the point and box.
    //
    template <bool V>
    inline fbbox_t<V> operator+(const fvec3& Point, const fbbox_t<V>& BBox) noexcept
    {
        return BBox + Point;
    }
}
