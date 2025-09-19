#pragma once
namespace xmath
{
    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Creates a rectangle with all components set to zero.
    //
    // Returns:
    // The zero rectangle.
    //
    constexpr irect irect::fromZero(void) noexcept
    {
        return irect(ivec2(0, 0), ivec2(0, 0));
    }

    //------------------------------------------------------------------------------
    // fromIdentity
    //------------------------------------------------------------------------------
    //
    // Creates an identity rectangle, typically from (0,0) to (1,1) or similar; here assuming minimal non-empty.
    //
    // Returns:
    // The identity rectangle.
    //
    constexpr irect irect::fromIdentity(void) noexcept
    {
        return irect(ivec2(0, 0), ivec2(1, 1));
    }

    //------------------------------------------------------------------------------
    // fromPoint
    //------------------------------------------------------------------------------
    //
    // Creates a degenerate rectangle from a single point (min == max).
    //
    // Parameters:
    // Point - The point.
    //
    // Returns:
    // The rectangle.
    //
    inline irect irect::fromPoint(const ivec2& Point) noexcept
    {
        return irect(Point, Point);
    }

    //------------------------------------------------------------------------------
    // fromMinMax
    //------------------------------------------------------------------------------
    //
    // Creates a rectangle from min and max corners.
    //
    // Parameters:
    // Min - The minimum corner.
    // Max - The maximum corner.
    //
    // Returns:
    // The rectangle.
    //
    inline irect irect::fromMinMax(const ivec2& Min, const ivec2& Max) noexcept
    {
        return irect(Min, Max);
    }

    //------------------------------------------------------------------------------
    // fromCenterRadius
    //------------------------------------------------------------------------------
    //
    // Creates a rectangle from center and radius (symmetric extent).
    //
    // Parameters:
    // Center - The center point.
    // Radius - The radius.
    //
    // Returns:
    // The rectangle.
    //
    inline irect irect::fromCenterRadius(const ivec2& Center, int Radius) noexcept
    {
        return irect(Center - ivec2(Radius), Center + ivec2(Radius));
    }

    //------------------------------------------------------------------------------
    // fromVerts
    //------------------------------------------------------------------------------
    //
    // Creates a rectangle enclosing a span of vertices.
    //
    // Parameters:
    // Verts - The span of vertices.
    //
    // Returns:
    // The enclosing rectangle.
    //
    inline irect irect::fromVerts(const std::span<const ivec2> Verts) noexcept
    {
        if (Verts.empty()) return fromZero();
        ivec2 min = Verts[0], max = Verts[0];
        for (const auto& v : Verts.subspan(1)) {
            min = ivec2::Min(min, v);
            max = ivec2::Max(max, v);
        }
        return irect(min, max);
    }

    //------------------------------------------------------------------------------
    // Intersect (static, two rects)
    //------------------------------------------------------------------------------
    //
    // Checks if two rectangles intersect.
    //
    // Parameters:
    // Rect1 - First rectangle.
    // Rect2 - Second rectangle.
    //
    // Returns:
    // True if they intersect.
    //
    inline bool irect::Intersect(const irect& Rect1, const irect& Rect2) noexcept
    {
        return !(Rect1.m_Max.m_X < Rect2.m_Min.m_X || Rect1.m_Max.m_Y < Rect2.m_Min.m_Y ||
            Rect1.m_Min.m_X > Rect2.m_Max.m_X || Rect1.m_Min.m_Y > Rect2.m_Max.m_Y);
    }

    //------------------------------------------------------------------------------
    // Intersect (static, rect and point)
    //------------------------------------------------------------------------------
    //
    // Checks if a rectangle contains a point.
    //
    // Parameters:
    // Rect - The rectangle.
    // Point - The point.
    //
    // Returns:
    // True if contains.
    //
    inline bool irect::Intersect(const irect& Rect, const ivec2& Point) noexcept
    {
        return Point.m_X >= Rect.m_Min.m_X && Point.m_X <= Rect.m_Max.m_X &&
            Point.m_Y >= Rect.m_Min.m_Y && Point.m_Y <= Rect.m_Max.m_Y;
    }

    //------------------------------------------------------------------------------
    // Intersect (static, rect and line segment)
    //------------------------------------------------------------------------------
    //
    // Checks if a rectangle intersects a line segment, computing entry parameter T.
    //
    // Parameters:
    // Rect - The rectangle.
    // T - Output entry parameter (0-1 if intersects).
    // P0 - Line start.
    // P1 - Line end.
    //
    // Returns:
    // True if intersects.
    //
    inline bool irect::Intersect(const irect& Rect, float& T, const ivec2& P0, const ivec2& P1) noexcept
    {
        // Slab method for 2D AABB-line intersection.
        ivec2 dir = P1 - P0;
        if (dir.m_X == 0 && dir.m_Y == 0) return Intersect(Rect, P0);
        fvec2 inv_dir(1.0f / dir.m_X, 1.0f / dir.m_Y);
        fvec2 t_near = (fvec2(static_cast<float>(Rect.m_Min.m_X), static_cast<float>(Rect.m_Min.m_Y)) - fvec2(static_cast<float>(P0.m_X), static_cast<float>(P0.m_Y))) * inv_dir;
        fvec2 t_far =  (fvec2(static_cast<float>(Rect.m_Max.m_X), static_cast<float>(Rect.m_Max.m_Y)) - fvec2(static_cast<float>(P0.m_X), static_cast<float>(P0.m_Y))) * inv_dir;
        if (std::isnan(t_near.m_X) || std::isnan(t_near.m_Y) || std::isnan(t_far.m_X) || std::isnan(t_far.m_Y)) return false;
        if (t_near.m_X > t_far.m_X) std::swap(t_near.m_X, t_far.m_X);
        if (t_near.m_Y > t_far.m_Y) std::swap(t_near.m_Y, t_far.m_Y);
        float tmin = std::max(t_near.m_X, t_near.m_Y);
        float tmax = std::min(t_far.m_X, t_far.m_Y);
        if (tmin > tmax) return false;
        if (tmax < 0.0f) return false;
        T = std::max(0.0f, tmin);
        if (T > 1.0f) return false;
        return true;
    }

    //------------------------------------------------------------------------------
    // Intersect (static, rect and verts)
    //------------------------------------------------------------------------------
    //
    // Checks if a rectangle intersects any of the vertices.
    //
    // Parameters:
    // Rect - The rectangle.
    // Verts - The span of vertices.
    //
    // Returns:
    // True if any vertex intersects.
    //
    inline bool irect::Intersect(const irect& Rect, const std::span<const ivec2> Verts) noexcept
    {
        size_t n = Verts.size();
        if (n < 2) 
        {  // Handle points only
            for (const auto& v : Verts) 
            {
                if (Intersect(Rect, v)) return true;
            }
            return false;
        }

        // Check vertices
        for (const auto& v : Verts) 
        {
            if (Intersect(Rect, v)) return true;
        }

        // Check edges
        for (size_t i = 0; i < n; ++i) 
        {
            const ivec2& p0 = Verts[i];
            const ivec2& p1 = Verts[(i + 1) % n];
            float t;  // Dummy, since we don't need it
            if (Intersect(Rect, t, p0, p1)) return true;
        }
        return false;
    }

    //------------------------------------------------------------------------------
    // Intersection (static)
    //------------------------------------------------------------------------------
    //
    // Computes the intersection of two rectangles.
    //
    // Parameters:
    // Rect1 - First rectangle.
    // Rect2 - Second rectangle.
    //
    // Returns:
    // The intersecting rectangle (empty if none).
    //
    inline irect irect::Intersection(const irect& Rect1, const irect& Rect2) noexcept
    {
        ivec2 min = ivec2::Max(Rect1.m_Min, Rect2.m_Min);
        ivec2 max = ivec2::Min(Rect1.m_Max, Rect2.m_Max);
        if (min.m_X > max.m_X || min.m_Y > max.m_Y) return fromZero();
        return irect(min, max);
    }

    //------------------------------------------------------------------------------
    // irect (constructor from point)
    //------------------------------------------------------------------------------
    //
    // Constructs a degenerate rectangle from a point.
    //
    constexpr irect::irect(const ivec2& Point) noexcept : m_Min(Point), m_Max(Point) {}

    //------------------------------------------------------------------------------
    // irect (constructor from min max)
    //------------------------------------------------------------------------------
    //
    // Constructs from min and max corners.
    //
    constexpr irect::irect(const ivec2& Min, const ivec2& Max) noexcept : m_Min(Min), m_Max(Max) {}

    //------------------------------------------------------------------------------
    // irect (constructor from verts)
    //------------------------------------------------------------------------------
    //
    // Constructs enclosing verts.
    //
    inline irect::irect(const std::span<const ivec2> Verts) noexcept
    {
        *this = fromVerts(Verts);
    }

    //------------------------------------------------------------------------------
    // irect (constructor from center radius)
    //------------------------------------------------------------------------------
    //
    // Constructs from center and radius.
    //
    constexpr irect::irect(const ivec2& Center, int Radius) noexcept : m_Min(Center - ivec2(Radius)), m_Max(Center + ivec2(Radius)) {}

    //------------------------------------------------------------------------------
    // irect (constructor from array)
    //------------------------------------------------------------------------------
    //
    // Constructs from array [min.m_X, min.m_Y, max.m_X, max.m_Y].
    //
    constexpr irect::irect(const std::array<int, 4>& Conversion) noexcept : m_Elements(Conversion) {}

    //------------------------------------------------------------------------------
    // irect (constructor from L T R B)
    //------------------------------------------------------------------------------
    //
    // Constructs from left, top, right, bottom.
    //
    constexpr irect::irect(int L, int T, int R, int B) noexcept : m_Min(L, T), m_Max(R, B) {}

    //------------------------------------------------------------------------------
    // operator std::array<int,4>
    //------------------------------------------------------------------------------
    //
    // Converts to array [min.m_X, min.m_Y, max.m_X, max.m_Y].
    //
    constexpr irect::operator std::array<int, 4>(void) const noexcept
    {
        return { m_Min.m_X, m_Min.m_Y, m_Max.m_X, m_Max.m_Y };
    }

    //------------------------------------------------------------------------------
    // operator std::string
    //------------------------------------------------------------------------------
    //
    // Converts to string representation.
    //
    inline irect::operator std::string(void) const noexcept
    {
        return toString();
    }

    //------------------------------------------------------------------------------
    // toString
    //------------------------------------------------------------------------------
    //
    // Gets string representation, e.g., "[min max]".
    //
    // Returns:
    // The string.
    //
    inline std::string irect::toString(void) const noexcept
    {
        return std::format("[{}, {}]", m_Min.ToString(), m_Max.ToString());
    }

    //------------------------------------------------------------------------------
    // operator<<
    //------------------------------------------------------------------------------
    //
    // Streams to ostream.
    //
    // Parameters:
    // Os - The stream.
    // Rect - The rectangle.
    //
    // Returns:
    // The stream.
    //
    inline std::ostream& operator<<(std::ostream& Os, const irect& Rect) noexcept
    {
        return Os << Rect.toString();
    }

    //------------------------------------------------------------------------------
    // Intersect (member, other rect)
    //------------------------------------------------------------------------------
    //
    // Checks intersection with another rectangle.
    //
    // Parameters:
    // Other - The other rectangle.
    //
    // Returns:
    // True if intersects.
    //
    inline bool irect::Intersect(const irect& Other) const noexcept
    {
        return Intersect(*this, Other);
    }

    //------------------------------------------------------------------------------
    // Intersect (member, point)
    //------------------------------------------------------------------------------
    //
    // Checks if contains point.
    //
    // Parameters:
    // Point - The point.
    //
    // Returns:
    // True if contains.
    //
    inline bool irect::Intersect(const ivec2& Point) const noexcept
    {
        return Intersect(*this, Point);
    }

    //------------------------------------------------------------------------------
    // Intersect (member, line)
    //------------------------------------------------------------------------------
    //
    // Checks intersection with line segment.
    //
    // Parameters:
    // T - Output entry T.
    // P0 - Start.
    // P1 - End.
    //
    // Returns:
    // True if intersects.
    //
    inline bool irect::Intersect(float& T, const ivec2& P0, const ivec2& P1) const noexcept
    {
        return Intersect(*this, T, P0, P1);
    }

    //------------------------------------------------------------------------------
    // Intersect (member, verts)
    //------------------------------------------------------------------------------
    //
    // Checks intersection with verts.
    //
    // Parameters:
    // Verts - The verts.
    //
    // Returns:
    // True if any intersects.
    //
    inline bool irect::Intersect(const std::span<ivec2> Verts) const noexcept
    {
        return Intersect(*this, Verts);
    }

    //------------------------------------------------------------------------------
    // Intersection (member)
    //------------------------------------------------------------------------------
    //
    // Computes intersection with other.
    //
    // Parameters:
    // Other - The other.
    //
    // Returns:
    // The intersection.
    //
    inline irect irect::Intersection(const irect& Other) const noexcept
    {
        return Intersection(*this, Other);
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if all components are finite (always true for ints).
    //
    // Returns:
    // True.
    //
    inline bool irect::isFinite(void) const noexcept
    {
        return true;  // Ints are always finite.
    }

    //------------------------------------------------------------------------------
    // isValid
    //------------------------------------------------------------------------------
    //
    // Checks if min <= max.
    //
    // Returns:
    // True if valid.
    //
    inline bool irect::isValid(void) const noexcept
    {
        return m_Min.m_X <= m_Max.m_X && m_Min.m_Y <= m_Max.m_Y;
    }

    //------------------------------------------------------------------------------
    // getSize
    //------------------------------------------------------------------------------
    //
    // Gets size (max - min).
    //
    // Returns:
    // The size.
    //
    inline ivec2 irect::getSize(void) const noexcept
    {
        return m_Max - m_Min;
    }

    //------------------------------------------------------------------------------
    // getCenter
    //------------------------------------------------------------------------------
    //
    // Gets center ((min + max)/2).
    //
    // Returns:
    // The center.
    //
    inline ivec2 irect::getCenter(void) const noexcept
    {
        return (m_Min + m_Max) / 2;
    }

    //------------------------------------------------------------------------------
    // getRadius
    //------------------------------------------------------------------------------
    //
    // Gets half-diagonal length.
    //
    // Returns:
    // The radius.
    //
    inline float irect::getRadius(void) const noexcept
    {
        return getSize().Length() / 2;
    }

    //------------------------------------------------------------------------------
    // getRadiusSquared
    //------------------------------------------------------------------------------
    //
    // Gets squared half-diagonal.
    //
    // Returns:
    // The squared radius.
    //
    inline long long irect::getRadiusSquared(void) const noexcept
    {
        return getSize().LengthSquared() / 4;
    }

    //------------------------------------------------------------------------------
    // getArea
    //------------------------------------------------------------------------------
    //
    // Gets area (width * height).
    //
    // Returns:
    // The area.
    //
    inline int irect::getArea(void) const noexcept
    {
        ivec2 size = getSize();
        return size.m_X * size.m_Y;
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if all components in [Min, Max].
    //
    // Parameters:
    // Min - Lower bound.
    // Max - Upper bound.
    //
    // Returns:
    // True if in range.
    //
    inline bool irect::isInRange(int Min, int Max) const noexcept
    {
        return m_Min.m_X >= Min && m_Min.m_Y >= Min && m_Max.m_X <= Max && m_Max.m_Y <= Max;
    }

    //------------------------------------------------------------------------------
    // Contains
    //------------------------------------------------------------------------------
    //
    // Checks if contains another rectangle.
    //
    // Parameters:
    // Other - The other.
    //
    // Returns:
    // True if contains.
    //
    inline bool irect::Contains(const irect& Other) const noexcept
    {
        return Other.m_Min.m_X >= m_Min.m_X && Other.m_Min.m_Y >= m_Min.m_Y &&
               Other.m_Max.m_X <= m_Max.m_X && Other.m_Max.m_Y <= m_Max.m_Y;
    }

    //------------------------------------------------------------------------------
    // ContainsPoint
    //------------------------------------------------------------------------------
    //
    // Checks if contains point.
    //
    // Parameters:
    // Point - The point.
    //
    // Returns:
    // True if contains.
    //
    inline bool irect::ContainsPoint(const ivec2& Point) const noexcept
    {
        return Intersect(Point);
    }

    //------------------------------------------------------------------------------
    // getClosestPoint
    //------------------------------------------------------------------------------
    //
    // Gets closest point on rectangle to given point.
    //
    // Parameters:
    // Point - The point.
    //
    // Returns:
    // The closest point.
    //
    inline ivec2 irect::getClosestPoint(const ivec2& Point) const noexcept
    {
        return ivec2(std::clamp(Point.m_X, m_Min.m_X, m_Max.m_X), std::clamp(Point.m_Y, m_Min.m_Y, m_Max.m_Y));
    }

    //------------------------------------------------------------------------------
    // getVerts
    //------------------------------------------------------------------------------
    //
    // Fills span with 4 corners (bottom-left, bottom-right, top-right, top-left).
    //
    // Parameters:
    // Dst - Span of at least 4 ivec2.
    //
    inline void irect::getVerts(std::span<ivec2> Dst) const noexcept
    {
        if (Dst.size() < 4) return;
        Dst[0] = ivec2(m_Min.m_X, m_Min.m_Y);
        Dst[1] = ivec2(m_Max.m_X, m_Min.m_Y);
        Dst[2] = ivec2(m_Max.m_X, m_Max.m_Y);
        Dst[3] = ivec2(m_Min.m_X, m_Max.m_Y);
    }

    //------------------------------------------------------------------------------
    // getWidth
    //------------------------------------------------------------------------------
    //
    // Gets width (max.m_X - min.m_X).
    //
    // Returns:
    // The width.
    //
    inline int irect::getWidth(void) const noexcept
    {
        return m_Max.m_X - m_Min.m_X;
    }

    //------------------------------------------------------------------------------
    // getHeight
    //------------------------------------------------------------------------------
    //
    // Gets height (max.m_Y - min.m_Y).
    //
    // Returns:
    // The height.
    //
    inline int irect::getHeight(void) const noexcept
    {
        return m_Max.m_Y - m_Min.m_Y;
    }

    //------------------------------------------------------------------------------
    // Difference
    //------------------------------------------------------------------------------
    //
    // Computes difference metric, e.g., symmetric area difference.
    //
    // Parameters:
    // Other - The other.
    //
    // Returns:
    // The difference.
    //
    inline float irect::Difference(const irect& Other) const noexcept
    {
        irect inter = Intersection(Other);
        int area1 = getArea() - inter.getArea();
        int area2 = Other.getArea() - inter.getArea();
        return static_cast<float>(area1 + area2) / std::max(1, getArea() + Other.getArea() - inter.getArea());
    }

    //------------------------------------------------------------------------------
    // isEmpty
    //------------------------------------------------------------------------------
    //
    // Checks if width or height <=0.
    //
    // Returns:
    // True if empty.
    //
    inline bool irect::isEmpty(void) const noexcept
    {
        return getWidth() <= 0 || getHeight() <= 0;
    }

    //------------------------------------------------------------------------------
    // setupFromPoint
    //------------------------------------------------------------------------------
    //
    // Sets to degenerate from point.
    //
    // Parameters:
    // Point - The point.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setupFromPoint(const ivec2& Point) noexcept
    {
        m_Min = m_Max = Point;
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupFromMinMax
    //------------------------------------------------------------------------------
    //
    // Sets from min max.
    //
    // Parameters:
    // Min - Min.
    // Max - Max.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setupFromMinMax(const ivec2& Min, const ivec2& Max) noexcept
    {
        m_Min = Min;
        m_Max = Max;
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupFromCenterRadius
    //------------------------------------------------------------------------------
    //
    // Sets from center radius.
    //
    // Parameters:
    // Center - Center.
    // Radius - Radius.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setupFromCenterRadius(const ivec2& Center, int Radius) noexcept
    {
        m_Min = Center - ivec2(Radius);
        m_Max = Center + ivec2(Radius);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupFromVerts
    //------------------------------------------------------------------------------
    //
    // Sets to enclose verts.
    //
    // Parameters:
    // Verts - Verts.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setupFromVerts(const std::span<const ivec2> Verts) noexcept
    {
        *this = fromVerts(Verts);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setup
    //------------------------------------------------------------------------------
    //
    // Sets from L T R B.
    //
    // Parameters:
    // L - Left.
    // T - Top.
    // R - Right.
    // B - Bottom.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setup(int L, int T, int R, int B) noexcept
    {
        m_Min = ivec2(L, T);
        m_Max = ivec2(R, B);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setZero
    //------------------------------------------------------------------------------
    //
    // Sets to zero.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setZero(void) noexcept
    {
        m_Min = m_Max = ivec2(0);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setIdentity
    //------------------------------------------------------------------------------
    //
    // Sets to identity.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setIdentity(void) noexcept
    {
        *this = fromIdentity();
        return *this;
    }

    //------------------------------------------------------------------------------
    // AddVerts
    //------------------------------------------------------------------------------
    //
    // Expands to include verts.
    //
    // Parameters:
    // Verts - Verts.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::AddVerts(const std::span<const ivec2> Verts) noexcept
    {
        for (const auto& v : Verts) {
            m_Min = ivec2::Min(m_Min, v);
            m_Max = ivec2::Max(m_Max, v);
        }
        return *this;
    }

    //------------------------------------------------------------------------------
    // Inflate
    //------------------------------------------------------------------------------
    //
    // Expands by delta (negative deflates).
    //
    // Parameters:
    // Delta - Amount to inflate.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::Inflate(const ivec2& Delta) noexcept
    {
        m_Min -= Delta;
        m_Max += Delta;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Deflate
    //------------------------------------------------------------------------------
    //
    // Deflates by delta.
    //
    // Parameters:
    // Delta - Amount to deflate.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::Deflate(const ivec2& Delta) noexcept
    {
        m_Min += Delta;
        m_Max -= Delta;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Translate
    //------------------------------------------------------------------------------
    //
    // Translates by delta.
    //
    // Parameters:
    // Delta - Translation.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::Translate(const ivec2& Delta) noexcept
    {
        m_Min += Delta;
        m_Max += Delta;
        return *this;
    }

    //------------------------------------------------------------------------------
    // setWidth
    //------------------------------------------------------------------------------
    //
    // Sets width, keeping min, adjusting max.
    //
    // Parameters:
    // W - New width.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setWidth(int W) noexcept
    {
        m_Max.m_X = m_Min.m_X + W;
        return *this;
    }

    //------------------------------------------------------------------------------
    // setHeight
    //------------------------------------------------------------------------------
    //
    // Sets height, keeping min, adjusting max.
    //
    // Parameters:
    // H - New height.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setHeight(int H) noexcept
    {
        m_Max.m_Y = m_Min.m_Y + H;
        return *this;
    }

    //------------------------------------------------------------------------------
    // setSize
    //------------------------------------------------------------------------------
    //
    // Sets size, keeping min, adjusting max.
    //
    // Parameters:
    // W - Width.
    // H - Height.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setSize(int W, int H) noexcept
    {
        m_Max = m_Min + ivec2(W, H);
        return *this;
    }

    //------------------------------------------------------------------------------
    // setMax
    //------------------------------------------------------------------------------
    //
    // Sets to maximum bounds (INT_MIN to INT_MAX).
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::setMax(void) noexcept
    {
        m_Min = ivec2(std::numeric_limits<int>::min());
        m_Max = ivec2(std::numeric_limits<int>::max());
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator+= (rect)
    //------------------------------------------------------------------------------
    //
    // Unions with other.
    //
    // Parameters:
    // Other - Other.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::operator+=(const irect& Other) noexcept
    {
        m_Min = ivec2::Min(m_Min, Other.m_Min);
        m_Max = ivec2::Max(m_Max, Other.m_Max);
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator+= (point)
    //------------------------------------------------------------------------------
    //
    // Adds point (expands).
    //
    // Parameters:
    // Point - Point.
    //
    // Returns:
    // Reference to self.
    //
    inline irect& irect::operator+=(const ivec2& Point) noexcept
    {
        m_Min = ivec2::Min(m_Min, Point);
        m_Max = ivec2::Max(m_Max, Point);
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator+ (rect)
    //------------------------------------------------------------------------------
    //
    // Returns union.
    //
    // Parameters:
    // Other - Other.
    //
    // Returns:
    // The union.
    //
    inline irect irect::operator+(const irect& Other) const noexcept
    {
        return irect(ivec2::Min(m_Min, Other.m_Min), ivec2::Max(m_Max, Other.m_Max));
    }

    //------------------------------------------------------------------------------
    // operator+ (point)
    //------------------------------------------------------------------------------
    //
    // Returns expanded by point.
    //
    // Parameters:
    // Point - Point.
    //
    // Returns:
    // The result.
    //
    inline irect irect::operator+(const ivec2& Point) const noexcept
    {
        return irect(ivec2::Min(m_Min, Point), ivec2::Max(m_Max, Point));
    }

    //------------------------------------------------------------------------------
    // operator==
    //------------------------------------------------------------------------------
    //
    // Checks equality.
    //
    // Parameters:
    // Other - Other.
    //
    // Returns:
    // True if equal.
    //
    inline bool irect::operator==(const irect& Other) const noexcept
    {
        return m_Min == Other.m_Min && m_Max == Other.m_Max;
    }

    //------------------------------------------------------------------------------
    // operator!=
    //------------------------------------------------------------------------------
    //
    // Checks inequality.
    //
    // Parameters:
    // Other - Other.
    //
    // Returns:
    // True if not equal.
    //
    inline bool irect::operator!=(const irect& Other) const noexcept
    {
        return !(*this == Other);
    }

    //------------------------------------------------------------------------------
    // operator[] (const)
    //------------------------------------------------------------------------------
    //
    // Accesses component (0:min.m_X,1:min.m_Y,2:max.m_X,3:max.m_Y).
    //
    // Parameters:
    // Index - Index.
    //
    // Returns:
    // The value.
    //
    inline int irect::operator[](int Index) const noexcept
    {
        return m_Elements[Index];
    }

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses component (mutable).
    //
    // Parameters:
    // Index - Index.
    //
    // Returns:
    // Reference.
    //
    inline int& irect::operator[](int Index) noexcept
    {
        return m_Elements[Index];
    }

    //------------------------------------------------------------------------------
    // operator+ (friend, point + rect)
    //------------------------------------------------------------------------------
    //
    // Returns rect expanded by point.
    //
    // Parameters:
    // Point - Point.
    // Rect - Rect.
    //
    // Returns:
    // The result.
    //
    irect operator+(const ivec2& Point, const irect& Rect) noexcept
    {
        return Rect + Point;
    }
}