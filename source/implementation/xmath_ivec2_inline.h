#pragma once
namespace xmath
{
    //------------------------------------------------------------------------------
    // ivec2
    //------------------------------------------------------------------------------
    constexpr ivec2::ivec2(int x, int y) noexcept
        : m_X(x), m_Y(y)
    {
    }

    //------------------------------------------------------------------------------
    // ivec2
    //------------------------------------------------------------------------------
    //
    // Constructs a vector with all components set to the same value.
    //
    // Parameters:
    // value - The value for all components.
    //
    constexpr ivec2::ivec2(int value) noexcept
        : m_X(value), m_Y(value)
    {
    }

    //------------------------------------------------------------------------------
    // ivec2
    //------------------------------------------------------------------------------
    //
    // Constructs a vector from a span of ints (at least 2 elements).
    //
    // Parameters:
    // Span - The span containing at least x and y values.
    //
    constexpr ivec2::ivec2(std::span<int> Span) noexcept
    {
        if (Span.size() >= 2)
        {
            m_X = Span[0];
            m_Y = Span[1];
        }
    }

    //------------------------------------------------------------------------------
    // ivec2
    //------------------------------------------------------------------------------
    //
    // Constructs a vector from an array of doubles, casting to int.
    //
    // Parameters:
    // Conversion - The array with x and y values.
    //
    constexpr ivec2::ivec2(const std::array<double, 2>& Conversion) noexcept
    {
        m_X = static_cast<int>(Conversion[0]);
        m_Y = static_cast<int>(Conversion[1]);
    }

    //------------------------------------------------------------------------------
    // ivec2
    //------------------------------------------------------------------------------
    //
    // Constructs a vector from an array of floats, casting to int.
    //
    // Parameters:
    // Conversion - The array with x and y values.
    //
    constexpr ivec2::ivec2(const std::array<float, 2>& Conversion) noexcept
    {
        m_X = static_cast<int>(Conversion[0]);
        m_Y = static_cast<int>(Conversion[1]);
    }

    //------------------------------------------------------------------------------
    // operator std::array<int, 2>
    //------------------------------------------------------------------------------
    //
    // Conversion to array of ints.
    //
    // Returns:
    // An array with x and y as ints.
    //
    constexpr ivec2::operator std::array<int, 2>(void) const noexcept
    {
        return { m_X, m_Y };
    }

    //------------------------------------------------------------------------------
    // operator std::array<float, 2>
    //------------------------------------------------------------------------------
    //
    // Conversion to array of floats.
    //
    // Returns:
    // An array with x and y as floats.
    //
    constexpr ivec2::operator std::array<float, 2>(void) const noexcept
    {
        return { static_cast<float>(m_X), static_cast<float>(m_Y) };
    }

    //------------------------------------------------------------------------------
    // operator std::string
    //------------------------------------------------------------------------------
    //
    // Conversion to string representation.
    //
    // Returns:
    // A string in "(x, y)" format.
    //
    inline ivec2::operator std::string(void) const noexcept
    {
        return ToString();
    }

    //------------------------------------------------------------------------------
    // ToString
    //------------------------------------------------------------------------------
    //
    // Returns a string representation of the vector.
    //
    // Returns:
    // A string in "(x, y)" format.
    //
    inline std::string ivec2::ToString(void) const noexcept
    {
        return std::format("({}, {})", m_X, m_Y);
    }

    //------------------------------------------------------------------------------
    // operator<<
    //------------------------------------------------------------------------------
    //
    // Overloads the stream output operator to print the vector in "(x, y)" format.
    //
    // Parameters:
    // os - The output stream.
    // vec - The vector to print.
    //
    // Returns:
    // Reference to the output stream.
    //
    inline std::ostream& operator<< (std::ostream& os, const ivec2& vec) noexcept
    {
        return os << '(' << vec.m_X << ", " << vec.m_Y << ')';
    }

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Returns a vector with all components set to 0.
    //
    constexpr ivec2 ivec2::fromZero(void) noexcept
    {
        return ivec2(0, 0);
    }

    //------------------------------------------------------------------------------
    // fromOne
    //------------------------------------------------------------------------------
    //
    // Returns a vector with all components set to 1.
    //
    constexpr ivec2 ivec2::fromOne(void) noexcept
    {
        return ivec2(1, 1);
    }

    //------------------------------------------------------------------------------
    // fromUnit
    //------------------------------------------------------------------------------
    //
    // Returns a vector with all components set to 1.
    //
    constexpr ivec2 ivec2::fromUnit(void) noexcept
    {
        return ivec2(1, 1);
    }

    //------------------------------------------------------------------------------
    // fromUp
    //------------------------------------------------------------------------------
    //
    // Returns the up direction vector (0, 1).
    //
    constexpr ivec2 ivec2::fromUp(void) noexcept
    {
        return ivec2(0, 1);
    }

    //------------------------------------------------------------------------------
    // fromDown
    //------------------------------------------------------------------------------
    //
    // Returns the down direction vector (0, -1).
    //
    constexpr ivec2 ivec2::fromDown(void) noexcept
    {
        return ivec2(0, -1);
    }

    //------------------------------------------------------------------------------
    // fromLeft
    //------------------------------------------------------------------------------
    //
    // Returns the left direction vector (-1, 0).
    //
    constexpr ivec2 ivec2::fromLeft(void) noexcept
    {
        return ivec2(-1, 0);
    }

    //------------------------------------------------------------------------------
    // fromRight
    //------------------------------------------------------------------------------
    //
    // Returns the right direction vector (1, 0).
    //
    constexpr ivec2 ivec2::fromRight(void) noexcept
    {
        return ivec2(1, 0);
    }

    //------------------------------------------------------------------------------
    // setupMax
    //------------------------------------------------------------------------------
    //
    // Sets all components to the maximum int value.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::setupMax(void) noexcept
    {
        m_X = std::numeric_limits<int>::max();
        m_Y = std::numeric_limits<int>::max();
        return *this;
    }

    //------------------------------------------------------------------------------
    // setupZero
    //------------------------------------------------------------------------------
    //
    // Sets all components to 0.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::setupZero(void) noexcept
    {
        m_X = 0;
        m_Y = 0;
        return *this;
    }

    //------------------------------------------------------------------------------
    // setup
    //------------------------------------------------------------------------------
    //
    // Sets the x and y components.
    //
    // Parameters:
    // x - The x component.
    // y - The y component.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::setup(int x, int y) noexcept
    {
        m_X = x;
        m_Y = y;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes the dot product of two vectors.
    //
    // Parameters:
    // a - The first vector.
    // b - The second vector.
    //
    // Returns:
    // The dot product.
    //
    constexpr long long ivec2::Dot(const ivec2& a, const ivec2& b) noexcept
    {
        return static_cast<long long>(a.m_X) * b.m_X + static_cast<long long>(a.m_Y) * b.m_Y;
    }

    //------------------------------------------------------------------------------
    // Min
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise minimum of two vectors.
    //
    // Parameters:
    // a - The first vector.
    // b - The second vector.
    //
    // Returns:
    // A vector with the minimum components.
    //
    constexpr ivec2 ivec2::Min(const ivec2& a, const ivec2& b) noexcept
    {
        return ivec2(std::min(a.m_X, b.m_X), std::min(a.m_Y, b.m_Y));
    }

    //------------------------------------------------------------------------------
    // Max
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise maximum of two vectors.
    //
    // Parameters:
    // a - The first vector.
    // b - The second vector.
    //
    // Returns:
    // A vector with the maximum components.
    //
    constexpr ivec2 ivec2::Max(const ivec2& a, const ivec2& b) noexcept
    {
        return ivec2(std::max(a.m_X, b.m_X), std::max(a.m_Y, b.m_Y));
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean distance between two vectors.
    //
    // Parameters:
    // a - The first vector.
    // b - The second vector.
    //
    // Returns:
    // The distance.
    //
    inline float ivec2::Distance(const ivec2& a, const ivec2& b) noexcept
    {
        float dx = static_cast<float>(a.m_X - b.m_X);
        float dy = static_cast<float>(a.m_Y - b.m_Y);
        return std::sqrt(dx * dx + dy * dy);
    }

    //------------------------------------------------------------------------------
    // Cross
    //------------------------------------------------------------------------------
    //
    // Computes the 2D cross product (determinant) of two vectors.
    //
    // Parameters:
    // a - The first vector.
    // b - The second vector.
    //
    // Returns:
    // The cross product.
    //
    constexpr long long ivec2::Cross(const ivec2& a, const ivec2& b) noexcept
    {
        return static_cast<long long>(a.m_X) * b.m_Y - static_cast<long long>(a.m_Y) * b.m_X;
    }

    //------------------------------------------------------------------------------
    // Dot
    //------------------------------------------------------------------------------
    //
    // Computes the dot product with another vector.
    //
    // Parameters:
    // a - The other vector.
    //
    // Returns:
    // The dot product.
    //
    constexpr long long ivec2::Dot(const ivec2& a) const noexcept
    {
        return Dot(*this, a);
    }

    //------------------------------------------------------------------------------
    // Min
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise minimum with another vector.
    //
    // Parameters:
    // a - The other vector.
    //
    // Returns:
    // A vector with the minimum components.
    //
    constexpr ivec2 ivec2::Min(const ivec2& a) const noexcept
    {
        return Min(*this, a);
    }

    //------------------------------------------------------------------------------
    // Max
    //------------------------------------------------------------------------------
    //
    // Computes the component-wise maximum with another vector.
    //
    // Parameters:
    // a - The other vector.
    //
    // Returns:
    // A vector with the maximum components.
    //
    constexpr ivec2 ivec2::Max(const ivec2& a) const noexcept
    {
        return Max(*this, a);
    }

    //------------------------------------------------------------------------------
    // Distance
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean distance to another vector.
    //
    // Parameters:
    // a - The other vector.
    //
    // Returns:
    // The distance.
    //
    inline float ivec2::Distance(const ivec2& a) const noexcept
    {
        return Distance(*this, a);
    }

    //------------------------------------------------------------------------------
    // Length
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean length of the vector.
    //
    // Returns:
    // The length.
    //
    inline float ivec2::Length(void) const noexcept
    {
        return std::sqrt(static_cast<float>(LengthSquared()));
    }

    //------------------------------------------------------------------------------
    // LengthSq
    //------------------------------------------------------------------------------
    //
    // Computes the squared Euclidean length of the vector.
    //
    // Returns:
    // The squared length.
    //
    constexpr long long ivec2::LengthSquared(void) const noexcept
    {
        return static_cast<long long>(m_X) * m_X + static_cast<long long>(m_Y) * m_Y;
    }

    //------------------------------------------------------------------------------
    // LengthManhattan
    //------------------------------------------------------------------------------
    //
    // Computes the Manhattan length of the vector.
    //
    // Returns:
    // The Manhattan length.
    //
    constexpr int ivec2::LengthManhattan(void) const noexcept
    {
        return xmath::Abs(m_X) + xmath::Abs(m_Y);
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if all components are within the given range.
    //
    // Parameters:
    // min - The minimum value.
    // max - The maximum value.
    //
    // Returns:
    // True if all components are in [min, max].
    //
    constexpr bool ivec2::isInRange(int min, int max) const noexcept
    {
        return m_X >= min && m_X <= max && m_Y >= min && m_Y <= max;
    }

    //------------------------------------------------------------------------------
    // Equals
    //------------------------------------------------------------------------------
    //
    // Checks if this vector equals another.
    //
    // Parameters:
    // other - The other vector.
    //
    // Returns:
    // True if equal.
    //
    constexpr bool ivec2::Equals(const ivec2& other) const noexcept
    {
        return m_X == other.m_X && m_Y == other.m_Y;
    }

    //------------------------------------------------------------------------------
    // AbsCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with absolute values.
    //
    // Returns:
    // The absolute vector.
    //
    constexpr ivec2 ivec2::AbsCopy(void) const noexcept
    {
        return ivec2(xmath::Abs(m_X), xmath::Abs(m_Y));
    }

    //------------------------------------------------------------------------------
    // Abs
    //------------------------------------------------------------------------------
    //
    // Applies absolute value to this vector.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::Abs(void) noexcept
    {
        m_X = xmath::Abs(m_X);
        m_Y = xmath::Abs(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // SignCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with sign values (1, -1, or 0).
    //
    // Returns:
    // The sign vector.
    //
    constexpr ivec2 ivec2::SignCopy(void) const noexcept
    {
        auto sign = [](int val) { return (val > 0) ? 1 : (val < 0 ? -1 : 0); };
        return ivec2(sign(m_X), sign(m_Y));
    }

    //------------------------------------------------------------------------------
    // Sign
    //------------------------------------------------------------------------------
    //
    // Applies sign to this vector (1, -1, or 0).
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::Sign(void) noexcept
    {
        auto sign = [](int& val) { val = (val > 0) ? 1 : (val < 0 ? -1 : 0); };
        sign(m_X);
        sign(m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ModCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy with modulo applied component-wise.
    //
    // Parameters:
    // divisor - The divisor.
    //
    // Returns:
    // The modulo vector.
    //
    inline ivec2 ivec2::ModCopy(int divisor) const noexcept
    {
        return ivec2(m_X % divisor, m_Y % divisor);
    }

    //------------------------------------------------------------------------------
    // Mod
    //------------------------------------------------------------------------------
    //
    // Applies modulo to this vector component-wise.
    //
    // Parameters:
    // divisor - The divisor.
    //
    // Returns:
    // Reference to this vector.
    //
    inline ivec2& ivec2::Mod(int divisor) noexcept
    {
        m_X %= divisor;
        m_Y %= divisor;
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClampCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy clamped component-wise to the range.
    //
    // Parameters:
    // min_val - The minimum value.
    // max_val - The maximum value.
    //
    // Returns:
    // The clamped vector.
    //
    inline ivec2 ivec2::ClampCopy(int min_val, int max_val) const noexcept
    {
        return ivec2(std::clamp(m_X, min_val, max_val), std::clamp(m_Y, min_val, max_val));
    }

    //------------------------------------------------------------------------------
    // Clamp
    //------------------------------------------------------------------------------
    //
    // Clamps this vector component-wise to the range.
    //
    // Parameters:
    // min_val - The minimum value.
    // max_val - The maximum value.
    //
    // Returns:
    // Reference to this vector.
    //
    inline ivec2& ivec2::Clamp(int min_val, int max_val) noexcept
    {
        m_X = std::clamp(m_X, min_val, max_val);
        m_Y = std::clamp(m_Y, min_val, max_val);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClampCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy clamped component-wise between min and max vectors.
    //
    // Parameters:
    // min - The min vector.
    // max - The max vector.
    //
    // Returns:
    // The clamped vector.
    //
    inline ivec2 ivec2::ClampCopy(const ivec2& min, const ivec2& max) const noexcept
    {
        return ivec2(std::clamp(m_X, min.m_X, max.m_X), std::clamp(m_Y, min.m_Y, max.m_Y));
    }

    //------------------------------------------------------------------------------
    // Clamp
    //------------------------------------------------------------------------------
    //
    // Clamps this vector component-wise between min and max vectors.
    //
    // Parameters:
    // min - The min vector.
    // max - The max vector.
    //
    // Returns:
    // Reference to this vector.
    //
    inline ivec2& ivec2::Clamp(const ivec2& min, const ivec2& max) noexcept
    {
        m_X = std::clamp(m_X, min.m_X, max.m_X);
        m_Y = std::clamp(m_Y, min.m_Y, max.m_Y);
        return *this;
    }

    //------------------------------------------------------------------------------
    // DistanceSquare
    //------------------------------------------------------------------------------
    //
    // Computes the squared Euclidean distance to another vector.
    //
    // Parameters:
    // v - The other vector.
    //
    // Returns:
    // The squared distance.
    //
    constexpr long long ivec2::DistanceSquare(const ivec2& v) const noexcept
    {
        long long dx = static_cast<long long>(m_X - v.m_X);
        long long dy = static_cast<long long>(m_Y - v.m_Y);
        return dx * dx + dy * dy;
    }

    //------------------------------------------------------------------------------
    // GridSnap
    //------------------------------------------------------------------------------
    //
    // Snaps this vector to the grid.
    //
    // Parameters:
    // gridX - The grid size in x.
    // gridY - The grid size in y.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::GridSnap(int gridX, int gridY) noexcept
    {
        m_X = (m_X / gridX) * gridX;
        m_Y = (m_Y / gridY) * gridY;
        return *this;
    }

    //------------------------------------------------------------------------------
    // Perp
    //------------------------------------------------------------------------------
    //
    // Returns the perpendicular vector (-y, x).
    //
    // Returns:
    // The perpendicular vector.
    //
    constexpr ivec2 ivec2::Perp(void) const noexcept
    {
        return ivec2(-m_Y, m_X);
    }

    //------------------------------------------------------------------------------
    // WhichSideOfLine
    //------------------------------------------------------------------------------
    //
    // Determines which side of the line V0-V1 this point is on.
    //
    // Parameters:
    // V0 - The start of the line.
    // V1 - The end of the line.
    //
    // Returns:
    // Positive if left, negative if right, zero if on line.
    //
    constexpr long long ivec2::WhichSideOfLine(const ivec2& V0, const ivec2& V1) const noexcept
    {
        long long dx1 = static_cast<long long>(V1.m_X - V0.m_X);
        long long dy1 = static_cast<long long>(V1.m_Y - V0.m_Y);
        long long dx2 = static_cast<long long>(m_X - V0.m_X);
        long long dy2 = static_cast<long long>(m_Y - V0.m_Y);
        return dx1 * dy2 - dy1 * dx2;
    }

    //------------------------------------------------------------------------------
    // x
    //------------------------------------------------------------------------------
    //
    // Returns the x component.
    //
    // Returns:
    // The x value.
    //
    constexpr int ivec2::x(void) const noexcept
    {
        return m_X;
    }

    //------------------------------------------------------------------------------
    // y
    //------------------------------------------------------------------------------
    //
    // Returns the y component.
    //
    // Returns:
    // The y value.
    //
    constexpr int ivec2::y(void) const noexcept
    {
        return m_Y;
    }

    //------------------------------------------------------------------------------
    // xx
    //------------------------------------------------------------------------------
    //
    // Returns a vector with (x, x).
    //
    // Returns:
    // The swizzled vector.
    //
    constexpr ivec2 ivec2::xx(void) const noexcept
    {
        return ivec2(m_X, m_X);
    }

    //------------------------------------------------------------------------------
    // xy
    //------------------------------------------------------------------------------
    //
    // Returns a vector with (x, y).
    //
    // Returns:
    // The swizzled vector.
    //
    constexpr ivec2 ivec2::xy(void) const noexcept
    {
        return ivec2(m_X, m_Y);
    }

    //------------------------------------------------------------------------------
    // yx
    //------------------------------------------------------------------------------
    //
    // Returns a vector with (y, x).
    //
    // Returns:
    // The swizzled vector.
    //
    constexpr ivec2 ivec2::yx(void) const noexcept
    {
        return ivec2(m_Y, m_X);
    }

    //------------------------------------------------------------------------------
    // yy
    //------------------------------------------------------------------------------
    //
    // Returns a vector with (y, y).
    //
    // Returns:
    // The swizzled vector.
    //
    constexpr ivec2 ivec2::yy(void) const noexcept
    {
        return ivec2(m_Y, m_Y);
    }

    //------------------------------------------------------------------------------
    // operator+
    //------------------------------------------------------------------------------
    //
    // Adds two vectors.
    //
    // Parameters:
    // other - The vector to add.
    //
    // Returns:
    // The sum vector.
    //
    constexpr ivec2 ivec2::operator+(const ivec2& other) const noexcept
    {
        return ivec2(m_X + other.m_X, m_Y + other.m_Y);
    }

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Subtracts two vectors.
    //
    // Parameters:
    // other - The vector to subtract.
    //
    // Returns:
    // The difference vector.
    //
    constexpr ivec2 ivec2::operator-(const ivec2& other) const noexcept
    {
        return ivec2(m_X - other.m_X, m_Y - other.m_Y);
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Multiplies by a scalar.
    //
    // Parameters:
    // scalar - The scalar.
    //
    // Returns:
    // The scaled vector.
    //
    constexpr ivec2 ivec2::operator*(int scalar) const noexcept
    {
        return ivec2(m_X * scalar, m_Y * scalar);
    }

    //------------------------------------------------------------------------------
    // operator/
    //------------------------------------------------------------------------------
    //
    // Divides by a scalar (integer division).
    //
    // Parameters:
    // scalar - The scalar.
    //
    // Returns:
    // The divided vector.
    //
    constexpr ivec2 ivec2::operator/(int scalar) const noexcept
    {
        return ivec2(m_X / scalar, m_Y / scalar);
    }

    //------------------------------------------------------------------------------
    // operator+=
    //------------------------------------------------------------------------------
    //
    // Adds another vector to this.
    //
    // Parameters:
    // other - The vector to add.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::operator+=(const ivec2& other) noexcept
    {
        m_X += other.m_X;
        m_Y += other.m_Y;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator-=
    //------------------------------------------------------------------------------
    //
    // Subtracts another vector from this.
    //
    // Parameters:
    // other - The vector to subtract.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::operator-=(const ivec2& other) noexcept
    {
        m_X -= other.m_X;
        m_Y -= other.m_Y;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator*=
    //------------------------------------------------------------------------------
    //
    // Multiplies this by a scalar.
    //
    // Parameters:
    // scalar - The scalar.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::operator*=(int scalar) noexcept
    {
        m_X *= scalar;
        m_Y *= scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator/=
    //------------------------------------------------------------------------------
    //
    // Divides this by a scalar (integer division).
    //
    // Parameters:
    // scalar - The scalar.
    //
    // Returns:
    // Reference to this vector.
    //
    constexpr ivec2& ivec2::operator/=(int scalar) noexcept
    {
        m_X /= scalar;
        m_Y /= scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator==
    //------------------------------------------------------------------------------
    //
    // Checks equality with another vector.
    //
    // Parameters:
    // other - The other vector.
    //
    // Returns:
    // True if equal.
    //
    constexpr bool ivec2::operator==(const ivec2& other) const noexcept
    {
        return m_X == other.m_X && m_Y == other.m_Y;
    }

    //------------------------------------------------------------------------------
    // operator!=
    //------------------------------------------------------------------------------
    //
    // Checks inequality with another vector.
    //
    // Parameters:
    // other - The other vector.
    //
    // Returns:
    // True if not equal.
    //
    constexpr bool ivec2::operator!=(const ivec2& other) const noexcept
    {
        return !(*this == other);
    }

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses a component by index.
    //
    // Parameters:
    // index - The index (0 for x, 1 for y).
    //
    // Returns:
    // The component value.
    //
    constexpr int ivec2::operator[](std::int32_t index) const noexcept
    {
        return m_Elements[static_cast<size_t>(index)];
    }

    //------------------------------------------------------------------------------
    // operator[]
    //------------------------------------------------------------------------------
    //
    // Accesses a component by index for modification.
    //
    // Parameters:
    // index - The index (0 for x, 1 for y).
    //
    // Returns:
    // Reference to the component.
    //
    constexpr int& ivec2::operator[](std::int32_t index) noexcept
    {
        return m_Elements[static_cast<size_t>(index)];
    }

    //------------------------------------------------------------------------------
    // operator*
    //------------------------------------------------------------------------------
    //
    // Multiplies a scalar by a vector.
    //
    // Parameters:
    // scalar - The scalar.
    // v - The vector.
    //
    // Returns:
    // The scaled vector.
    //
    constexpr ivec2 operator*(int scalar, const ivec2& v) noexcept
    {
        return v * scalar;
    }

    //------------------------------------------------------------------------------
    // operator-
    //------------------------------------------------------------------------------
    //
    // Unary negation of the vector.
    //
    // Parameters:
    // v - The vector.
    //
    // Returns:
    // The negated vector.
    //
    constexpr ivec2 operator-(const ivec2& v) noexcept
    {
        return ivec2(-v.m_X, -v.m_Y);
    }
}