namespace xmath
{
    //------------------------------------------------------------------------------
    // Functions Implementations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // ModAngle
    //------------------------------------------------------------------------------
    //
    // Returns the equivalent angle in the range [0, 2PI) radians (or [0, 360) degrees).
    //
    // Params:
    //  Angle - Input angle in radians (can be negative or larger than 2PI).
    //
    // Returns:
    //  Normalized angle >= 0 and < 2PI.
    //
    // Notes:
    //  Uses std::fmod for modulo operation; handles negative values by adding full circle.
    //  Assumes input is finite; asserts validity to prevent NaN/inf propagation.
    //
    template<std::floating_point T>
    inline radian_t<T> ModAngle(radian_t<T> Angle) noexcept
    {
        assert(isValid(Angle.m_Value));
        constexpr radian_t<T> fullCircle{ std::numbers::pi_v<T> *T{ 2 } };
        Angle.m_Value = std::fmod(Angle.m_Value, fullCircle.m_Value);
        if (Angle.m_Value < T{ 0 }) Angle.m_Value += fullCircle.m_Value;
        return Angle;
    }

    //------------------------------------------------------------------------------
    // ModAngle2
    //------------------------------------------------------------------------------
    //
    // Returns the equivalent angle in the range [-PI, PI) radians (or [-180, 180) degrees).
    //
    // Params:
    //  Angle - Input angle in radians.
    //
    // Returns:
    //  Normalized angle >= -PI and < PI.
    //
    // Notes:
    //  Builds on ModAngle by shifting range; useful for signed angle differences.
    //  Assumes input is finite; asserts validity.
    //
    template<std::floating_point T>
    inline radian_t<T> ModAngle2(radian_t<T> Angle) noexcept
    {
        assert(isValid(Angle.m_Value));
        Angle = ModAngle(Angle + radian_t<T>{ std::numbers::pi_v<T> });
        return Angle - radian_t<T>{ std::numbers::pi_v<T> };
    }

    //------------------------------------------------------------------------------
    // MinAngleDiff
    //------------------------------------------------------------------------------
    //
    // Computes the smallest angular difference between two angles.
    //
    // Params:
    //  Angle1 - First angle in radians.
    //  Angle2 - Second angle in radians.
    //
    // Returns:
    //  Difference in [-PI, PI) radians; positive if Angle1 > Angle2 in shortest path.
    //
    // Notes:
    //  Uses ModAngle2 for normalization; ideal for interpolation or comparisons.
    //  Assumes inputs are finite; asserts validity.
    //
    template<std::floating_point T>
    inline radian_t<T> MinAngleDiff(radian_t<T> Angle1, radian_t<T> Angle2) noexcept
    {
        assert(isValid(Angle1.m_Value) && isValid(Angle2.m_Value));
        return ModAngle2(Angle1 - Angle2);
    }

    //------------------------------------------------------------------------------
    // LerpAngle
    //------------------------------------------------------------------------------
    //
    // Linearly interpolates between two angles using the shortest path.
    //
    // Params:
    //  t - Interpolation factor [0, 1].
    //  Angle1 - Starting angle in radians.
    //  Angle2 - Ending angle in radians.
    //
    // Returns:
    //  Interpolated angle along shortest arc.
    //
    // Notes:
    //  Uses MinAngleDiff/ModAngle for correct wrapping; t clamped implicitly via lerp.
    //  Assumes inputs finite and t in [0,1]; asserts angle validity.
    //
    template<std::floating_point T>
    inline radian_t<T> LerpAngle(T t, radian_t<T> Angle1, radian_t<T> Angle2) noexcept
    {
        assert(isValid(Angle1.m_Value) && isValid(Angle2.m_Value) && isInRange(t, T{ 0 }, T{ 1 }));
        return ModAngle(Angle1 + t * ModAngle(Angle2 - Angle1));
    }

    //------------------------------------------------------------------------------
    // SinCos
    //------------------------------------------------------------------------------
    //
    // Computes sine and cosine of an angle in one call (optimized).
    //
    // Params:
    //  Angle - Input angle in radians.
    //  S - Output sine value.
    //  C - Output cosine value.
    //
    // Notes:
    //  Uses std::sincos for efficiency; better than separate sin/cos calls.
    //  Assumes Angle finite; asserts validity; no return value.
    //
    template<std::floating_point T>
    inline void SinCos(radian_t<T> Angle, T& S, T& C) noexcept
    {
        assert(isValid(Angle.m_Value));
        S = Sin(Angle);
        C = Cos(Angle);
    }

    //------------------------------------------------------------------------------
    // Sin
    //------------------------------------------------------------------------------
    //
    // Computes the sine of an angle.
    //
    // Params:
    //  m_X - Input angle in radians.
    //
    // Returns:
    //  Sine value in [-1, 1].
    //
    // Notes:
    //  Range: Full real line; periodic every 2PI.
    //  Assumes finite input; asserts validity.
    //
    template<std::floating_point T>
    inline T Sin(radian_t<T> x) noexcept
    {
        assert(isValid(x.m_Value));
        return std::sin(x.m_Value);
    }

    //------------------------------------------------------------------------------
    // Cos
    //------------------------------------------------------------------------------
    //
    // Computes the cosine of an angle.
    //
    // Params:
    //  m_X - Input angle in radians.
    //
    // Returns:
    //  Cosine value in [-1, 1].
    //
    // Notes:
    //  Range: Full real line; periodic every 2PI.
    //  Assumes finite input; asserts validity.
    //
    template<std::floating_point T>
    inline T Cos(radian_t<T> x) noexcept
    {
        assert(isValid(x.m_Value));
        return std::cos(x.m_Value);
    }

    //------------------------------------------------------------------------------
    // Tan
    //------------------------------------------------------------------------------
    //
    // Computes the tangent of an angle.
    //
    // Params:
    //  m_X - Input angle in radians.
    //
    // Returns:
    //  Tangent value (undefined at odd multiples of PI/2).
    //
    // Notes:
    //  Range: Full real line except discontinuities; periodic every PI.
    //  Assumes finite input and not at poles; asserts validity.
    //
    template<std::floating_point T>
    inline T Tan(radian_t<T> x) noexcept
    {
        assert(isValid(x.m_Value));
        return std::tan(x.m_Value);
    }

    //------------------------------------------------------------------------------
    // Asin
    //------------------------------------------------------------------------------
    //
    // Computes the arcsine (inverse sine) of a value.
    //
    // Params:
    //  m_X - Input in [-1, 1].
    //
    // Returns:
    //  Angle in [-PI/2, PI/2] radians.
    //
    // Notes:
    //  Domain: [-1, 1]; NaN outside.
    //  Assumes m_X in domain; asserts validity and range.
    //
    template<std::floating_point T>
    inline radian_t<T> Asin(T x) noexcept
    {
        assert(isValid(x) && isInRange(x, T{ -1 }, T{ 1 }));
        return radian_t<T>{ std::asin(x) };
    }

    //------------------------------------------------------------------------------
    // Atan
    //------------------------------------------------------------------------------
    //
    // Computes the arctangent (inverse tan) of a value.
    //
    // Params:
    //  m_X - Input in [-1, 1].
    //
    // Returns:
    //  Angle in [-PI/2, PI/2] radians.
    //
    // Notes:
    //  Domain: [-1, 1]; NaN outside.
    //  Assumes m_X in domain; asserts validity and range.
    //
    template<std::floating_point T>
    inline radian_t<T> Atan(T x) noexcept
    {
        assert(isValid(x));
        return radian_t<T>{ std::atan(x) };
    }

    //------------------------------------------------------------------------------
    // Atan2
    //------------------------------------------------------------------------------
    //
    // atan2(y, x) returns the angle angle in radians between the positive x-axis and the point (x, y). The result is in the range [-pi,pi]
    //
    // Params:
    //  y: The ordinate (vertical component) of the point.
    //  x: The abscissa (horizontal component) of the point.
//
    // Returns:
    //  Angle in [-PI/2, PI/2] radians.
    //
    // Notes:
    //  Domain: [-1, 1]; NaN outside.
    //  Assumes m_X in domain; asserts validity and range.
    //
    template<std::floating_point T>
    inline radian_t<T> Atan2(T y, T x) noexcept
    {
        assert(isValid(x) && isValid(y));
        return radian_t<T>{ std::atan2(y, x) };
    }

    //------------------------------------------------------------------------------
    // Acos
    //------------------------------------------------------------------------------
    //
    // Computes the arccosine (inverse cosine) of a value.
    //
    // Params:
    //  m_X - Input in [-1, 1].
    //
    // Returns:
    //  Angle in [0, PI] radians.
    //
    // Notes:
    //  Domain: [-1, 1]; NaN outside.
    //  Assumes m_X in domain; asserts validity and range.
    //
    template<std::floating_point T>
    inline radian_t<T> Acos(T x) noexcept
    {
        assert(isValid(x) && isInRange(x, T{ -1 }, T{ 1 }));
        return radian_t<T>{ std::acos(x) };
    }
}