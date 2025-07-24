namespace xmath
{
    //------------------------------------------------------------------------------
    // Functions Implementations
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Exp
    //------------------------------------------------------------------------------
    //
    // Computes e raised to the power of m_X (exponential).
    //
    // Params:
    //  m_X - Input exponent.
    //
    // Returns:
    //  e^m_X (>0 for all real m_X).
    //
    // Notes:
    //  Range: Full real line; always positive.
    //  Assumes finite input; asserts validity.
    //
    template<std::floating_point T>
    inline T Exp(T x) noexcept
    {
        assert(isValid(x));
        return std::exp(x);
    }

    //------------------------------------------------------------------------------
    // Pow
    //------------------------------------------------------------------------------
    //
    // Computes a raised to the power of b.
    //
    // Params:
    //  a - Base (non-negative if b fractional).
    //  b - Exponent.
    //
    // Returns:
    //  a^b.
    //
    // Notes:
    //  Domain: a > 0 for non-integer b; handles negatives for integer b.
    //  Assumes valid inputs; asserts validity.
    //
    template<std::floating_point T>
    inline T Pow(T a, T b) noexcept
    {
        assert(isValid(a) && isValid(b));
        return std::pow(a, b);
    }

    //------------------------------------------------------------------------------
    // FMod
    //------------------------------------------------------------------------------
    //
    // Computes floating-point modulo (remainder) of m_X / m_Y.
    //
    // Params:
    //  m_X - Dividend.
    //  m_Y - Divisor (non-zero).
    //
    // Returns:
    //  Remainder with sign of m_X.
    //
    // Notes:
    //  Differs from % for negatives; m_Y != 0.
    //  Assumes finite/non-zero m_Y; asserts validity.
    //
    template<std::floating_point T>
    inline T FMod(T x, T y) noexcept
    {
        assert(isValid(x) && isValid(y) && y != T{ 0 });
        return std::fmod(x, y);
    }

    //------------------------------------------------------------------------------
    // ModFX
    //------------------------------------------------------------------------------
    //
    // Separates floating-point number into fractional and integer parts.
    //
    // Params:
    //  m_X - Input number.
    //  m_Y - Output integer part.
    //
    // Returns:
    //  Fractional part (sign of m_X, magnitude <1).
    //
    // Notes:
    //  m_Y gets floor(m_X) for positive, ceil(m_X) for negative (per std::modf).
    //  Assumes finite m_X; asserts validity.
    //
    template<std::floating_point T>
    inline T ModFX(T x, T& y) noexcept
    {
        assert(isValid(x));
        return std::modf(x, &y);
    }

    //------------------------------------------------------------------------------
    // Log
    //------------------------------------------------------------------------------
    //
    // Computes natural logarithm (base e) of m_X.
    //
    // Params:
    //  m_X - Input > 0.
    //
    // Returns:
    //  ln(m_X).
    //
    // Notes:
    //  Domain: (0, +inf); -inf at 0, NaN <=0.
    //  Assumes m_X > 0 and finite; asserts validity and positivity.
    //
    template<std::floating_point T>
    inline T Log(T x) noexcept
    {
        assert(isValid(x) && x > T{ 0 });
        return std::log(x);
    }

    //------------------------------------------------------------------------------
    // Log2
    //------------------------------------------------------------------------------
    //
    // Computes logarithm base 2 of m_X.
    //
    // Params:
    //  m_X - Input > 0.
    //
    // Returns:
    //  log2(m_X).
    //
    // Notes:
    //  Domain: (0, +inf); uses std::log2 for efficiency.
    //  Assumes m_X > 0 and finite; asserts validity and positivity.
    //
    template<std::floating_point T>
    inline T Log2(T x) noexcept
    {
        assert(isValid(x) && x > T{ 0 });
        return std::log2(x);
    }

    //------------------------------------------------------------------------------
    // Log10
    //------------------------------------------------------------------------------
    //
    // Computes logarithm base 10 of m_X.
    //
    // Params:
    //  m_X - Input > 0.
    //
    // Returns:
    //  log10(m_X).
    //
    // Notes:
    //  Domain: (0, +inf); uses std::log10.
    //  Assumes m_X > 0 and finite; asserts validity and positivity.
    //
    template<std::floating_point T>
    inline T Log10(T x) noexcept
    {
        assert(isValid(x) && x > T{ 0 });
        return std::log10(x);
    }

    //------------------------------------------------------------------------------
    // i2f
    //------------------------------------------------------------------------------
    //
    // Converts integer-like to float (trivial cast).
    //
    // Params:
    //  i - Input value.
    //
    // Returns:
    //  Static cast to T (floating-point).
    //
    // Notes:
    //  No range checks; assumes representable.
    //  No assert needed as trivial.
    //
    template<std::floating_point T>
    constexpr T i2f(T i) noexcept { return i; }

    //------------------------------------------------------------------------------
    // f2i
    //------------------------------------------------------------------------------
    //
    // Converts float to int32_t (truncates).
    //
    // Params:
    //  f - Input float.
    //
    // Returns:
    //  Truncated integer.
    //
    // Notes:
    //  Uses static_cast; no rounding.
    //  Assumes finite and in int32 range; asserts validity.
    //
    template<std::floating_point T>
    constexpr std::int32_t f2i(T f) noexcept
    {
        assert(isValid(f));
        return static_cast<std::int32_t>(f);
    }

    //------------------------------------------------------------------------------
    // FSel
    //------------------------------------------------------------------------------
    //
    // Branchless select: b if a >= 0, else c.
    //
    // Params:
    //  a - Condition value.
    //  b - Value if true.
    //  c - Value if false.
    //
    // Returns:
    //  Selected value.
    //
    // Notes:
    //  Useful for perf in hot paths; assumes finite a.
    //  Asserts validity of a.
    //
    template <typename T>
    constexpr T FSel(T a, T b, T c) noexcept
    {
        assert(isValid(a));
        return a >= T{ 0 } ? b : c;
    }

    //------------------------------------------------------------------------------
    // Sqr
    //------------------------------------------------------------------------------
    //
    // Computes square of m_X.
    //
    // Params:
    //  m_X - Input value.
    //
    // Returns:
    //  m_X * m_X.
    //
    // Notes:
    //  Overflow possible for large m_X; no checks.
    //  Assumes finite m_X; asserts validity.
    //
    template <typename T>
    constexpr T Sqr(T x) noexcept
    {
        assert(isValid(x));
        return x * x;
    }

    //------------------------------------------------------------------------------
    // Sqrt
    //------------------------------------------------------------------------------
    //
    // Computes square root of m_X.
    //
    // Params:
    //  m_X - Non-negative input.
    //
    // Returns:
    //  sqrt(m_X) >= 0.
    //
    // Notes:
    //  Domain: [0, +inf); NaN for <0.
    //  Assumes m_X >= 0 and finite; asserts validity and non-negative.
    //
    template<std::floating_point T>
    inline T Sqrt(T x) noexcept
    {
        assert(isValid(x) && x >= T{ 0 });
        return std::sqrt(x);
    }

    //------------------------------------------------------------------------------
    // InvSqrt
    //------------------------------------------------------------------------------
    //
    // Computes approximate inverse square root (1/sqrt(m_X)).
    //
    // Params:
    //  m_X - Positive input >0.
    //
    // Returns:
    //  1/sqrt(m_X).
    //
    // Notes:
    //  Fast approx via bit hack + Newton iterations; accuracy ~1e-5.
    //  Assumes m_X >0 and finite; asserts validity and positivity.
    //
    template<std::floating_point T>
    inline T InvSqrt(T x) noexcept
    {
        assert(isValid(x) && x > T{ 0 });
        const T xhalf = T{ 0.5 } *x;
        auto i = std::bit_cast<std::conditional_t<std::is_same_v<T, float>, std::int32_t, std::int64_t>>(x);
        i = (std::is_same_v<T, float> ? 0x5f3759df : 0x5fe6eb50c7b537a9LL) - (i >> 1);
        x = std::bit_cast<T>(i);
        x *= (T{ 1.5 } - xhalf * x * x);
        x *= (T{ 1.5 } - xhalf * x * x);
        return x;
    }

    //------------------------------------------------------------------------------
    // Min
    //------------------------------------------------------------------------------
    //
    // Returns the minimum of two values.
    //
    // Params:
    //  a - First value.
    //  b - Second value.
    //
    // Returns:
    //  Smaller of a and b.
    //
    // Notes:
    //  Comparable types; assumes valid inputs.
    //  No assert as trivial comparison.
    //
    template <typename T1, typename T2>
    constexpr auto Min(T1 a, T2 b) noexcept -> decltype(a + b)
    {
        return a < b ? a : b;
    }

    //------------------------------------------------------------------------------
    // Max
    //------------------------------------------------------------------------------
    //
    // Returns the maximum of two values.
    //
    // Params:
    //  a - First value.
    //  b - Second value.
    //
    // Returns:
    //  Larger of a and b.
    //
    // Notes:
    //  Comparable types; assumes valid inputs.
    //  No assert as trivial comparison.
    //
    template <typename T1, typename T2>
    constexpr auto Max(T1 a, T2 b) noexcept -> decltype(a + b)
    {
        return a > b ? a : b;
    }

    //------------------------------------------------------------------------------
    // isFinite
    //------------------------------------------------------------------------------
    //
    // Checks if a floating-point value is finite (not NaN or infinity).
    //
    // Params:
    //  x - Input value.
    //
    // Returns:
    //  true if x is finite; false otherwise.
    //
    // Notes:
    //  Uses std::isfinite from <cmath>.
    //  Assumes T is a floating-point type.
    //
    //------------------------------------------------------------------------------
    template <std::floating_point T>
    constexpr bool isFinite(T x) noexcept
    {
        return std::isfinite(x);
    }

    //------------------------------------------------------------------------------
    // CopySign
    //------------------------------------------------------------------------------
    //
    // Returns a value with the magnitude of x and the sign of y.
    //
    // Params:
    //  x - The value whose magnitude is used.
    //  y - The value whose sign is used.
    //
    // Returns:
    //  A value with the magnitude of x and the sign of y.
    //
    // Notes:
    //  Uses std::fabs from <cmath> to get the magnitude.
    //  Assumes T is a floating-point type.
    //
    //------------------------------------------------------------------------------
    template <std::floating_point T>
    constexpr T CopySign(T x, T y) noexcept
    {
        return xmath::Abs(x) * (y >= T(0) ? T(1) : T(-1));
    }

    //------------------------------------------------------------------------------
    // FEqual
    //------------------------------------------------------------------------------
    //
    // Checks if two floats are equal within tolerance.
    //
    // Params:
    //  f0 - First float.
    //  f1 - Second float.
    //  tol - Tolerance (default flt_tol_v).
    //
    // Returns:
    //  True if |f0 - f1| < tol.
    //
    // Notes:
    //  For floating-point comparison; tol >0.
    //  Assumes finite inputs; asserts validity.
    //
    template<std::floating_point T>
    inline bool FEqual(T f0, T f1, T tol) noexcept
    {
        assert(isValid(f0) && isValid(f1) && tol > T{ 0 });
        const T diff = f0 - f1;
        return (diff > -tol) && (diff < tol);
    }

    //------------------------------------------------------------------------------
    // FLess
    //------------------------------------------------------------------------------
    //
    // Checks if f0 is less than f1 within tolerance.
    //
    // Params:
    //  f0 - First float.
    //  f1 - Second float.
    //  tol - Tolerance (default flt_tol_v).
    //
    // Returns:
    //  True if (f0 - f1) < tol.
    //
    // Notes:
    //  Fuzzy < for floats; tol >0.
    //  Assumes finite inputs; asserts validity.
    //
    template<std::floating_point T>
    constexpr bool FLess(T f0, T f1, T tol) noexcept
    {
        assert(isValid(f0) && isValid(f1) && tol > T{ 0 });
        return (f0 - f1) < tol;
    }

    //------------------------------------------------------------------------------
    // FGreater
    //------------------------------------------------------------------------------
    //
    // Checks if f0 is greater than f1 within tolerance.
    //
    // Params:
    //  f0 - First float.
    //  f1 - Second float.
    //  tol - Tolerance (default flt_tol_v).
    //
    // Returns:
    //  True if (f0 - f1) > tol.
    //
    // Notes:
    //  Fuzzy > for floats; tol >0.
    //  Assumes finite inputs; asserts validity.
    //
    template<std::floating_point T>
    constexpr bool FGreater(T f0, T f1, T tol) noexcept
    {
        assert(isValid(f0) && isValid(f1) && tol > T{ 0 });
        return (f0 - f1) > tol;
    }

    //------------------------------------------------------------------------------
    // Sign
    //------------------------------------------------------------------------------
    //
    // Returns true if m_X >= 0 (positive or zero).
    //
    // Params:
    //  m_X - Input value.
    //
    // Returns:
    //  True for non-negative.
    //
    // Notes:
    //  Simple sign check; no zero special case.
    //  Assumes valid m_X; asserts if floating.
    //
    template <typename T>
    constexpr bool Sign(T x) noexcept
    {
        if constexpr (std::floating_point<T>) assert(isValid(x));
        return x >= T{ 0 };
    }

    //------------------------------------------------------------------------------
    // LRound
    //------------------------------------------------------------------------------
    //
    // Rounds float to nearest integer (towards positive infinity on tie).
    //
    // Params:
    //  m_X - Input float.
    //
    // Returns:
    //  Rounded int32_t.
    //
    // Notes:
    //  Uses std::round; handles negatives correctly.
    //  Assumes finite m_X in int32 range; asserts validity.
    //
    template<std::floating_point T>
    inline std::int32_t LRound(T x) noexcept
    {
        assert(isValid(x) && isInRange(x, T{ std::numeric_limits<T>::min() }, T{ std::numeric_limits<T>::max() }));
        return static_cast<std::int32_t>(std::round(x));
    }

    //------------------------------------------------------------------------------
    // Round
    //------------------------------------------------------------------------------
    //
    // Rounds a to nearest multiple of b.
    //
    // Params:
    //  a - Value to round.
    //  b - Rounding step (>0).
    //
    // Returns:
    //  Rounded value.
    //
    // Notes:
    //  Handles negatives; b !=0.
    //  Assumes finite/valid; asserts b >0 and validity.
    //
    template<std::floating_point T>
    inline T Round(T a, T b) noexcept
    {
        assert(isValid(a) && isValid(b) && b > T{ 0 });
        const T quotient = a / b;
        return (quotient < T{ 0 } ? std::ceil(quotient - T{ 0.5 }) : std::floor(quotient + T{ 0.5 })) * b;
    }

    //------------------------------------------------------------------------------
    // Round
    //------------------------------------------------------------------------------
    //
    // Rounds a floating-point value to the nearest integer.
    //
    // Params:
    //  m_X - Input value.
    //
    // Returns:
    //  Nearest integer to m_X, rounding halfway cases away from zero.
    //
    // Notes:
    //  Uses std::round; handles floating-point types.
    //  Asserts validity of input.
    //
    template <std::floating_point T>
    constexpr T Round(T x) noexcept
    {
        assert(isValid(x)); // Assumes isValid checks for finiteness, etc.
        return std::round(x);
    }


    //------------------------------------------------------------------------------
    // Ceil
    //------------------------------------------------------------------------------
    //
    // Computes smallest integer not less than m_X.
    //
    // Params:
    //  m_X - Input float.
    //
    // Returns:
    //  Ceiling value.
    //
    // Notes:
    //  Towards +inf; e.g., ceil(-1.1) = -1.
    //  Assumes finite m_X; asserts validity.
    //
    template<std::floating_point T>
    inline T Ceil(T x) noexcept
    {
        assert(isValid(x));
        return std::ceil(x);
    }

    //------------------------------------------------------------------------------
    // Floor
    //------------------------------------------------------------------------------
    //
    // Computes largest integer not greater than m_X.
    //
    // Params:
    //  m_X - Input float.
    //
    // Returns:
    //  Floor value.
    //
    // Notes:
    //  Towards -inf; e.g., floor(-1.1) = -2.
    //  Assumes finite m_X; asserts validity.
    //
    template<std::floating_point T>
    inline T Floor(T x) noexcept
    {
        assert(isValid(x));
        return std::floor(x);
    }

    //------------------------------------------------------------------------------
    // isInRange
    //------------------------------------------------------------------------------
    //
    // Checks if m_X is between min and max (inclusive).
    //
    // Params:
    //  m_X - Value to check.
    //  min - Lower bound.
    //  max - Upper bound (min <= max assumed).
    //
    // Returns:
    //  True if min <= m_X <= max.
    //
    // Notes:
    //  Comparable types; no floating-point tol.
    //  No assert as trivial; assumes min <= max.
    //
    template <typename T>
    constexpr bool isInRange(T x, T min, T max) noexcept
    {
        return (min <= x) && (x <= max);
    }

    //------------------------------------------------------------------------------
    // Range
    //------------------------------------------------------------------------------
    //
    // Clamps m_X to [min, max].
    //
    // Params:
    //  m_X - Value to clamp.
    //  min - Lower bound.
    //  max - Upper bound (min <= max assumed).
    //
    // Returns:
    //  Clamped value.
    //
    // Notes:
    //  Uses std::clamp in C++17+; min <= max.
    //  No assert; trivial clamp.
    //
    template <typename T>
    constexpr T Range(T x, T min, T max) noexcept
    {
        return (x < min) ? min : (x > max) ? max : x;
    }

    //------------------------------------------------------------------------------
    // Abs
    //------------------------------------------------------------------------------
    //
    // Computes absolute value of m_X.
    //
    // Params:
    //  m_X - Input value.
    //
    // Returns:
    //  |m_X| (non-negative).
    //
    // Notes:
    //  Uses std::abs; handles signed types.
    //  Assumes finite for floats; asserts validity if floating.
    //
    template <typename T>
    constexpr T Abs(T x) noexcept
    {
        if constexpr (std::floating_point<T>) assert(isValid(x));
        return std::abs(x);
    }

    //------------------------------------------------------------------------------
    // Clamp
    //------------------------------------------------------------------------------
    //
    // Clamps a value between a minimum and maximum bound.
    //
    // Params:
    //  m_X     - Input value.
    //  m_Low   - Minimum bound.
    //  m_High  - Maximum bound.
    //
    // Returns:
    //  m_X clamped to [m_Low, m_High].
    //
    // Notes:
    //  Assumes m_Low <= m_High.
    //  Works with types supporting comparison operators.
    //
    template <typename T>
    constexpr T Clamp(const T& value, const T& low, const T& high) noexcept
    {
        assert(!(high < low)); // Ensure valid bounds
        return (value < low) ? low : (value > high) ? high : value;
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Linear interpolation between a and b.
    //
    // Params:
    //  t - Factor [0,1].
    //  a - Start value.
    //  b - End value.
    //
    // Returns:
    //  a + t*(b - a).
    //
    // Notes:
    //  Uses std::lerp; t clamped implicitly.
    //  Assumes finite a/b/t, t in [0,1]; asserts validity.
    //
    template <typename T>
    constexpr T Lerp(float t, T a, T b) noexcept
    {
        assert(isValid(t) && isValid(a) && isValid(b) && isInRange(t, 0.0f, 1.0f));
        return a + t * (b - a);
    }

    //------------------------------------------------------------------------------
    // Trunc
    //------------------------------------------------------------------------------
    //
    // Truncates the fractional part of a floating-point number.
    //
    // Params:
    //  x - Input value.
    //
    // Returns:
    //  The integer part of x, as a floating-point value.
    //
    // Notes:
    //  Behaves like std::trunc.
    //  Assumes finite x; asserts validity.
    //
    //------------------------------------------------------------------------------
    template <std::floating_point T>
    constexpr T Trunc(T x) noexcept
    {
        assert(isValid(x));
        return x < T(0) ? std::ceil(x) : std::floor(x);
    }

    //------------------------------------------------------------------------------
    // isValid
    //------------------------------------------------------------------------------
    //
    // Checks if floating-point m_X is finite (not NaN/inf).
    //
    // Params:
    //  m_X - Input float.
    //
    // Returns:
    //  True if finite.
    //
    // Notes:
    //  Uses std::isfinite; essential for safety in math ops.
    //  No assert as check function.
    //
    template<std::floating_point T>
    inline bool isValid(T x) noexcept
    {
        return std::isfinite(x);
    }

    //------------------------------------------------------------------------------
    // SolvedQuadraticRoots
    //------------------------------------------------------------------------------
    //
    // Solves quadratic equation ax^2 + bx + c = 0 for real roots.
    //
    // Params:
    //  root1 - First root output.
    //  root2 - Second root output.
    //  a - Coefficient a (non-zero for quadratic).
    //  b - Coefficient b.
    //  c - Coefficient c.
    //
    // Returns:
    //  True if real roots found (discriminant >=0).
    //
    // Notes:
    //  Handles degenerate cases (linear/const); stable computation to avoid cancellation.
    //  Assumes finite coeffs; asserts validity; roots unordered.
    //
    inline bool SolvedQuadraticRoots(float& root1, float& root2, float a, float b, float c) noexcept
    {
        assert(isValid(a) && isValid(b) && isValid(c));
        if (a == 0.0f)
        {
            if (b == 0.0f) return false;
            root1 = root2 = -c / b;
            return true;
        }
        if (b == 0.0f)
        {
            const float neg_c_over_a = -c / a;
            if (neg_c_over_a < 0.0f) return false;
            root1 = Sqrt(neg_c_over_a);
            root2 = -root1;
            return true;
        }
        const float part1 = b / (2.0f * a);
        const float part2 = c / (a * part1 * part1);
        const float one_minus_part2 = 1.0f - part2;
        if (one_minus_part2 < 0.0f) return false;
        const float part3 = -1.0f - Sqrt(one_minus_part2);
        root1 = part1 * part3;
        root2 = (part1 * part2) / part3;
        return true;
    }
}