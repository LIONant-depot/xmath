#ifndef XMATH_IMATH_H
    #error "You should not include this header directly, just need to include xmath_imath.h"
#endif
#pragma once

namespace xmath
{
    //------------------------------------------------------------------------------
    // ivec2
    //------------------------------------------------------------------------------
    //
    // 2D vector class for integers without SIMD optimization.
    //
    // Notes:
    //  This class does not initialize its memory if default constructed. 
    //  Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    //  Immutable/copy ops suffixed _copy.
    //  Use const for safety; asserts in inline for validity (components within int range, but always finite).
    //  Targets C++20: constexpr where possible, no unicode, specific spacing.
    //
    struct ivec2
    {
        union
        {
            std::array<int, 2>      m_Elements;
            struct
            {
                int m_X, m_Y;
            };
        };

        // Constructors
        constexpr                       ivec2                           (void)                                              noexcept = default;
        constexpr                       ivec2                           (ivec2&&)                                           noexcept = default;
        constexpr                       ivec2                           (const ivec2&)                                      noexcept = default;
        constexpr                       ivec2                           (int x, int y)                                      noexcept;
        constexpr                       ivec2                           (int value)                                         noexcept;
        constexpr                       ivec2                           (std::span<int> Span)                               noexcept;
        constexpr                       ivec2                           (const std::array<double, 2>& Conversion)           noexcept;
        constexpr                       ivec2                           (const std::array<float, 2>& Conversion)            noexcept;

        // Assignment and conversion operators
        constexpr           ivec2&      operator =                      (const ivec2&)                                      noexcept = default;
        constexpr                       operator std::array<int, 2>     (void)                                  const       noexcept;
        constexpr                       operator std::array<float, 2>   (void)                                  const       noexcept;
        inline                          operator std::string(void)                                              const       noexcept;
        inline              std::string ToString                        (void)                                  const       noexcept;
        inline friend std::ostream&     operator<<                      (std::ostream& os, const ivec2& vec)                noexcept;

        // Static properties
        static constexpr    ivec2       fromZero                        (void)                                              noexcept;
        static constexpr    ivec2       fromOne                         (void)                                              noexcept;
        static constexpr    ivec2       fromUnit                        (void)                                              noexcept;
        static constexpr    ivec2       fromUp                          (void)                                              noexcept;
        static constexpr    ivec2       fromDown                        (void)                                              noexcept;
        static constexpr    ivec2       fromLeft                        (void)                                              noexcept;
        static constexpr    ivec2       fromRight                       (void)                                              noexcept;

        // Setup functions
        constexpr           ivec2&      setupMax                        (void)                                              noexcept;
        constexpr           ivec2&      setupZero                       (void)                                              noexcept;
        constexpr           ivec2&      setup                           ( int x, int y )                                    noexcept;


        // Static methods
        static constexpr    long long   Dot                             (const ivec2& a, const ivec2& b)                    noexcept;
        static constexpr    ivec2       Min                             (const ivec2& a, const ivec2& b)                    noexcept;
        static constexpr    ivec2       Max                             (const ivec2& a, const ivec2& b)                    noexcept;
        static inline       float       Distance                        (const ivec2& a, const ivec2& b)                    noexcept;
        static constexpr    long long   Cross                           (const ivec2& a, const ivec2& b)                    noexcept;

        // Static methods as members
        constexpr           long long   Dot                             (const ivec2& a)                        const       noexcept;
        constexpr           ivec2       Min                             (const ivec2& a)                        const       noexcept;
        constexpr           ivec2       Max                             (const ivec2& a)                        const       noexcept;
        inline              float       Distance                        (const ivec2& a)                        const       noexcept;

        // Instance methods - Basic operations
        inline              float       Length                          (void)                                  const       noexcept;
        constexpr           long long   LengthSquared                   (void)                                  const       noexcept;
        constexpr           int         LengthManhattan                 (void)                                  const       noexcept;
        constexpr           bool        isInRange                       (int min, int max)                      const       noexcept;
        constexpr           bool        Equals                          (const ivec2& other)                    const       noexcept;

        // Instance methods - Component-wise math
        constexpr           ivec2       AbsCopy                         (void)                                  const       noexcept;
        constexpr           ivec2&      Abs                             (void)                                              noexcept;
        constexpr           ivec2       SignCopy                        (void)                                  const       noexcept;
        constexpr           ivec2&      Sign                            (void)                                              noexcept;
        inline              ivec2       ModCopy                         (int divisor)                           const       noexcept;
        inline              ivec2&      Mod                             (int divisor)                                       noexcept;
        inline              ivec2       ClampCopy                       (int min_val, int max_val)              const       noexcept;
        inline              ivec2&      Clamp                           (int min_val, int max_val)                          noexcept;
        inline              ivec2       ClampCopy                       (const ivec2& min, const ivec2& max)    const       noexcept;
        inline              ivec2&      Clamp                           (const ivec2& min, const ivec2& max)                noexcept;

        // Instance methods - Geometry
        constexpr           long long   DistanceSquare                  (const ivec2& v)                        const       noexcept;
        constexpr           ivec2&      GridSnap                        (int gridX, int gridY)                              noexcept;
        constexpr           ivec2       Perp                            (void)                                  const       noexcept;
        constexpr           long long   WhichSideOfLine                 (const ivec2& V0, const ivec2& V1)      const       noexcept;

        // Swizzle methods for int (HLSL-style, return copy with swizzled components)
        constexpr           int         x                               (void)                                  const       noexcept;
        constexpr           int         y                               (void)                                  const       noexcept;

        // Swizzle methods for ivec2 (HLSL-style, return copy with swizzled components)
        constexpr           ivec2       xx                              (void)                                  const       noexcept;
        constexpr           ivec2       xy                              (void)                                  const       noexcept;
        constexpr           ivec2       yx                              (void)                                  const       noexcept;
        constexpr           ivec2       yy                              (void)                                  const       noexcept;

        // Operator overloads
        constexpr           ivec2       operator+                       (const ivec2& other)                    const       noexcept;
        constexpr           ivec2       operator-                       (const ivec2& other)                    const       noexcept;
        constexpr           ivec2       operator*                       (int scalar)                            const       noexcept;
        constexpr           ivec2       operator/                       (int scalar)                            const       noexcept;
        constexpr           ivec2&      operator+=                      (const ivec2& other)                                noexcept;
        constexpr           ivec2&      operator-=                      (const ivec2& other)                                noexcept;
        constexpr           ivec2&      operator*=                      (int scalar)                                        noexcept;
        constexpr           ivec2&      operator/=                      (int scalar)                                        noexcept;
        constexpr           bool        operator==                      (const ivec2& other)                    const       noexcept;
        constexpr           bool        operator!=                      (const ivec2& other)                    const       noexcept;
        constexpr           int         operator[]                      (std::int32_t index)                    const       noexcept;
        constexpr           int&        operator[]                      (std::int32_t index)                                noexcept;

        // Friend operators
        friend constexpr    ivec2       operator*                       (int scalar, const ivec2& v)                        noexcept;
        friend constexpr    ivec2       operator-                       (const ivec2& v)                                    noexcept;
    };
}