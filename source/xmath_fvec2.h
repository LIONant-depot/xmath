#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
    #error "You must include xmath_fvector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // fvec2
    //------------------------------------------------------------------------------
    //
    // 2D vector class without SIMD optimization.
    //
    // Notes:
    //  This class does not initialize its memory if default constructed. 
    //  Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    //  Immutable/copy ops suffixed _copy.
    //  Use const for safety; asserts in inline for validity (isfinite components).
    //  Targets C++20: constexpr where possible, no unicode, specific spacing.
    //
    struct fvec2
    {
        union
        {
            std::array<float, 2>    m_Elements;
            struct
            {
                float m_X, m_Y;
            };
        };

        // Constructors
        constexpr                       fvec2                   (void)                                              noexcept = default;
        constexpr                       fvec2                   (fvec2&&)                                           noexcept = default;
        constexpr                       fvec2                   (const fvec2&)                                      noexcept = default;
        constexpr                       fvec2                   (float x, float y)                                  noexcept;
        constexpr                       fvec2                   (float value)                                       noexcept;
        constexpr                       fvec2                   (std::span<float> Span)                             noexcept;
        constexpr                       fvec2                   (const std::array<double,2>& Conversion )           noexcept;

        // Assignment and conversion operators
        constexpr           fvec2&      operator =              (const fvec2&)                                      noexcept = default;
        constexpr                       operator std::array<double,2> (void)                            const       noexcept;
        inline                          operator std::string    (void)                                  const       noexcept;
        std::string                     ToString                (void)                                  const       noexcept;
        inline friend std::ostream&     operator<<              (std::ostream& os, const fvec2& vec)                noexcept;

        // Static properties
        static constexpr    fvec2       fromZero                (void)                                              noexcept;
        static constexpr    fvec2       fromOne                 (void)                                              noexcept;
        static constexpr    fvec2       fromUnitX               (void)                                              noexcept;
        static constexpr    fvec2       fromUnitY               (void)                                              noexcept;
        static constexpr    fvec2       fromUp                  (void)                                              noexcept;
        static constexpr    fvec2       fromDown                (void)                                              noexcept;
        static constexpr    fvec2       fromLeft                (void)                                              noexcept;
        static constexpr    fvec2       fromRight               (void)                                              noexcept;
        static inline       fvec2       fromRandomUnitVector    (void)                                              noexcept;

        // Static methods
        static constexpr    float       Dot                     (const fvec2& a, const fvec2& b)                    noexcept;
        static constexpr    fvec2       Min                     (const fvec2& a, const fvec2& b)                    noexcept;
        static constexpr    fvec2       Max                     (const fvec2& a, const fvec2& b)                    noexcept;
        static constexpr    fvec2       Lerp                    (const fvec2& a, const fvec2& b, float t)           noexcept;
        static inline       float       Distance                (const fvec2& a, const fvec2& b)                    noexcept;
        static constexpr    float       Cross                   (const fvec2& a, const fvec2& b)                    noexcept;

        // Static methods as members
        constexpr           float       Dot                     (const fvec2& a)                        const       noexcept;
        constexpr           fvec2       Min                     (const fvec2& a)                        const       noexcept;
        constexpr           fvec2       Max                     (const fvec2& a)                        const       noexcept;
        constexpr           fvec2       Lerp                    (const fvec2& a, float t)               const       noexcept;
        inline              float       Distance                (const fvec2& a)                        const       noexcept;
        constexpr           float       Cross                   (const fvec2& a)                        const       noexcept;

        // Instance methods - Basic operations
        inline              float       Length                  (void)                                  const       noexcept;
        constexpr           float       LengthSq                (void)                                  const       noexcept;
        inline              fvec2       NormalizeCopy           (void)                                  const       noexcept;
        inline              fvec2&      Normalize               (void)                                              noexcept;
        inline              fvec2       NormalizeSafeCopy       (void)                                  const       noexcept;
        inline              fvec2&      NormalizeSafe           (void)                                              noexcept;
        inline              fvec2       LimitLengthCopy         (float MaxLength)                       const       noexcept;
        inline              bool        isFinite                (void)                                  const       noexcept;
        constexpr           bool        isInRange               (float min, float max)                  const       noexcept;
        inline              bool        Equals                  (const fvec2& other, float tolerance)   const       noexcept;

        // Instance methods - Component-wise math
        constexpr           fvec2       AbsCopy                 (void)                                  const       noexcept;
        constexpr           fvec2&      Abs                     (void)                                              noexcept;
        constexpr           fvec2       OneOverCopy             (void)                                  const       noexcept;
        constexpr           fvec2&      OneOver                 (void)                                              noexcept;
        inline              fvec2       SqrtCopy                ( void )                                const       noexcept;
        inline              fvec2&      Sqrt                    ( void )                                            noexcept;
        inline              fvec2       InvSqrtCopy             ( void )                                const       noexcept;
        inline              fvec2&      InvSqrt                 ( void )                                            noexcept;
        inline              fvec2       SignCopy                ( void )                                const       noexcept;
        inline              fvec2&      Sign                    ( void )                                            noexcept;
        inline              fvec2       FloorCopy               ( void )                                const       noexcept;
        inline              fvec2&      Floor                   ( void )                                            noexcept;
        inline              fvec2       CeilCopy                ( void )                                const       noexcept;
        inline              fvec2&      Ceil                    ( void )                                            noexcept;
        inline              fvec2       FractCopy               ( void )                                const       noexcept;
        inline              fvec2&      Fract                   ( void )                                            noexcept;
        inline              fvec2       RoundCopy               ( void )                                const       noexcept;
        inline              fvec2&      Round                   ( void )                                            noexcept;
        inline              fvec2       TruncCopy               ( void )                                const       noexcept;
        inline              fvec2&      Trunc                   ( void )                                            noexcept;
        inline              fvec2       ModCopy                 ( float divisor )                       const       noexcept;
        inline              fvec2&      Mod                     ( float divisor )                                   noexcept;
        inline              fvec2       ClampCopy               ( float min_val, float max_val)         const       noexcept;
        inline              fvec2&      Clamp                   ( float min_val, float max_val)                     noexcept;
        inline              fvec2       ClampCopy               ( const fvec2& min, const fvec2& max )  const       noexcept;
        inline              fvec2&      Clamp                   ( const fvec2& min, const fvec2& max )              noexcept;
        inline              fvec2       Step                    ( float edge)                           const       noexcept;
        inline              fvec2       SmoothStep              ( float edge0, float edge1 )            const       noexcept;
        inline              fvec2       LogCopy                 ( void )                                const       noexcept;
        inline              fvec2&      Log                     ( void )                                            noexcept;
        inline              fvec2       Log2Copy                ( void )                                const       noexcept;
        inline              fvec2&      Log2                    ( void )                                            noexcept;
        inline              fvec2       PowCopy                 (float exp)                             const       noexcept;
        inline              fvec2&      Pow                     (float exp)                                         noexcept;
        inline              fvec2       SinCopy                 ( void )                                const       noexcept;
        inline              fvec2&      Sin                     ( void )                                            noexcept;
        inline              fvec2       CosCopy                 ( void )                                const       noexcept;
        inline              fvec2&      Cos                     ( void )                                            noexcept;
        inline              fvec2       TanCopy                 ( void )                                const       noexcept;
        inline              fvec2&      Tan                     ( void )                                            noexcept;
        inline              fvec2       AsinCopy                ( void )                                const       noexcept;
        inline              fvec2&      Asin                    ( void )                                            noexcept;
        inline              fvec2       AcosCopy                ( void )                                const       noexcept;
        inline              fvec2&      Acos                    ( void )                                            noexcept;
        inline              fvec2       AtanCopy                ( void )                                const       noexcept;
        inline              fvec2&      Atan                    ( void )                                            noexcept;
        inline              fvec2       Atan2Copy               (const fvec2& x)                        const       noexcept;
        inline              fvec2&      Atan2                   ( const fvec2& x )                                  noexcept;
        inline              radian      SignedAngleBetween      (const fvec2& v)                        const       noexcept;

        // Component-wise exponential functions
        inline              fvec2       ExpCopy                 ( void )                                const       noexcept;
        inline              fvec2&      Exp                     ( void )                                            noexcept;

        // Instance methods - Geometry
        constexpr           fvec2       Reflection              (const fvec2& normal)                   const       noexcept;
        constexpr           float       DistanceSquare          (const fvec2& v)                        const       noexcept;
        inline              radian      AngleBetween            (const fvec2& v)                        const       noexcept;
        constexpr           fvec2&      GridSnap                (float gridX, float gridY)                          noexcept;
        constexpr           fvec2       Perp                    (void)                                  const       noexcept;
        constexpr           float       WhichSideOfLine         ( const fvec2& V0, const fvec2& V1 )    const       noexcept;
        inline              fvec2       ClosestPointInLine      ( const fvec2& V0, const fvec2& V1 )    const       noexcept;
        inline              fvec2       ClosestPointInLineSegment( const fvec2& V0, const fvec2& V1 )   const       noexcept;
        inline              fvec2       RotateCopy              (radian angle)                          const       noexcept;
        inline              fvec2&      Rotate                  (radian angle)                                      noexcept;
        inline              fvec2       ProjectCopy             (const fvec2& onto)                     const       noexcept;
        inline              fvec2&      Project                 (const fvec2& onto)                                 noexcept;

        // Swizzle methods for float (HLSL-style, return copy with swizzled components)
        constexpr           float       x                       (void)                                  const       noexcept;
        constexpr           float       y                       (void)                                  const       noexcept;

        // Swizzle methods for fvec2 (HLSL-style, return copy with swizzled components)
        constexpr           fvec2       xx                      (void)                                  const       noexcept;
        constexpr           fvec2       xy                      (void)                                  const       noexcept;
        constexpr           fvec2       yx                      (void)                                  const       noexcept;
        constexpr           fvec2       yy                      (void)                                  const       noexcept;

        // Swizzle methods for fvec3 (HLSL-style, return copy with swizzled components)
        constexpr           fvec3       xxx                     (void)                                  const       noexcept;
        constexpr           fvec3       xxy                     (void)                                  const       noexcept;
        constexpr           fvec3       xyx                     (void)                                  const       noexcept;
        constexpr           fvec3       xyy                     (void)                                  const       noexcept;
        constexpr           fvec3       yxx                     (void)                                  const       noexcept;
        constexpr           fvec3       yxy                     (void)                                  const       noexcept;
        constexpr           fvec3       yyx                     (void)                                  const       noexcept;
        constexpr           fvec3       yyy                     (void)                                  const       noexcept;

        // Swizzle methods for fvec4 (HLSL-style, return copy with swizzled components)
        inline              fvec4       xxxx                    (void)                                  const       noexcept;
        inline              fvec4       xxxy                    (void)                                  const       noexcept;
        inline              fvec4       xxyx                    (void)                                  const       noexcept;
        inline              fvec4       xxyy                    (void)                                  const       noexcept;
        inline              fvec4       xyxx                    (void)                                  const       noexcept;
        inline              fvec4       xyxy                    (void)                                  const       noexcept;
        inline              fvec4       xyyx                    (void)                                  const       noexcept;
        inline              fvec4       xyyy                    (void)                                  const       noexcept;
        inline              fvec4       yxxx                    (void)                                  const       noexcept;
        inline              fvec4       yxxy                    (void)                                  const       noexcept;
        inline              fvec4       yxyx                    (void)                                  const       noexcept;
        inline              fvec4       yxyy                    (void)                                  const       noexcept;
        inline              fvec4       yyxx                    (void)                                  const       noexcept;
        inline              fvec4       yyxy                    (void)                                  const       noexcept;
        inline              fvec4       yyyx                    (void)                                  const       noexcept;
        inline              fvec4       yyyy                    (void)                                  const       noexcept;

        // Operator overloads
        constexpr           fvec2       operator+               (const fvec2& other)                    const       noexcept;
        constexpr           fvec2       operator-               (const fvec2& other)                    const       noexcept;
        constexpr           fvec2       operator*               (float scalar)                          const       noexcept;
        constexpr           fvec2       operator/               (float scalar)                          const       noexcept;
        constexpr           fvec2&      operator+=              (const fvec2& other)                                noexcept;
        constexpr           fvec2&      operator-=              (const fvec2& other)                                noexcept;
        constexpr           fvec2&      operator*=              (float scalar)                                      noexcept;
        constexpr           fvec2&      operator/=              (float scalar)                                      noexcept;
        constexpr           bool        operator==              (const fvec2& other)                    const       noexcept;
        constexpr           bool        operator!=              (const fvec2& other)                    const       noexcept;
        constexpr           float       operator[]              (std::int32_t index)                    const       noexcept;
        constexpr           float&      operator[]              (std::int32_t index)                                noexcept;

        // Friend operators
        friend constexpr    fvec2       operator*               (float scalar, const fvec2& v)                      noexcept;
        friend constexpr    fvec2       operator-               (const fvec2& v)                                    noexcept;
    };
}
