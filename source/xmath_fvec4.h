#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
    #error "You must include xmath_vector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // fvec4
    //------------------------------------------------------------------------------
    //
    // 4D vector class with SIMD optimization (SSE).
    //
    // Notes:
    //  This class does not initialize its memory if default constructed. 
    //  Aligned to 16 bytes for SIMD.
    //  Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    //  Immutable/copy ops suffixed _copy.
    //  Use const for safety; asserts in inline for validity (isfinite components).
    //  Targets C++20: constexpr where possible, no unicode, specific spacing.
    //
    struct alignas(16) fvec4
    {
        union
        {
            floatx4                 m_XYZW;
            std::array<float, 4>    m_Elements;
            struct
            {
                float m_X, m_Y, m_Z, m_W;
            };
        };

        // Constructors
        constexpr                       fvec4                       (void)                                              noexcept = default;
        constexpr                       fvec4                       (fvec4&&)                                           noexcept = default;
        constexpr                       fvec4                       (const fvec4&)                                      noexcept = default;
        constexpr                       fvec4&                      operator = (const fvec4&)                           noexcept = default;
        inline                          fvec4                       (float x, float y, float z, float w)                noexcept;
        inline                          fvec4                       (float value)                                       noexcept;
        inline                          fvec4                       (const fvec3& other, float w)                       noexcept;
        inline                          fvec4                       (float x, const fvec3& other)                       noexcept;
        constexpr explicit              fvec4                       (const floatx4& reg)                                noexcept;
        inline                          fvec4                       (const fvec2& xy, const fvec2& zw)                  noexcept;
        inline                          fvec4                       (std::span<float> Span)                             noexcept;

        // Assignment and conversion operators
        constexpr                       operator std::array<double,4>(void)                                 const       noexcept;
        inline                          operator std::string        (void)                                  const       noexcept;
        std::string                     ToString                    (void)                                  const       noexcept;
        inline friend std::ostream&     operator<<                  (std::ostream& os, const fvec4& vec)                noexcept;

        // Static properties
        static constexpr    fvec4       fromZero                    (void)                                              noexcept;
        static constexpr    fvec4       fromOne                     (void)                                              noexcept;
        static constexpr    fvec4       fromUnitX                   (void)                                              noexcept;
        static constexpr    fvec4       fromUnitY                   (void)                                              noexcept;
        static constexpr    fvec4       fromUnitZ                   (void)                                              noexcept;
        static constexpr    fvec4       fromUnitW                   (void)                                              noexcept;
        static constexpr    fvec4       fromUp                      (void)                                              noexcept;
        static constexpr    fvec4       fromDown                    (void)                                              noexcept;
        static constexpr    fvec4       fromLeft                    (void)                                              noexcept;
        static constexpr    fvec4       fromRight                   (void)                                              noexcept;
        static constexpr    fvec4       fromForward                 (void)                                              noexcept;
        static constexpr    fvec4       fromBack                    (void)                                              noexcept;
        static inline       fvec4       fromRandomUnitVector        (void)                                              noexcept;

        // Static methods
        static inline       float       Dot                         (const fvec4& a, const fvec4& b)                    noexcept;
        static inline       fvec4       Min                         (const fvec4& a, const fvec4& b)                    noexcept;
        static inline       fvec4       Max                         (const fvec4& a, const fvec4& b)                    noexcept;
        static inline       fvec4       Lerp                        (const fvec4& a, const fvec4& b, float t)           noexcept;
        static inline       float       Distance                    (const fvec4& a, const fvec4& b)                    noexcept;

        // Static methods as members
        inline              float       Dot                         (const fvec4& a)                            const   noexcept;
        inline              fvec4       Min                         (const fvec4& a)                            const   noexcept;
        inline              fvec4       Max                         (const fvec4& a)                            const   noexcept;
        inline              fvec4       Lerp                        (const fvec4& a, float t)                   const   noexcept;
        inline              float       Distance                    (const fvec4& a)                            const   noexcept;

        // Instance methods
        inline              float       Length                      (void)                                      const   noexcept;
        inline              float       LengthSq                    (void)                                      const   noexcept;
        inline              fvec4       LimitLengthCopy             (float MaxLength)                           const   noexcept;
        inline              fvec4       NormalizeCopy               (void)                                      const   noexcept;
        inline              fvec4&      Normalize                   (void)                                              noexcept;
        inline              fvec4       HomogeneousCopy             (void)                                      const   noexcept;
        inline              fvec4&      Homogenize                  (void)                                              noexcept;
        inline              fvec4       NormalizeSafeCopy           (void)                                      const   noexcept;
        inline              fvec4&      NormalizeSafe               (void)                                              noexcept;
        inline              bool        isFinite                    (void)                                      const   noexcept;
        inline              bool        isInRange                   (float min, float max)                      const   noexcept;
        inline              fvec4       Reflection                  (const fvec4& normal)                       const   noexcept;
        inline              float       DistanceSquare              (const fvec4& v)                            const   noexcept;
        inline              radian      AngleBetween                (const fvec4& v)                            const   noexcept;
        inline              fvec4&      GridSnap                    (float gridX, float gridY, float gridZ, float gridW)noexcept;

        // Instance methods - Component-wise math
        inline              fvec4       AbsCopy                     (void)                                      const   noexcept;
        inline              fvec4&      Abs                         (void)                                              noexcept;
        inline              fvec4       OneOverCopy                 (void)                                      const   noexcept;
        inline              fvec4&      OneOver                     (void)                                              noexcept;
        inline              fvec4       SqrtCopy                    ( void )                                    const   noexcept;
        inline              fvec4&      Sqrt                        ( void )                                            noexcept;
        inline              fvec4       InvSqrtFastCopy             ( void )                                    const   noexcept;
        inline              fvec4       InvSqrtCopy                 ( void )                                    const   noexcept;
        inline              fvec4&      InvSqrt                     ( void )                                            noexcept;
        inline              fvec4&      InvSqrtFast                 ( void )                                            noexcept;
        inline              fvec4       SignCopy                    ( void )                                    const   noexcept;
        inline              fvec4&      Sign                        ( void )                                            noexcept;
        inline              fvec4       FloorCopy                   ( void )                                    const   noexcept;
        inline              fvec4&      Floor                       ( void )                                            noexcept;
        inline              fvec4       CeilCopy                    ( void )                                    const   noexcept;
        inline              fvec4&      Ceil                        ( void )                                            noexcept;
        inline              fvec4       FractCopy                   ( void )                                    const   noexcept;
        inline              fvec4&      Fract                       ( void )                                            noexcept;
        inline              fvec4       RoundCopy                   ( void )                                    const   noexcept;
        inline              fvec4&      Round                       ( void )                                            noexcept;
        inline              fvec4       TruncCopy                   ( void )                                    const   noexcept;
        inline              fvec4&      Trunc                       ( void )                                            noexcept;
        inline              fvec4       ModCopy                     ( float divisor )                           const   noexcept;
        inline              fvec4&      Mod                         ( float divisor )                                   noexcept;
        inline              fvec4       ClampCopy                   ( float min_val, float max_val)             const   noexcept;
        inline              fvec4&      Clamp                       ( float min_val, float max_val)                     noexcept;
        inline              fvec4       ClampCopy                   ( const fvec4& min, const fvec4& max)       const   noexcept;
        inline              fvec4&      Clamp                       ( const fvec4& min, const fvec4& max)               noexcept;
        inline              fvec4       Step                        ( float edge)                               const   noexcept;
        inline              fvec4       SmoothStep                  ( float edge0, float edge1 )                const   noexcept;
        inline              fvec4       LogCopy                     ( void )                                    const   noexcept;
        inline              fvec4&      Log                         ( void )                                            noexcept;
        inline              fvec4       Log2Copy                    ( void )                                    const   noexcept;
        inline              fvec4&      Log2                        ( void )                                            noexcept;
        inline              fvec4       PowCopy                     (float exp)                                 const   noexcept;
        inline              fvec4&      Pow                         (float exp)                                         noexcept;
        inline              fvec4       SinCopy                     ( void )                                    const   noexcept;
        inline              fvec4&      Sin                         ( void )                                            noexcept;
        inline              fvec4       CosCopy                     ( void )                                    const   noexcept;
        inline              fvec4&      Cos                         ( void )                                            noexcept;
        inline              fvec4       TanCopy                     ( void )                                    const   noexcept;
        inline              fvec4&      Tan                         ( void )                                            noexcept;
        inline              fvec4       AsinCopy                    ( void )                                    const   noexcept;
        inline              fvec4&      Asin                        ( void )                                            noexcept;
        inline              fvec4       AcosCopy                    ( void )                                    const   noexcept;
        inline              fvec4&      Acos                        ( void )                                            noexcept;
        inline              fvec4       AtanCopy                    ( void )                                    const   noexcept;
        inline              fvec4&      Atan                        ( void )                                            noexcept;
        inline              fvec4       Atan2Copy                   (const fvec4& x)                            const   noexcept;
        inline              fvec4&      Atan2                       (const fvec4& x)                                    noexcept;

        // Swizzle methods for float (HLSL-style, return copy with swizzled components)
        inline              float       x                           (void)                                      const   noexcept;
        inline              float       y                           (void)                                      const   noexcept;
        inline              float       z                           (void)                                      const   noexcept;
        inline              float       w                           (void)                                      const   noexcept;

        // Swizzle methods for fvec2 (HLSL-style, return copy with swizzled components)
        inline              fvec2       xx                          (void)                                      const   noexcept;
        inline              fvec2       xy                          (void)                                      const   noexcept;
        inline              fvec2       xz                          (void)                                      const   noexcept;
        inline              fvec2       xw                          (void)                                      const   noexcept;
        inline              fvec2       yx                          (void)                                      const   noexcept;
        inline              fvec2       yy                          (void)                                      const   noexcept;
        inline              fvec2       yz                          (void)                                      const   noexcept;
        inline              fvec2       yw                          (void)                                      const   noexcept;
        inline              fvec2       zx                          (void)                                      const   noexcept;
        inline              fvec2       zy                          (void)                                      const   noexcept;
        inline              fvec2       zz                          (void)                                      const   noexcept;
        inline              fvec2       zw                          (void)                                      const   noexcept;
        inline              fvec2       wx                          (void)                                      const   noexcept;
        inline              fvec2       wy                          (void)                                      const   noexcept;
        inline              fvec2       wz                          (void)                                      const   noexcept;
        inline              fvec2       ww                          (void)                                      const   noexcept;

        // Swizzle methods for fvec3 (HLSL-style, return copy with swizzled components)
        inline              fvec3       xxx                         (void)                                      const   noexcept;
        inline              fvec3       xxy                         (void)                                      const   noexcept;
        inline              fvec3       xxz                         (void)                                      const   noexcept;
        inline              fvec3       xxw                         (void)                                      const   noexcept;
        inline              fvec3       xyx                         (void)                                      const   noexcept;
        inline              fvec3       xyy                         (void)                                      const   noexcept;
        inline              fvec3       xyz                         (void)                                      const   noexcept;
        inline              fvec3       xyw                         (void)                                      const   noexcept;
        inline              fvec3       xzx                         (void)                                      const   noexcept;
        inline              fvec3       xzy                         (void)                                      const   noexcept;
        inline              fvec3       xzz                         (void)                                      const   noexcept;
        inline              fvec3       xzw                         (void)                                      const   noexcept;
        inline              fvec3       xwx                         (void)                                      const   noexcept;
        inline              fvec3       xwy                         (void)                                      const   noexcept;
        inline              fvec3       xwz                         (void)                                      const   noexcept;
        inline              fvec3       xww                         (void)                                      const   noexcept;
        inline              fvec3       yxx                         (void)                                      const   noexcept;
        inline              fvec3       yxy                         (void)                                      const   noexcept;
        inline              fvec3       yxz                         (void)                                      const   noexcept;
        inline              fvec3       yxw                         (void)                                      const   noexcept;
        inline              fvec3       yyx                         (void)                                      const   noexcept;
        inline              fvec3       yyy                         (void)                                      const   noexcept;
        inline              fvec3       yyz                         (void)                                      const   noexcept;
        inline              fvec3       yyw                         (void)                                      const   noexcept;
        inline              fvec3       yzx                         (void)                                      const   noexcept;
        inline              fvec3       yzy                         (void)                                      const   noexcept;
        inline              fvec3       yzz                         (void)                                      const   noexcept;
        inline              fvec3       yzw                         (void)                                      const   noexcept;
        inline              fvec3       ywx                         (void)                                      const   noexcept;
        inline              fvec3       ywy                         (void)                                      const   noexcept;
        inline              fvec3       ywz                         (void)                                      const   noexcept;
        inline              fvec3       yww                         (void)                                      const   noexcept;
        inline              fvec3       zxx                         (void)                                      const   noexcept;
        inline              fvec3       zxy                         (void)                                      const   noexcept;
        inline              fvec3       zxz                         (void)                                      const   noexcept;
        inline              fvec3       zxw                         (void)                                      const   noexcept;
        inline              fvec3       zyx                         (void)                                      const   noexcept;
        inline              fvec3       zyy                         (void)                                      const   noexcept;
        inline              fvec3       zyz                         (void)                                      const   noexcept;
        inline              fvec3       zyw                         (void)                                      const   noexcept;
        inline              fvec3       zzx                         (void)                                      const   noexcept;
        inline              fvec3       zzy                         (void)                                      const   noexcept;
        inline              fvec3       zzz                         (void)                                      const   noexcept;
        inline              fvec3       zzw                         (void)                                      const   noexcept;
        inline              fvec3       zwx                         (void)                                      const   noexcept;
        inline              fvec3       zwy                         (void)                                      const   noexcept;
        inline              fvec3       zwz                         (void)                                      const   noexcept;
        inline              fvec3       zww                         (void)                                      const   noexcept;
        inline              fvec3       wxx                         (void)                                      const   noexcept;
        inline              fvec3       wxy                         (void)                                      const   noexcept;
        inline              fvec3       wxz                         (void)                                      const   noexcept;
        inline              fvec3       wxw                         (void)                                      const   noexcept;
        inline              fvec3       wyx                         (void)                                      const   noexcept;
        inline              fvec3       wyy                         (void)                                      const   noexcept;
        inline              fvec3       wyz                         (void)                                      const   noexcept;
        inline              fvec3       wyw                         (void)                                      const   noexcept;
        inline              fvec3       wzx                         (void)                                      const   noexcept;
        inline              fvec3       wzy                         (void)                                      const   noexcept;
        inline              fvec3       wzz                         (void)                                      const   noexcept;
        inline              fvec3       wzw                         (void)                                      const   noexcept;
        inline              fvec3       wwx                         (void)                                      const   noexcept;
        inline              fvec3       wwy                         (void)                                      const   noexcept;
        inline              fvec3       wwz                         (void)                                      const   noexcept;
        inline              fvec3       www                         (void)                                      const   noexcept;

        // Swizzle methods for fvec4 (HLSL-style, return copy with swizzled components)
        inline              fvec4       xxxx                        (void)                                      const   noexcept;
        inline              fvec4       xxxy                        (void)                                      const   noexcept;
        inline              fvec4       xxxz                        (void)                                      const   noexcept;
        inline              fvec4       xxxw                        (void)                                      const   noexcept;
        inline              fvec4       xxyx                        (void)                                      const   noexcept;
        inline              fvec4       xxyy                        (void)                                      const   noexcept;
        inline              fvec4       xxyz                        (void)                                      const   noexcept;
        inline              fvec4       xxyw                        (void)                                      const   noexcept;
        inline              fvec4       xxzx                        (void)                                      const   noexcept;
        inline              fvec4       xxzy                        (void)                                      const   noexcept;
        inline              fvec4       xxzz                        (void)                                      const   noexcept;
        inline              fvec4       xxzw                        (void)                                      const   noexcept;
        inline              fvec4       xxwx                        (void)                                      const   noexcept;
        inline              fvec4       xxwy                        (void)                                      const   noexcept;
        inline              fvec4       xxwz                        (void)                                      const   noexcept;
        inline              fvec4       xxww                        (void)                                      const   noexcept;
        inline              fvec4       xyxx                        (void)                                      const   noexcept;
        inline              fvec4       xyxy                        (void)                                      const   noexcept;
        inline              fvec4       xyxz                        (void)                                      const   noexcept;
        inline              fvec4       xyxw                        (void)                                      const   noexcept;
        inline              fvec4       xyyx                        (void)                                      const   noexcept;
        inline              fvec4       xyyy                        (void)                                      const   noexcept;
        inline              fvec4       xyyz                        (void)                                      const   noexcept;
        inline              fvec4       xyyw                        (void)                                      const   noexcept;
        inline              fvec4       xyzx                        (void)                                      const   noexcept;
        inline              fvec4       xyzy                        (void)                                      const   noexcept;
        inline              fvec4       xyzz                        (void)                                      const   noexcept;
        inline              fvec4       xyzw                        (void)                                      const   noexcept;
        inline              fvec4       xywx                        (void)                                      const   noexcept;
        inline              fvec4       xywy                        (void)                                      const   noexcept;
        inline              fvec4       xywz                        (void)                                      const   noexcept;
        inline              fvec4       xyww                        (void)                                      const   noexcept;
        inline              fvec4       xzxx                        (void)                                      const   noexcept;
        inline              fvec4       xzxy                        (void)                                      const   noexcept;
        inline              fvec4       xzxz                        (void)                                      const   noexcept;
        inline              fvec4       xzxw                        (void)                                      const   noexcept;
        inline              fvec4       xzyx                        (void)                                      const   noexcept;
        inline              fvec4       xzyy                        (void)                                      const   noexcept;
        inline              fvec4       xzyz                        (void)                                      const   noexcept;
        inline              fvec4       xzyw                        (void)                                      const   noexcept;
        inline              fvec4       xzzx                        (void)                                      const   noexcept;
        inline              fvec4       xzzy                        (void)                                      const   noexcept;
        inline              fvec4       xzzz                        (void)                                      const   noexcept;
        inline              fvec4       xzzw                        (void)                                      const   noexcept;
        inline              fvec4       xzwx                        (void)                                      const   noexcept;
        inline              fvec4       xzwy                        (void)                                      const   noexcept;
        inline              fvec4       xzwz                        (void)                                      const   noexcept;
        inline              fvec4       xzww                        (void)                                      const   noexcept;
        inline              fvec4       xwxx                        (void)                                      const   noexcept;
        inline              fvec4       xwxy                        (void)                                      const   noexcept;
        inline              fvec4       xwxz                        (void)                                      const   noexcept;
        inline              fvec4       xwxw                        (void)                                      const   noexcept;
        inline              fvec4       xwyx                        (void)                                      const   noexcept;
        inline              fvec4       xwyy                        (void)                                      const   noexcept;
        inline              fvec4       xwyz                        (void)                                      const   noexcept;
        inline              fvec4       xwyw                        (void)                                      const   noexcept;
        inline              fvec4       xwzx                        (void)                                      const   noexcept;
        inline              fvec4       xwzy                        (void)                                      const   noexcept;
        inline              fvec4       xwzz                        (void)                                      const   noexcept;
        inline              fvec4       xwzw                        (void)                                      const   noexcept;
        inline              fvec4       xwwx                        (void)                                      const   noexcept;
        inline              fvec4       xwwy                        (void)                                      const   noexcept;
        inline              fvec4       xwwz                        (void)                                      const   noexcept;
        inline              fvec4       xwww                        (void)                                      const   noexcept;
        inline              fvec4       yxxx                        (void)                                      const   noexcept;
        inline              fvec4       yxxy                        (void)                                      const   noexcept;
        inline              fvec4       yxxz                        (void)                                      const   noexcept;
        inline              fvec4       yxxw                        (void)                                      const   noexcept;
        inline              fvec4       yxyx                        (void)                                      const   noexcept;
        inline              fvec4       yxyy                        (void)                                      const   noexcept;
        inline              fvec4       yxyz                        (void)                                      const   noexcept;
        inline              fvec4       yxyw                        (void)                                      const   noexcept;
        inline              fvec4       yxzx                        (void)                                      const   noexcept;
        inline              fvec4       yxzy                        (void)                                      const   noexcept;
        inline              fvec4       yxzz                        (void)                                      const   noexcept;
        inline              fvec4       yxzw                        (void)                                      const   noexcept;
        inline              fvec4       yxwx                        (void)                                      const   noexcept;
        inline              fvec4       yxwy                        (void)                                      const   noexcept;
        inline              fvec4       yxwz                        (void)                                      const   noexcept;
        inline              fvec4       yxww                        (void)                                      const   noexcept;
        inline              fvec4       yyxx                        (void)                                      const   noexcept;
        inline              fvec4       yyxy                        (void)                                      const   noexcept;
        inline              fvec4       yyxz                        (void)                                      const   noexcept;
        inline              fvec4       yyxw                        (void)                                      const   noexcept;
        inline              fvec4       yyyx                        (void)                                      const   noexcept;
        inline              fvec4       yyyy                        (void)                                      const   noexcept;
        inline              fvec4       yyyz                        (void)                                      const   noexcept;
        inline              fvec4       yyyw                        (void)                                      const   noexcept;
        inline              fvec4       yyzx                        (void)                                      const   noexcept;
        inline              fvec4       yyzy                        (void)                                      const   noexcept;
        inline              fvec4       yyzz                        (void)                                      const   noexcept;
        inline              fvec4       yyzw                        (void)                                      const   noexcept;
        inline              fvec4       yywx                        (void)                                      const   noexcept;
        inline              fvec4       yywy                        (void)                                      const   noexcept;
        inline              fvec4       yywz                        (void)                                      const   noexcept;
        inline              fvec4       yyww                        (void)                                      const   noexcept;
        inline              fvec4       yzxx                        (void)                                      const   noexcept;
        inline              fvec4       yzxy                        (void)                                      const   noexcept;
        inline              fvec4       yzxz                        (void)                                      const   noexcept;
        inline              fvec4       yzxw                        (void)                                      const   noexcept;
        inline              fvec4       yzyx                        (void)                                      const   noexcept;
        inline              fvec4       yzyy                        (void)                                      const   noexcept;
        inline              fvec4       yzyz                        (void)                                      const   noexcept;
        inline              fvec4       yzyw                        (void)                                      const   noexcept;
        inline              fvec4       yzzx                        (void)                                      const   noexcept;
        inline              fvec4       yzzy                        (void)                                      const   noexcept;
        inline              fvec4       yzzz                        (void)                                      const   noexcept;
        inline              fvec4       yzzw                        (void)                                      const   noexcept;
        inline              fvec4       yzwx                        (void)                                      const   noexcept;
        inline              fvec4       yzwy                        (void)                                      const   noexcept;
        inline              fvec4       yzwz                        (void)                                      const   noexcept;
        inline              fvec4       yzww                        (void)                                      const   noexcept;
        inline              fvec4       ywxx                        (void)                                      const   noexcept;
        inline              fvec4       ywxy                        (void)                                      const   noexcept;
        inline              fvec4       ywxz                        (void)                                      const   noexcept;
        inline              fvec4       ywxw                        (void)                                      const   noexcept;
        inline              fvec4       ywyx                        (void)                                      const   noexcept;
        inline              fvec4       ywyy                        (void)                                      const   noexcept;
        inline              fvec4       ywyz                        (void)                                      const   noexcept;
        inline              fvec4       ywyw                        (void)                                      const   noexcept;
        inline              fvec4       ywzx                        (void)                                      const   noexcept;
        inline              fvec4       ywzy                        (void)                                      const   noexcept;
        inline              fvec4       ywzz                        (void)                                      const   noexcept;
        inline              fvec4       ywzw                        (void)                                      const   noexcept;
        inline              fvec4       ywwx                        (void)                                      const   noexcept;
        inline              fvec4       ywwy                        (void)                                      const   noexcept;
        inline              fvec4       ywwz                        (void)                                      const   noexcept;
        inline              fvec4       ywww                        (void)                                      const   noexcept;
        inline              fvec4       zxxx                        (void)                                      const   noexcept;
        inline              fvec4       zxxy                        (void)                                      const   noexcept;
        inline              fvec4       zxxz                        (void)                                      const   noexcept;
        inline              fvec4       zxxw                        (void)                                      const   noexcept;
        inline              fvec4       zxyx                        (void)                                      const   noexcept;
        inline              fvec4       zxyy                        (void)                                      const   noexcept;
        inline              fvec4       zxyz                        (void)                                      const   noexcept;
        inline              fvec4       zxyw                        (void)                                      const   noexcept;
        inline              fvec4       zxzx                        (void)                                      const   noexcept;
        inline              fvec4       zxzy                        (void)                                      const   noexcept;
        inline              fvec4       zxzz                        (void)                                      const   noexcept;
        inline              fvec4       zxzw                        (void)                                      const   noexcept;
        inline              fvec4       zxwx                        (void)                                      const   noexcept;
        inline              fvec4       zxwy                        (void)                                      const   noexcept;
        inline              fvec4       zxwz                        (void)                                      const   noexcept;
        inline              fvec4       zxww                        (void)                                      const   noexcept;
        inline              fvec4       zyxx                        (void)                                      const   noexcept;
        inline              fvec4       zyxy                        (void)                                      const   noexcept;
        inline              fvec4       zyxz                        (void)                                      const   noexcept;
        inline              fvec4       zyxw                        (void)                                      const   noexcept;
        inline              fvec4       zyyx                        (void)                                      const   noexcept;
        inline              fvec4       zyyy                        (void)                                      const   noexcept;
        inline              fvec4       zyyz                        (void)                                      const   noexcept;
        inline              fvec4       zyyw                        (void)                                      const   noexcept;
        inline              fvec4       zyzx                        (void)                                      const   noexcept;
        inline              fvec4       zyzy                        (void)                                      const   noexcept;
        inline              fvec4       zyzz                        (void)                                      const   noexcept;
        inline              fvec4       zyzw                        (void)                                      const   noexcept;
        inline              fvec4       zywx                        (void)                                      const   noexcept;
        inline              fvec4       zywy                        (void)                                      const   noexcept;
        inline              fvec4       zywz                        (void)                                      const   noexcept;
        inline              fvec4       zyww                        (void)                                      const   noexcept;
        inline              fvec4       zzxx                        (void)                                      const   noexcept;
        inline              fvec4       zzxy                        (void)                                      const   noexcept;
        inline              fvec4       zzxz                        (void)                                      const   noexcept;
        inline              fvec4       zzxw                        (void)                                      const   noexcept;
        inline              fvec4       zzyx                        (void)                                      const   noexcept;
        inline              fvec4       zzyy                        (void)                                      const   noexcept;
        inline              fvec4       zzyz                        (void)                                      const   noexcept;
        inline              fvec4       zzyw                        (void)                                      const   noexcept;
        inline              fvec4       zzzx                        (void)                                      const   noexcept;
        inline              fvec4       zzzy                        (void)                                      const   noexcept;
        inline              fvec4       zzzz                        (void)                                      const   noexcept;
        inline              fvec4       zzzw                        (void)                                      const   noexcept;
        inline              fvec4       zzwx                        (void)                                      const   noexcept;
        inline              fvec4       zzwy                        (void)                                      const   noexcept;
        inline              fvec4       zzwz                        (void)                                      const   noexcept;
        inline              fvec4       zzww                        (void)                                      const   noexcept;
        inline              fvec4       zwxx                        (void)                                      const   noexcept;
        inline              fvec4       zwxy                        (void)                                      const   noexcept;
        inline              fvec4       zwxz                        (void)                                      const   noexcept;
        inline              fvec4       zwxw                        (void)                                      const   noexcept;
        inline              fvec4       zwyx                        (void)                                      const   noexcept;
        inline              fvec4       zwyy                        (void)                                      const   noexcept;
        inline              fvec4       zwyz                        (void)                                      const   noexcept;
        inline              fvec4       zwyw                        (void)                                      const   noexcept;
        inline              fvec4       zwzx                        (void)                                      const   noexcept;
        inline              fvec4       zwzy                        (void)                                      const   noexcept;
        inline              fvec4       zwzz                        (void)                                      const   noexcept;
        inline              fvec4       zwzw                        (void)                                      const   noexcept;
        inline              fvec4       zwwx                        (void)                                      const   noexcept;
        inline              fvec4       zwwy                        (void)                                      const   noexcept;
        inline              fvec4       zwwz                        (void)                                      const   noexcept;
        inline              fvec4       zwww                        (void)                                      const   noexcept;
        inline              fvec4       wxxx                        (void)                                      const   noexcept;
        inline              fvec4       wxxy                        (void)                                      const   noexcept;
        inline              fvec4       wxxz                        (void)                                      const   noexcept;
        inline              fvec4       wxxw                        (void)                                      const   noexcept;
        inline              fvec4       wxyx                        (void)                                      const   noexcept;
        inline              fvec4       wxyy                        (void)                                      const   noexcept;
        inline              fvec4       wxyz                        (void)                                      const   noexcept;
        inline              fvec4       wxyw                        (void)                                      const   noexcept;
        inline              fvec4       wxzx                        (void)                                      const   noexcept;
        inline              fvec4       wxzy                        (void)                                      const   noexcept;
        inline              fvec4       wxzz                        (void)                                      const   noexcept;
        inline              fvec4       wxzw                        (void)                                      const   noexcept;
        inline              fvec4       wxwx                        (void)                                      const   noexcept;
        inline              fvec4       wxwy                        (void)                                      const   noexcept;
        inline              fvec4       wxwz                        (void)                                      const   noexcept;
        inline              fvec4       wxww                        (void)                                      const   noexcept;
        inline              fvec4       wyxx                        (void)                                      const   noexcept;
        inline              fvec4       wyxy                        (void)                                      const   noexcept;
        inline              fvec4       wyxz                        (void)                                      const   noexcept;
        inline              fvec4       wyxw                        (void)                                      const   noexcept;
        inline              fvec4       wyyx                        (void)                                      const   noexcept;
        inline              fvec4       wyyy                        (void)                                      const   noexcept;
        inline              fvec4       wyyz                        (void)                                      const   noexcept;
        inline              fvec4       wyyw                        (void)                                      const   noexcept;
        inline              fvec4       wyzx                        (void)                                      const   noexcept;
        inline              fvec4       wyzy                        (void)                                      const   noexcept;
        inline              fvec4       wyzz                        (void)                                      const   noexcept;
        inline              fvec4       wyzw                        (void)                                      const   noexcept;
        inline              fvec4       wywx                        (void)                                      const   noexcept;
        inline              fvec4       wywy                        (void)                                      const   noexcept;
        inline              fvec4       wywz                        (void)                                      const   noexcept;
        inline              fvec4       wyww                        (void)                                      const   noexcept;
        inline              fvec4       wzxx                        (void)                                      const   noexcept;
        inline              fvec4       wzxy                        (void)                                      const   noexcept;
        inline              fvec4       wzxz                        (void)                                      const   noexcept;
        inline              fvec4       wzxw                        (void)                                      const   noexcept;
        inline              fvec4       wzyx                        (void)                                      const   noexcept;
        inline              fvec4       wzyy                        (void)                                      const   noexcept;
        inline              fvec4       wzyz                        (void)                                      const   noexcept;
        inline              fvec4       wzyw                        (void)                                      const   noexcept;
        inline              fvec4       wzzx                        (void)                                      const   noexcept;
        inline              fvec4       wzzy                        (void)                                      const   noexcept;
        inline              fvec4       wzzz                        (void)                                      const   noexcept;
        inline              fvec4       wzzw                        (void)                                      const   noexcept;
        inline              fvec4       wzwx                        (void)                                      const   noexcept;
        inline              fvec4       wzwy                        (void)                                      const   noexcept;
        inline              fvec4       wzwz                        (void)                                      const   noexcept;
        inline              fvec4       wzww                        (void)                                      const   noexcept;
        inline              fvec4       wwxx                        (void)                                      const   noexcept;
        inline              fvec4       wwxy                        (void)                                      const   noexcept;
        inline              fvec4       wwxz                        (void)                                      const   noexcept;
        inline              fvec4       wwxw                        (void)                                      const   noexcept;
        inline              fvec4       wwyx                        (void)                                      const   noexcept;
        inline              fvec4       wwyy                        (void)                                      const   noexcept;
        inline              fvec4       wwyz                        (void)                                      const   noexcept;
        inline              fvec4       wwyw                        (void)                                      const   noexcept;
        inline              fvec4       wwzx                        (void)                                      const   noexcept;
        inline              fvec4       wwzy                        (void)                                      const   noexcept;
        inline              fvec4       wwzz                        (void)                                      const   noexcept;
        inline              fvec4       wwzw                        (void)                                      const   noexcept;
        inline              fvec4       wwwx                        (void)                                      const   noexcept;
        inline              fvec4       wwwy                        (void)                                      const   noexcept;
        inline              fvec4       wwwz                        (void)                                      const   noexcept;
        inline              fvec4       wwww                        (void)                                      const   noexcept;

        // Operator overloads
        inline              fvec4       operator+                   (const fvec4& other)                        const   noexcept;
        inline              fvec4       operator-                   (const fvec4& other)                        const   noexcept;
        inline              fvec4       operator*                   (float scalar)                              const   noexcept;
        inline              fvec4       operator/                   (float scalar)                              const   noexcept;
        inline              fvec4&      operator+=                  (const fvec4& other)                                noexcept;
        inline              fvec4&      operator-=                  (const fvec4& other)                                noexcept;
        inline              fvec4&      operator*=                  (float scalar)                                      noexcept;
        inline              fvec4&      operator/=                  (float scalar)                                      noexcept;
        inline              bool        operator==                  (const fvec4& other)                        const   noexcept;
        inline              bool        operator!=                  (const fvec4& other)                        const   noexcept;
        inline              float       operator[]                  (std::int32_t index)                        const   noexcept;
        inline              float&      operator[]                  (std::int32_t index)                                noexcept;

        // Friend operators
        friend inline       fvec4       operator*                   (float scalar, const fvec4& v)                      noexcept;
        friend inline       fvec4       operator-                   (const fvec4& v)                                    noexcept;
    };
}