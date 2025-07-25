#pragma once
#ifndef XMATH_FLOAT_VECTOR_H
    #error "You must include xmath_fvector.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // supported data structure for fvec3_t:
    //------------------------------------------------------------------------------
    namespace details::f_vec3
    {
        struct alignas(16) simd_data
        {
            union
            {
                floatx4             m_XYZW;
                std::array<float,3> m_Elements;
                struct
                {
                    float m_X, m_Y, m_Z, m_W;
                };
            };
        };

        struct cpu_data
        {
            union
            {
                std::array<float, 3> m_Elements;
                struct
                {
                    float m_X, m_Y, m_Z;
                };
            };
        };
    }

    //------------------------------------------------------------------------------
    // fvec3_t
    //------------------------------------------------------------------------------
    //
    // 3D vector class with SIMD optimization (SSE).
    //
    // Notes:
    //  This class does not initialize its memory if default constructed. 
    //  Aligned to 16 bytes for SIMD; pads with m_W, assume trash for W.
    //  Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    //  Immutable/copy ops suffixed _copy.
    //  Use const for safety; asserts in inline for validity (isfinite components).
    //  Targets C++20: constexpr\consteval where possible.
    //
    template <bool T_USE_SIMD_V >
    struct fvec3_t : std::conditional_t< T_USE_SIMD_V, details::f_vec3::simd_data, details::f_vec3::cpu_data >
    {
        // Constructors
        constexpr                       fvec3_t                (void)                                               noexcept = default;
        constexpr                       fvec3_t                (float x, float y, float z)                          noexcept;
        constexpr                       fvec3_t                (float value)                                        noexcept;
        inline                          fvec3_t                (radian pitch, radian yaw)                           noexcept;
        constexpr explicit              fvec3_t                (const floatx4& reg)                                 noexcept requires T_USE_SIMD_V;
        constexpr                       fvec3_t                (const fvec3_t<!T_USE_SIMD_V>& other)                noexcept;
        inline                          fvec3_t                (const fvec2& other, float z )                       noexcept;
        inline                          fvec3_t                (float x, const fvec2& other )                       noexcept;
        inline                          fvec3_t                (std::span<float> Span )                             noexcept;
        constexpr                       fvec3_t                (const std::array<double,3>& Conversion )            noexcept;

        // Assignment and conversion operators
        constexpr                       operator std::array<double,3> ()                                const       noexcept;
        inline                          operator std::string    ()                                      const       noexcept;
        std::string                     ToString                (void)                                  const       noexcept;
        template<bool V>
        inline friend std::ostream&     operator<<              (std::ostream& os, const fvec3_t<V>& vec)           noexcept;

        // Static properties
        static consteval    fvec3_t     fromZero                (void)                                              noexcept;
        static consteval    fvec3_t     fromOne                 (void)                                              noexcept;
        static consteval    fvec3_t     fromUp                  (void)                                              noexcept;
        static consteval    fvec3_t     fromDown                (void)                                              noexcept;
        static consteval    fvec3_t     fromLeft                (void)                                              noexcept;
        static consteval    fvec3_t     fromRight               (void)                                              noexcept;
        static consteval    fvec3_t     fromForward             (void)                                              noexcept;
        static consteval    fvec3_t     fromBack                (void)                                              noexcept;
        static inline       fvec3_t     fromRandomUnitVector    (void)                                              noexcept;

        // Static methods
        static inline       float       Dot                     (const fvec3_t& a, const fvec3_t& b)                noexcept;
        static inline       fvec3_t     Cross                   (const fvec3_t& a, const fvec3_t& b)                noexcept;
        static inline       fvec3_t     Min                     (const fvec3_t& a, const fvec3_t& b)                noexcept;
        static inline       fvec3_t     Max                     (const fvec3_t& a, const fvec3_t& b)                noexcept;
        static inline       fvec3_t     Lerp                    (const fvec3_t& a, const fvec3_t& b, float t)       noexcept;
        static inline       float       Distance                (const fvec3_t& a, const fvec3_t& b)                noexcept;

        // Static methods as members
        inline              float       Dot                     (const fvec3_t& a)                      const       noexcept;
        inline              fvec3_t     Cross                   (const fvec3_t& a)                      const       noexcept;
        inline              fvec3_t     Min                     (const fvec3_t& a)                      const       noexcept;
        inline              fvec3_t     Max                     (const fvec3_t& a)                      const       noexcept;
        inline              fvec3_t     Lerp                    (const fvec3_t& a, float t)             const       noexcept;
        inline              float       Distance                (const fvec3_t& a)                      const       noexcept;
        inline              float       DistanceSquare          (const fvec3_t& v)                      const       noexcept;

        // Instance methods - Basic operations
        inline              float       Length                  (void)                                  const       noexcept;
        inline              float       LengthSq                (void)                                  const       noexcept;
        inline              fvec3_t     NormalizeCopy           (void)                                  const       noexcept;
        inline              fvec3_t&    Normalize               (void)                                              noexcept;
        inline              fvec3_t     NormalizeSafeCopy       (void)                                  const       noexcept;
        inline              fvec3_t&    NormalizeSafe           (void)                                              noexcept;
        inline              fvec3_t     LimitLengthCopy         (float MaxLength)                       const       noexcept;
        inline              bool        isFinite                (void)                                  const       noexcept;
        inline              bool        isInRange               (float min, float max)                  const       noexcept;
        inline              bool        Equals                  (const fvec3_t& other, float tolerance) const       noexcept;

        // Instance methods - Component-wise math
        inline              fvec3_t     AbsCopy                 (void)                                  const       noexcept;
        inline              fvec3_t&    Abs                     (void)                                              noexcept;
        inline              fvec3_t     OneOverCopy             (void)                                  const       noexcept;
        inline              fvec3_t&    OneOver                 (void)                                              noexcept;
        inline              fvec3_t     SqrtCopy                ( void )                                const       noexcept;
        inline              fvec3_t&    Sqrt                    ( void )                                            noexcept;
        inline              fvec3_t     InvSqrtCopy             ( void )                                const       noexcept;
        inline              fvec3_t&    InvSqrt                 ( void )                                            noexcept;
        inline              fvec3_t     SignCopy                ( void )                                const       noexcept;
        inline              fvec3_t&    Sign                    ( void )                                            noexcept;
        inline              fvec3_t     FloorCopy               ( void )                                const       noexcept;
        inline              fvec3_t&    Floor                   ( void )                                            noexcept;
        inline              fvec3_t     CeilCopy                ( void )                                const       noexcept;
        inline              fvec3_t&    Ceil                    ( void )                                            noexcept;
        inline              fvec3_t     FractCopy               ( void )                                const       noexcept;
        inline              fvec3_t&    Fract                   ( void )                                            noexcept;
        inline              fvec3_t     RoundCopy               ( void )                                const       noexcept;
        inline              fvec3_t&    Round                   ( void )                                            noexcept;
        inline              fvec3_t     TruncCopy               ( void )                                const       noexcept;
        inline              fvec3_t&    Trunc                   ( void )                                            noexcept;
        inline              fvec3_t     ModCopy                 ( float divisor )                       const       noexcept;
        inline              fvec3_t&    Mod                     ( float divisor )                                   noexcept;
        inline              fvec3_t     ClampCopy               ( float min_val, float max_val)         const       noexcept;
        inline              fvec3_t&    Clamp                   ( float min_val, float max_val)                     noexcept;
        inline              fvec3_t     ClampCopy               ( const fvec3_t& min, const fvec3_t& max) const     noexcept;
        inline              fvec3_t&    Clamp                   ( const fvec3_t& min, const fvec3_t& max)           noexcept;
        inline              fvec3_t     Step                    ( float edge)                           const       noexcept;
        inline              fvec3_t     SmoothStep              ( float edge0, float edge1 )            const       noexcept;
        inline              fvec3_t     LogCopy                 ( void )                                const       noexcept;
        inline              fvec3_t&    Log                     ( void )                                            noexcept;
        inline              fvec3_t     Log2Copy                ( void )                                const       noexcept;
        inline              fvec3_t&    Log2                    ( void )                                            noexcept;
        inline              fvec3_t     PowCopy                 (float exp)                             const       noexcept;
        inline              fvec3_t&    Pow                     (float exp)                                         noexcept;
        inline              fvec3_t&    Exp                     ( void )                                            noexcept;
        inline              fvec3_t     ExpCopy                 (void )                                 const       noexcept;
        inline              fvec3_t     SinCopy                 ( void )                                const       noexcept;
        inline              fvec3_t&    Sin                     ( void )                                            noexcept;
        inline              fvec3_t     CosCopy                 ( void )                                const       noexcept;
        inline              fvec3_t&    Cos                     ( void )                                            noexcept;
        inline              fvec3_t     TanCopy                 ( void )                                const       noexcept;
        inline              fvec3_t&    Tan                     ( void )                                            noexcept;
        inline              fvec3_t     AsinCopy                ( void )                                const       noexcept;
        inline              fvec3_t&    Asin                    ( void )                                            noexcept;
        inline              fvec3_t     AcosCopy                ( void )                                const       noexcept;
        inline              fvec3_t&    Acos                    ( void )                                            noexcept;
        inline              fvec3_t     AtanCopy                ( void )                                const       noexcept;
        inline              fvec3_t&    Atan                    ( void )                                            noexcept;
        inline              fvec3_t     Atan2Copy               (const fvec3_t& x)                      const       noexcept;
        inline              fvec3_t&    Atan2                   ( const fvec3_t& x )                                noexcept;

        // Instance methods - Geometry
        inline              radian      Pitch                   (void)                                  const       noexcept;
        inline              radian      Yaw                     (void)                                  const       noexcept;
        inline std::pair<radian, radian> PitchYaw               (void)                                  const       noexcept;
        inline std::pair<fvec3_t, radian> RotationTowards      (const fvec3_t& dest)                    const       noexcept;
        inline              radian      AngleBetween            (const fvec3_t& v)                      const       noexcept;
        inline              radian      SignedAngleBetween      (const fvec3_t& v)                      const       noexcept;
        inline              fvec3_t     VectorToLineSegment     (const fvec3_t& start, const fvec3_t& end) const    noexcept;
        inline              float       SquareDistToLineSeg     (const fvec3_t& start, const fvec3_t& end) const    noexcept;
        inline              fvec3_t     ClosestPointInLineSegment(const fvec3_t& start, const fvec3_t& end) const   noexcept;
        inline              float       ClosestPointToRectangle (const fvec3_t& p0, const fvec3_t& e0, const fvec3_t& e1, fvec3_t& outClosestPoint) const noexcept;
        inline              fvec3_t     RotateXCopy             (radian rx)                             const       noexcept;
        inline              fvec3_t     RotateYCopy             (radian ry)                             const       noexcept;
        inline              fvec3_t     RotateZCopy             (radian rz)                             const       noexcept;
        inline              fvec3_t&    RotateX                 (radian rx)                                         noexcept;
        inline              fvec3_t&    RotateY                 (radian ry)                                         noexcept;
        inline              fvec3_t&    RotateZ                 (radian rz)                                         noexcept;
        inline              fvec3_t     RotateCopy              (const radian3& r)                      const       noexcept;
        inline              fvec3_t&    Rotate                  (const radian3& r)                                  noexcept;
        inline              fvec3_t&    RotateInverse           (const radian3& r)                                  noexcept;
        inline              fvec3_t     RotateInverseCopy       (const radian3& r)                      const       noexcept;
        inline              fvec3_t     Reflection              (const fvec3_t& normal)                 const       noexcept;
        inline              fvec3_t&    GridSnap                (float gridX, float gridY, float gridZ)             noexcept;
        inline              bool        isRightHanded           (const fvec3_t& p1, const fvec3_t& p2)  const       noexcept;
        inline              fvec3_t     ProjectCopy             (const fvec3_t& onto)                   const       noexcept;
        inline              fvec3_t&    Project                 (const fvec3_t& onto)                               noexcept;
        inline              fvec3_t     Perpendicular           (const fvec3_t& normal)                 const       noexcept;
        inline              fvec3_t     ProjectOntoPlane        (const fvec3_t& normal)                 const       noexcept;
        inline              bool        isNearlyZero            (float tolerance = 1e-6f)               const       noexcept;
        inline              bool        isNormalized            (float tolerance = 1e-6f)               const       noexcept;
        inline              fvec3_t     MoveTowardsCopy         (const fvec3_t& target, float maxDistanceDelta) const noexcept;
        inline              fvec3_t&    MoveTowards             (const fvec3_t& target, float maxDistanceDelta)     noexcept;

        // Swizzle methods for float (HLSL-style, return copy with swizzled components)
        inline              float       x                       (void)                                  const       noexcept;
        inline              float       y                       (void)                                  const       noexcept;
        inline              float       z                       (void)                                  const       noexcept;

        // Swizzle methods for fvec2 (HLSL-style, return copy with swizzled components)
        inline              fvec2       xx                      (void)                                  const       noexcept;
        inline              fvec2       xy                      (void)                                  const       noexcept;
        inline              fvec2       xz                      (void)                                  const       noexcept;
        inline              fvec2       yx                      (void)                                  const       noexcept;
        inline              fvec2       yy                      (void)                                  const       noexcept;
        inline              fvec2       yz                      (void)                                  const       noexcept;
        inline              fvec2       zx                      (void)                                  const       noexcept;
        inline              fvec2       zy                      (void)                                  const       noexcept;
        inline              fvec2       zz                      (void)                                  const       noexcept;

        // Swizzle methods for fvec3 (HLSL-style, return copy with swizzled components)
        inline              fvec3_t     xxx                     ( void )                                const       noexcept;
        inline              fvec3_t     xxy                     ( void )                                const       noexcept;
        inline              fvec3_t     xxz                     ( void )                                const       noexcept;
        inline              fvec3_t     xyx                     ( void )                                const       noexcept;
        inline              fvec3_t     xyy                     ( void )                                const       noexcept;
        inline              fvec3_t     xyz                     ( void )                                const       noexcept;
        inline              fvec3_t     xzx                     ( void )                                const       noexcept;
        inline              fvec3_t     xzy                     ( void )                                const       noexcept;
        inline              fvec3_t     xzz                     ( void )                                const       noexcept;
        inline              fvec3_t     yxx                     ( void )                                const       noexcept;
        inline              fvec3_t     yxy                     ( void )                                const       noexcept;
        inline              fvec3_t     yxz                     ( void )                                const       noexcept;
        inline              fvec3_t     yyx                     ( void )                                const       noexcept;
        inline              fvec3_t     yyy                     ( void )                                const       noexcept;
        inline              fvec3_t     yyz                     ( void )                                const       noexcept;
        inline              fvec3_t     yzx                     ( void )                                const       noexcept;
        inline              fvec3_t     yzy                     ( void )                                const       noexcept;
        inline              fvec3_t     yzz                     ( void )                                const       noexcept;
        inline              fvec3_t     zxx                     ( void )                                const       noexcept;
        inline              fvec3_t     zxy                     ( void )                                const       noexcept;
        inline              fvec3_t     zxz                     ( void )                                const       noexcept;
        inline              fvec3_t     zyx                     ( void )                                const       noexcept;
        inline              fvec3_t     zyy                     ( void )                                const       noexcept;
        inline              fvec3_t     zyz                     ( void )                                const       noexcept;
        inline              fvec3_t     zzx                     ( void )                                const       noexcept;
        inline              fvec3_t     zzy                     ( void )                                const       noexcept;
        inline              fvec3_t     zzz                     ( void )                                const       noexcept;

        // Swizzle methods for fvec4 (HLSL-style, return copy with swizzled components)
        inline              fvec4       xxxx                    (void)                                  const       noexcept;
        inline              fvec4       xxxy                    (void)                                  const       noexcept;
        inline              fvec4       xxxz                    (void)                                  const       noexcept;
        inline              fvec4       xxyx                    (void)                                  const       noexcept;
        inline              fvec4       xxyy                    (void)                                  const       noexcept;
        inline              fvec4       xxyz                    (void)                                  const       noexcept;
        inline              fvec4       xxzx                    (void)                                  const       noexcept;
        inline              fvec4       xxzy                    (void)                                  const       noexcept;
        inline              fvec4       xxzz                    (void)                                  const       noexcept;
        inline              fvec4       xyxx                    (void)                                  const       noexcept;
        inline              fvec4       xyxy                    (void)                                  const       noexcept;
        inline              fvec4       xyxz                    (void)                                  const       noexcept;
        inline              fvec4       xyyx                    (void)                                  const       noexcept;
        inline              fvec4       xyyy                    (void)                                  const       noexcept;
        inline              fvec4       xyyz                    (void)                                  const       noexcept;
        inline              fvec4       xyzx                    (void)                                  const       noexcept;
        inline              fvec4       xyzy                    (void)                                  const       noexcept;
        inline              fvec4       xyzz                    (void)                                  const       noexcept;
        inline              fvec4       xzxx                    (void)                                  const       noexcept;
        inline              fvec4       xzxy                    (void)                                  const       noexcept;
        inline              fvec4       xzxz                    (void)                                  const       noexcept;
        inline              fvec4       xzyx                    (void)                                  const       noexcept;
        inline              fvec4       xzyy                    (void)                                  const       noexcept;
        inline              fvec4       xzyz                    (void)                                  const       noexcept;
        inline              fvec4       xzzx                    (void)                                  const       noexcept;
        inline              fvec4       xzzy                    (void)                                  const       noexcept;
        inline              fvec4       xzzz                    (void)                                  const       noexcept;
        inline              fvec4       yxxx                    (void)                                  const       noexcept;
        inline              fvec4       yxxy                    (void)                                  const       noexcept;
        inline              fvec4       yxxz                    (void)                                  const       noexcept;
        inline              fvec4       yxyx                    (void)                                  const       noexcept;
        inline              fvec4       yxyy                    (void)                                  const       noexcept;
        inline              fvec4       yxyz                    (void)                                  const       noexcept;
        inline              fvec4       yxzx                    (void)                                  const       noexcept;
        inline              fvec4       yxzy                    (void)                                  const       noexcept;
        inline              fvec4       yxzz                    (void)                                  const       noexcept;
        inline              fvec4       yyxx                    (void)                                  const       noexcept;
        inline              fvec4       yyxy                    (void)                                  const       noexcept;
        inline              fvec4       yyxz                    (void)                                  const       noexcept;
        inline              fvec4       yyyx                    (void)                                  const       noexcept;
        inline              fvec4       yyyy                    (void)                                  const       noexcept;
        inline              fvec4       yyyz                    (void)                                  const       noexcept;
        inline              fvec4       yyzx                    (void)                                  const       noexcept;
        inline              fvec4       yyzy                    (void)                                  const       noexcept;
        inline              fvec4       yyzz                    (void)                                  const       noexcept;
        inline              fvec4       yzxx                    (void)                                  const       noexcept;
        inline              fvec4       yzxy                    (void)                                  const       noexcept;
        inline              fvec4       yzxz                    (void)                                  const       noexcept;
        inline              fvec4       yzyx                    (void)                                  const       noexcept;
        inline              fvec4       yzyy                    (void)                                  const       noexcept;
        inline              fvec4       yzyz                    (void)                                  const       noexcept;
        inline              fvec4       yzzx                    (void)                                  const       noexcept;
        inline              fvec4       yzzy                    (void)                                  const       noexcept;
        inline              fvec4       yzzz                    (void)                                  const       noexcept;
        inline              fvec4       zxxx                    (void)                                  const       noexcept;
        inline              fvec4       zxxy                    (void)                                  const       noexcept;
        inline              fvec4       zxxz                    (void)                                  const       noexcept;
        inline              fvec4       zxyx                    (void)                                  const       noexcept;
        inline              fvec4       zxyy                    (void)                                  const       noexcept;
        inline              fvec4       zxyz                    (void)                                  const       noexcept;
        inline              fvec4       zxzx                    (void)                                  const       noexcept;
        inline              fvec4       zxzy                    (void)                                  const       noexcept;
        inline              fvec4       zxzz                    (void)                                  const       noexcept;
        inline              fvec4       zyxx                    (void)                                  const       noexcept;
        inline              fvec4       zyxy                    (void)                                  const       noexcept;
        inline              fvec4       zyxz                    (void)                                  const       noexcept;
        inline              fvec4       zyyx                    (void)                                  const       noexcept;
        inline              fvec4       zyyy                    (void)                                  const       noexcept;
        inline              fvec4       zyyz                    (void)                                  const       noexcept;
        inline              fvec4       zyzx                    (void)                                  const       noexcept;
        inline              fvec4       zyzy                    (void)                                  const       noexcept;
        inline              fvec4       zyzz                    (void)                                  const       noexcept;
        inline              fvec4       zzxx                    (void)                                  const       noexcept;
        inline              fvec4       zzxy                    (void)                                  const       noexcept;
        inline              fvec4       zzxz                    (void)                                  const       noexcept;
        inline              fvec4       zzyx                    (void)                                  const       noexcept;
        inline              fvec4       zzyy                    (void)                                  const       noexcept;
        inline              fvec4       zzyz                    (void)                                  const       noexcept;
        inline              fvec4       zzzx                    (void)                                  const       noexcept;
        inline              fvec4       zzzy                    (void)                                  const       noexcept;
        inline              fvec4       zzzz                    (void)                                  const       noexcept;

        // Operator overloads
        inline              fvec3_t     operator+               (const fvec3_t& other)                  const       noexcept;
        inline              fvec3_t     operator-               (const fvec3_t& other)                  const       noexcept;
        inline              fvec3_t     operator*               (const fvec3_t& other)                  const       noexcept;
        inline              fvec3_t     operator*               (float scalar)                          const       noexcept;
        inline              fvec3_t     operator/               (float scalar)                          const       noexcept;
        inline              fvec3_t&    operator+=              (const fvec3_t& other)                              noexcept;
        inline              fvec3_t&    operator-=              (const fvec3_t& other)                              noexcept;
        inline              fvec3_t&    operator*=              (float scalar)                                      noexcept;
        inline              fvec3_t&    operator/=              (float scalar)                                      noexcept;
        inline              bool        operator==              (const fvec3_t& other)                  const       noexcept;
        inline              bool        operator!=              (const fvec3_t& other)                  const       noexcept;
        inline              float       operator[]              (std::int32_t index)                    const       noexcept;
        inline              float&      operator[]              (std::int32_t index)                                noexcept;

        // Friend operators
        template<bool> friend fvec3_t   operator*               (float scalar, const fvec3_t& v)                    noexcept;
        template<bool> friend fvec3_t   operator-               (const fvec3_t& v)                                  noexcept;
    };
}
