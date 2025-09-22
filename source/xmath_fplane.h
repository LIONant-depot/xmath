#pragma once
#ifndef XMATH_FLOAT_SHAPES_H
    #error "You should not include this header directly, just need to include xmath_fshapes.h"
#endif

namespace xmath
{
    //------------------------------------------------------------------------------
    // supported data structure for fplane_t:
    //------------------------------------------------------------------------------
    namespace details::f_plane
    {
        struct alignas(16) simd_data
        {
            union
            {
                xmath::floatx4          m_XYZD;
                std::array<float,4>     m_Elements;
                xmath::fvec3            m_Normal;
                struct
                {
                    float m_X, m_Y, m_Z, m_D;
                };
            };
        };
        struct cpu_data
        {
            union
            {
                std::array<float, 4>    m_Elements;
                xmath::fvec3d           m_Normal;
                struct
                {
                    float m_X, m_Y, m_Z, m_D;
                };
            };
        };
    }

    //------------------------------------------------------------------------------
    // fplane_t
    //------------------------------------------------------------------------------
    //
    // 4D plane class (normal XYZ, offset D) with SIMD optimization (SSE).
    //
    // Notes:
    // This class does not initialize its memory if default constructed.
    // Aligned to 16 bytes for SIMD.
    // Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    // Immutable/copy ops suffixed _copy.
    // Use const for safety; asserts in inline for validity (isfinite components).
    // Targets C++20: constexpr/consteval where possible.
    // Assumes plane equation: X*x + Y*y + Z*z + D = 0, with normalized normal (X^2 + Y^2 + Z^2 = 1).
    //
    template <bool T_USE_SIMD_V >
    struct fplane_t : std::conditional_t< T_USE_SIMD_V, details::f_plane::simd_data, details::f_plane::cpu_data >
    {
        using parent_t = std::conditional_t< T_USE_SIMD_V, details::f_plane::simd_data, details::f_plane::cpu_data >;

        // Constructors
        constexpr                           fplane_t                    (void)                                                          noexcept = default;
        constexpr                           fplane_t                    (float X, float Y, float Z, float D)                            noexcept;
        constexpr                           fplane_t                    (const fvec3& Normal, float Distance)                           noexcept;
        inline                              fplane_t                    (const fvec3& Normal, const fvec3& Point)                       noexcept;
        inline                              fplane_t                    (const fvec3& P1, const fvec3& P2, const fvec3& P3)             noexcept;
        constexpr explicit                  fplane_t                    (const floatx4& Reg)                                            noexcept requires T_USE_SIMD_V;
        constexpr                           fplane_t                    (const fplane_t<!T_USE_SIMD_V>& Other)                          noexcept;
        constexpr                           fplane_t                    (const std::array<double, 4>& Conversion)                       noexcept;
        constexpr                           fplane_t                    (const std::array<float, 4>& Conversion)                        noexcept;
        // Assignment and conversion operators
        constexpr                           operator std::array<float, 4>   (void)                                              const   noexcept;
        constexpr                           operator std::array<double, 4>  (void)                                              const   noexcept;
        inline                              operator std::string            (void)                                              const   noexcept;
        inline              std::string     ToString                        (void)                                              const   noexcept;
        template <bool V >
        inline friend       std::ostream&   operator<<                  ( std::ostream& Os, const fplane_t<V>& Plane)                   noexcept;

        // Static properties
        static constexpr    fplane_t        fromZero                    (void)                                                          noexcept;
        static constexpr    fplane_t        fromIdentity                (void)                                                          noexcept;
        static inline       fplane_t        fromNormalDistance          (const fvec3& Normal, float Distance)                           noexcept;
        static inline       fplane_t        fromNormalPoint             (const fvec3& Normal, const fvec3& Point)                       noexcept;
        static inline       fplane_t        fromThreePoints             (const fvec3& P1, const fvec3& P2, const fvec3& P3)             noexcept;

        // Static methods
        static inline       bool            IntersectThreePlanes        (const fplane_t& P1, const fplane_t& P2, const fplane_t& P3, fvec3& Result) noexcept;
        static inline       float           Dot                         (const fplane_t& Plane, const fvec3& Point)                     noexcept;
        static inline       float           Distance                    (const fplane_t& Plane, const fvec3& Point)                     noexcept;
        static inline       std::int32_t    Side                        (const fplane_t& Plane, const fvec3& Point)                     noexcept;
        static inline       bool            IntersectLine               (const fplane_t& Plane, float& T, const fvec3& Start, const fvec3& Direction) noexcept;
        static inline       bool            IntersectLineSegment        (const fplane_t& Plane, float& T, const fvec3& P0, const fvec3& P1) noexcept;

        // Static methods as members
        inline              float           Dot                         (const fvec3& Point)                                    const   noexcept;
        inline              float           Distance                    (const fvec3& Point)                                    const   noexcept;
        inline              std::int32_t    Side                        (const fvec3& Point)                                    const   noexcept;
        inline              bool            IntersectLine               (float& T, const fvec3& Start, const fvec3& Direction)  const   noexcept;
        inline              bool            IntersectLineSegment        (float& T, const fvec3& P0, const fvec3& P1)            const   noexcept;

        // Instance methods - Basic operations
        inline              float           NormalLength                (void)                                                  const   noexcept;
        inline              float           NormalLengthSq              (void)                                                  const   noexcept;
        inline              fplane_t        NormalizeCopy               (void)                                                  const   noexcept;
        inline              fplane_t&       Normalize                   (void)                                                          noexcept;
        inline              fplane_t        NormalizeSafeCopy           (void)                                                  const   noexcept;
        inline              fplane_t&       NormalizeSafe               (void)                                                          noexcept;
        inline              bool            isFinite                    (void)                                                  const   noexcept;
        inline              bool            isNormalized                (float Tolerance = 1e-6f)                               const   noexcept;
        inline              bool            isNearlyZero                (float Tolerance = 1e-6f)                               const   noexcept;
        inline              bool            Equals                      (const fplane_t& Other, float Tolerance = 1e-6f)        const   noexcept;

        // Instance methods - Plane specifics
        inline              fvec3           Normal                      (void)                                                  const   noexcept;
        inline              float           Offset                      (void)                                                  const   noexcept;
        inline              fvec3           GetOrigin                   (void)                                                  const   noexcept;
        inline              fplane_t&       ComputeOffset               (const fvec3& Point)                                            noexcept;
        inline              void            DecomposeVector             (const fvec3& V, fvec3& Parallel, fvec3& Perpendicular) const   noexcept;
        inline              void            OrthoVectors                (fvec3& AxisA, fvec3& AxisB)                            const   noexcept;
        inline              fvec3           ReflectVector               (const fvec3& V)                                        const   noexcept;
        inline              fvec3           Project                     (const fvec3& Point)                                    const   noexcept;
        inline              bool            IsPointOnPlane              (const fvec3& Point, float Epsilon = 1e-5f)             const   noexcept;
        inline              bool            SameSide                    (const fvec3& P0, const fvec3& P1)                      const   noexcept;
        inline              bool            ClipNGon                    (fvec3* Dst, std::int32_t& DstCount, const fvec3* Src, std::int32_t SrcCount) const noexcept;

        // Instance methods - Setup operations
        inline              fplane_t&       setupFromNormalDistance     (const fvec3& Normal, float Distance)                           noexcept;
        inline              fplane_t&       setupFromNormalPoint        (const fvec3& Normal, const fvec3& Point)                       noexcept;
        inline              fplane_t&       setupFromThreePoints        (const fvec3& P1, const fvec3& P2, const fvec3& P3)             noexcept;
        inline              fplane_t        TranslateCopy               (const fvec3& Translation)                              const   noexcept;
        inline              fplane_t&       Translate                   (const fvec3& Translation)                                      noexcept;
        inline              fplane_t        FlipCopy                    (void)                                                  const   noexcept;
        inline              fplane_t&       Flip                        (void)                                                          noexcept;

        // Operator overloads
        inline              fplane_t        operator-                   (void)                                                  const   noexcept;
        inline              float           operator[]                  (std::int32_t Index)                                    const   noexcept;
        inline              float&          operator[]                  (std::int32_t Index)                                            noexcept;
        inline              bool            operator==                  (const fplane_t& Other)                                 const   noexcept;
        inline              bool            operator!=                  (const fplane_t& Other)                                 const   noexcept;

        // Friend operators
        template<bool V> friend fplane_t    operator*                   (const fmat4& M, const fplane_t& P)                             noexcept;
    };
}