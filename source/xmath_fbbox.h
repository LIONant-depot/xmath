#pragma once
#ifndef XMATH_FLOAT_SHAPES_H
    #error "You should not include this header directly, just need to include xmath_fshapes.h"
#endif

namespace xmath
{
    //------------------------------------------------------------------------------
    // supported data structure for fbbox_t:
    //------------------------------------------------------------------------------
    namespace details::f_bbox
    {
        struct alignas(16) simd_data
        {
            union
            {
                std::array<float, 4+4>  m_Elements;
                floatx4                 m_MinMax[2];
                struct
                {
                    fvec3 m_Min;
                    fvec3 m_Max;
                };
            };
        };
        struct cpu_data
        {
            union
            {
                std::array<float, 6>    m_Elements;
                struct
                {
                    fvec3d m_Min;
                    fvec3d m_Max;
                };
            };
        };
    }

    //------------------------------------------------------------------------------
    // fbbox_t
    //------------------------------------------------------------------------------
    //
    // 3D axis-aligned bounding box class with SIMD optimization (SSE).
    //
    // Notes:
    // This class does not initialize its memory if default constructed.
    // Aligned to 16 bytes for SIMD; non-SIMD version is compact.
    // Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    // Immutable/copy ops suffixed _copy.
    // Use const for safety; asserts in inline for validity (isfinite components).
    // Targets C++20: constexpr/consteval where possible.
    // Assumes min <= max for valid bounds; min/max are axis-aligned.
    //
    template <bool T_USE_SIMD_V>
    struct fbbox_t : std::conditional_t<T_USE_SIMD_V, details::f_bbox::simd_data, details::f_bbox::cpu_data>
    {
        using parent_t = std::conditional_t<T_USE_SIMD_V, details::f_bbox::simd_data, details::f_bbox::cpu_data>;

        // Constructors
        constexpr                           fbbox_t                     (void)                                                          noexcept = default;
        constexpr                           fbbox_t                     (const fvec3& Point)                                            noexcept;
        constexpr                           fbbox_t                     (const fvec3& Min, const fvec3& Max)                            noexcept;
        inline                              fbbox_t                     (const std::span<const fvec3> Verts )                           noexcept;
        constexpr                           fbbox_t                     (const fvec3& Center, float Radius)                             noexcept;
        constexpr explicit                  fbbox_t                     (const std::array<floatx4, 2>& MinMax)                          noexcept requires T_USE_SIMD_V;
        constexpr                           fbbox_t                     (const fbbox_t<!T_USE_SIMD_V>& Other)                           noexcept;
        constexpr                           fbbox_t                     (const std::array<float, 6>& Conversion)                        noexcept;

        // Assignment and conversion operators
        constexpr                           operator std::array<float, 6> (void)                                                const   noexcept;
        inline                              operator std::string        (void)                                                  const   noexcept;
        inline              std::string     toString                    (void)                                                  const   noexcept;
        template <bool V>
        inline friend       std::ostream&   operator<<                  (std::ostream& Os, const fbbox_t<V>& BBox)                      noexcept;

        // Static properties
        static constexpr    fbbox_t         fromZero                    (void)                                                          noexcept;
        static constexpr    fbbox_t         fromIdentity                (void)                                                          noexcept;
        static inline       fbbox_t         fromPoint                   (const fvec3& Point)                                            noexcept;
        static inline       fbbox_t         fromMinMax                  (const fvec3& Min, const fvec3& Max)                            noexcept;
        static inline       fbbox_t         fromCenterRadius            (const fvec3& Center, float Radius)                             noexcept;
        static inline       fbbox_t         fromVerts                   (const std::span<const fvec3> Verts )                           noexcept;

        // Static methods
        static inline       bool            Intersect                   (const fbbox_t& BBox1, const fbbox_t& BBox2)                    noexcept;
        static inline       bool            Intersect                   (const fbbox_t& BBox, const fvec3& Point)                       noexcept;
        static inline       bool            Intersect                   (const fbbox_t& BBox, const fplane& Plane)                      noexcept;
        static inline       bool            Intersect                   (const fbbox_t& BBox, float& T, const fvec3& P0, const fvec3& P1) noexcept;
        static inline       bool            Intersect                   (const fbbox_t& BBox, const std::span<const fvec3> Verts)       noexcept;
        static inline       bool            IntersectTriangle           (const fbbox_t& BBox, const fvec3& P0, const fvec3& P1, const fvec3& P2) noexcept;
        static inline       fbbox_t         Intersection                (const fbbox_t& BBox1, const fbbox_t& BBox2)                    noexcept;

        // Static methods as members
        inline              bool            Intersect                   (const fbbox_t& Other)                                  const   noexcept;
        inline              bool            Intersect                   (const fvec3& Point)                                    const   noexcept;
        inline              bool            Intersect                   (const fplane& Plane)                                   const   noexcept;
        inline              bool            Intersect                   (float& T, const fvec3& P0, const fvec3& P1)            const   noexcept;
        inline              bool            Intersect                   (const std::span<fvec3> Verts)                          const   noexcept;
        inline              bool            IntersectTriangle           (const fvec3& P0, const fvec3& P1, const fvec3& P2)     const   noexcept;
        inline              fbbox_t         Intersection                (const fbbox_t& Other)                                  const   noexcept;

        // Instance methods - Basic operations
        inline              bool            isFinite                    (void)                                                  const   noexcept;
        inline              bool            isValid                     (void)                                                  const   noexcept;
        inline              fvec3           getSize                     (void)                                                  const   noexcept;
        inline              fvec3           getCenter                   (void)                                                  const   noexcept;
        inline              float           getRadius                   (void)                                                  const   noexcept;
        inline              float           getRadiusSquared            (void)                                                  const   noexcept;
        inline              float           getSurfaceArea              (void)                                                  const   noexcept;
        inline              bool            isInRange                   (float Min, float Max)                                  const   noexcept;
        inline              bool            Contains                    (const fbbox_t& Other)                                  const   noexcept;
        inline              bool            ContainsPoint               (const fvec3& Point)                                    const   noexcept;
        inline              fvec3           getClosestPoint             (const fvec3& Point)                                    const   noexcept;
        inline              void            getVerts                    (std::span<fvec3> Dst)                                  const   noexcept;

        // Instance methods - Setup operations
        inline              fbbox_t&        setupFromPoint              (const fvec3& Point)                                            noexcept;
        inline              fbbox_t&        setupFromMinMax             (const fvec3& Min, const fvec3& Max)                            noexcept;
        inline              fbbox_t&        setupFromCenterRadius       (const fvec3& Center, float Radius)                             noexcept;
        inline              fbbox_t&        setupFromVerts              (const std::span<const fvec3> Verts)                            noexcept;
        inline              fbbox_t&        setZero                     (void)                                                          noexcept;
        inline              fbbox_t&        setIdentity                 (void)                                                          noexcept;
        inline              fbbox_t&        Inflate                     (const fvec3& Delta)                                            noexcept;
        inline              fbbox_t&        Deflate                     (const fvec3& Delta)                                            noexcept;
        inline              fbbox_t&        Translate                   (const fvec3& Delta)                                            noexcept;
        inline              fbbox_t&        Transform                   (const fmat4& Matrix)                                           noexcept;
        inline              fbbox_t         TransformCopy               (const fmat4& Matrix)                                   const   noexcept;

        // Operator overloads
        inline              fbbox_t&        operator+=                  (const fbbox_t& Other)                                          noexcept;
        inline              fbbox_t&        operator+=                  (const fvec3& Point)                                            noexcept;
        inline              fbbox_t         operator+                   (const fbbox_t& Other)                                  const   noexcept;
        inline              fbbox_t         operator+                   (const fvec3& Point)                                    const   noexcept;
        inline              bool            operator==                  (const fbbox_t& Other)                                  const   noexcept;
        inline              bool            operator!=                  (const fbbox_t& Other)                                  const   noexcept;
        inline              float           operator[]                  (std::int32_t Index)                                    const   noexcept;
        inline              float&          operator[]                  (std::int32_t Index)                                            noexcept;

        // Friend operators
        template<bool V> friend fbbox_t     operator+                   (const fvec3& Point, const fbbox_t& BBox)                       noexcept;
    };
}