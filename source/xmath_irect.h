#ifndef XMATH_IMATH_H
    #error "You should not include this header directly, just need to include xmath_imath.h"
#endif
#pragma once

namespace xmath
{
    //------------------------------------------------------------------------------
    // irect
    //------------------------------------------------------------------------------
    //
    // 2D axis-aligned integer bounding rectangle class.
    //
    // Notes:
    // This class does not initialize its memory if default constructed.
    // Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    // Immutable/copy ops suffixed _copy.
    // Use const for safety; asserts in inline for validity (min <= max).
    // Targets C++20: constexpr/consteval where possible.
    // Assumes min <= max for valid bounds; min/max are axis-aligned.
    //
    struct irect
    {
        union
        {
            std::array<int, 4> m_Elements;
            struct
            {
                ivec2 m_Min;
                ivec2 m_Max;
            };
        };

        // Constructors
        constexpr                               irect                       (void)                                                              noexcept = default;
        constexpr                               irect                       (const ivec2& Point)                                                noexcept;
        constexpr                               irect                       (const ivec2& Min, const ivec2& Max)                                noexcept;
        inline                                  irect                       (const std::span<const ivec2> Verts)                                noexcept;
        constexpr                               irect                       (const ivec2& Center, int Radius)                                   noexcept;
        constexpr                               irect                       (const std::array<int, 4>& Conversion)                              noexcept;
        constexpr                               irect                       (int L, int T, int R, int B)                                        noexcept;

        // Assignment and conversion operators
        constexpr           operator std::array<int, 4>                     (void)                                                      const   noexcept;
        inline              operator std::string                            (void)                                                      const   noexcept;
        inline              std::string         toString                    (void)                                                      const   noexcept;
        inline friend       std::ostream&       operator<<                  (std::ostream& Os, const irect& Rect)                               noexcept;

        // Static properties
        static constexpr    irect               fromZero                    (void)                                                              noexcept;
        static constexpr    irect               fromIdentity                (void)                                                              noexcept;
        static inline       irect               fromPoint                   (const ivec2& Point)                                                noexcept;
        static inline       irect               fromMinMax                  (const ivec2& Min, const ivec2& Max)                                noexcept;
        static inline       irect               fromCenterRadius            (const ivec2& Center, int Radius)                                   noexcept;
        static inline       irect               fromVerts                   (const std::span<const ivec2> Verts)                                noexcept;
        // Static methods
        static inline       bool                Intersect                   (const irect& Rect1, const irect& Rect2)                            noexcept;
        static inline       bool                Intersect                   (const irect& Rect, const ivec2& Point)                             noexcept;
        static inline       bool                Intersect                   (const irect& Rect, float& T, const ivec2& P0, const ivec2& P1)     noexcept;
        static inline       bool                Intersect                   (const irect& Rect, const std::span<const ivec2> Verts)             noexcept;
        static inline       irect               Intersection                (const irect& Rect1, const irect& Rect2)                            noexcept;

        // Static methods as members
        inline              bool                Intersect                   (const irect& Other)                                        const   noexcept;
        inline              bool                Intersect                   (const ivec2& Point)                                        const   noexcept;
        inline              bool                Intersect                   (float& T, const ivec2& P0, const ivec2& P1)                const   noexcept;
        inline              bool                Intersect                   (const std::span<ivec2> Verts)                              const   noexcept;
        inline              irect               Intersection                (const irect& Other)                                        const   noexcept;

        // Instance methods - Basic operations
        inline              bool                isFinite                    (void)                                                      const   noexcept;
        inline              bool                isValid                     (void)                                                      const   noexcept;
        inline              ivec2               getSize                     (void)                                                      const   noexcept;
        inline              ivec2               getCenter                   (void)                                                      const   noexcept;
        inline              float               getRadius                   (void)                                                      const   noexcept;
        inline              long long           getRadiusSquared            (void)                                                      const   noexcept;
        inline              int                 getArea                     (void)                                                      const   noexcept;
        inline              bool                isInRange                   (int Min, int Max)                                          const   noexcept;
        inline              bool                Contains                    (const irect& Other)                                        const   noexcept;
        inline              bool                ContainsPoint               (const ivec2& Point)                                        const   noexcept;
        inline              ivec2               getClosestPoint             (const ivec2& Point)                                        const   noexcept;
        inline              void                getVerts                    (std::span<ivec2> Dst)                                      const   noexcept;
        inline              int                 getWidth                    (void)                                                      const   noexcept;
        inline              int                 getHeight                   (void)                                                      const   noexcept;
        inline              float               Difference                  (const irect& Other)                                        const   noexcept;
        inline              bool                isEmpty                     (void)                                                      const   noexcept;

        // Instance methods - Setup operations
        inline              irect&              setupFromPoint              (const ivec2& Point)                                                noexcept;
        inline              irect&              setupFromMinMax             (const ivec2& Min, const ivec2& Max)                                noexcept;
        inline              irect&              setupFromCenterRadius       (const ivec2& Center, int Radius)                                   noexcept;
        inline              irect&              setupFromVerts              (const std::span<const ivec2> Verts)                                noexcept;
        inline              irect&              setup                       (int L, int T, int R, int B)                                        noexcept;
        inline              irect&              setZero                     (void)                                                              noexcept;
        inline              irect&              setIdentity                 (void)                                                              noexcept;

        inline              irect&              AddVerts                    (const std::span<const ivec2> Verts)                                noexcept;
        inline              irect&              Inflate                     (const ivec2& Delta)                                                noexcept;
        inline              irect&              Deflate                     (const ivec2& Delta)                                                noexcept;
        inline              irect&              Translate                   (const ivec2& Delta)                                                noexcept;
        inline              irect&              setWidth                    (int W)                                                             noexcept;
        inline              irect&              setHeight                   (int H)                                                             noexcept;
        inline              irect&              setSize                     (int W, int H)                                                      noexcept;
        inline              irect&              setMax                      (void)                                                              noexcept;

        // Operator overloads
        inline              irect&              operator+=                  (const irect& Other)                                                noexcept;
        inline              irect&              operator+=                  (const ivec2& Point)                                                noexcept;
        inline              irect               operator+                   (const irect& Other)                                        const   noexcept;
        inline              irect               operator+                   (const ivec2& Point)                                        const   noexcept;
        inline              bool                operator==                  (const irect& Other)                                        const   noexcept;
        inline              bool                operator!=                  (const irect& Other)                                        const   noexcept;
        inline              int                 operator[]                  (int Index)                                                 const   noexcept;
        inline              int&                operator[]                  (int Index)                                                         noexcept;

        // Friend operators
        inline friend       irect               operator+                   (const ivec2& Point, const irect& Rect)                             noexcept;
    };
}