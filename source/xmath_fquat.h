#pragma once
#ifndef XMATH_FLOAT_LINEAR_ALGEBRA_H
    #error "You should not include this header directly, just need to include xmath_flinear.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // supported data structure for fquat_t:
    //------------------------------------------------------------------------------
    namespace details::f_quat
    {
        struct alignas(16) simd_data
        {
            union
            {
                floatx4             m_XYZW;
                std::array<float,4> m_Elements;
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
                std::array<float, 4> m_Elements;
                struct
                {
                    float m_X, m_Y, m_Z, m_W;
                };
            };
        };
    }

    //------------------------------------------------------------------------------
    // fquat_t
    //------------------------------------------------------------------------------
    //
    // 4D quaternion class with SIMD optimization (SSE).
    //
    // Notes:
    //  This class does not initialize its memory if default constructed.
    //  Aligned to 16 bytes for SIMD.
    //  Prioritizes performance; mutable ops common (shorter names, chainable via & return).
    //  Immutable/copy ops suffixed _copy.
    //  Use const for safety; asserts in inline for validity (isfinite components).
    //  Targets C++20: constexpr\consteval where possible.
    //  No swizzle methods as quaternions require all components for meaningful operations.
    //  Assumes right-hand rule for rotations.
    //
    template <bool T_USE_SIMD_V >
    struct fquat_t : std::conditional_t< T_USE_SIMD_V, details::f_quat::simd_data, details::f_quat::cpu_data >
    {
        using parent_t = std::conditional_t< T_USE_SIMD_V, details::f_quat::simd_data, details::f_quat::cpu_data >;

        // Constructors
        constexpr                       fquat_t                (void)                                           noexcept = default;
        constexpr                       fquat_t                (float x, float y, float z, float w)             noexcept;
        constexpr                       fquat_t                (const fvec3& axis, radian angle)                noexcept;
        inline                          fquat_t                (const radian3& euler)                           noexcept;
        inline                          fquat_t                (const fvec3& from, const fvec3& to, const fvec3& up ) noexcept;
        inline                          fquat_t                (const fvec3& forward, const fvec3& up ) noexcept;
        constexpr explicit              fquat_t                (const floatx4& reg)                             noexcept requires T_USE_SIMD_V;
        constexpr                       fquat_t                (const fquat_t<!T_USE_SIMD_V>& other)            noexcept;
        constexpr                       fquat_t                (const std::array<double, 4>& conversion)        noexcept;

        // Assignment and conversion operators
        constexpr                       operator std::array<double, 4> ()                               const   noexcept;
        inline                          operator std::string    ()                                      const   noexcept;
        std::string                     ToString                (void)                                  const   noexcept;
        template <bool V >
        inline friend std::ostream&     operator<<              (std::ostream& os, const fquat_t<V>& quat)      noexcept;

        // Static properties
        static constexpr    fquat_t     fromIdentity            (void)                                          noexcept;
        static constexpr    fquat_t     fromZero                (void)                                          noexcept;
        static inline       fquat_t     fromAxisAngle           ( const fvec3& axis, radian angle )             noexcept;

        // Static methods
        static inline       float       Dot                     (const fquat_t& a, const fquat_t& b)            noexcept;
        static inline       fquat_t     Lerp                    (const fquat_t& a, const fquat_t& b, float t)   noexcept;
        static inline       fquat_t     LerpFast                (const fquat_t& Start, const fquat_t& End, float T) noexcept;
        static inline       fquat_t     LerpUnclamped           (const fquat_t& a, const fquat_t& b, float t )  noexcept;
        static inline       fquat_t     Slerp                   (const fquat_t& a, const fquat_t& b, float t)   noexcept;
        static inline       fquat_t     SlerpAccurate           (const fquat_t& a, const fquat_t& b, float t)   noexcept;
        static inline       fquat_t     SlerpUnclamped          (const fquat_t& a, const fquat_t& b, float t)   noexcept;
        static inline       fquat_t     Squad                   (const fquat_t& a, const fquat_t& a_tangent, const fquat_t& b, const fquat_t& b_tangent, float t) noexcept;
        static inline       radian      AngleBetween            (const fquat_t& a, const fquat_t& b)        noexcept;
        static inline       fquat_t     FromToRotation          (const fvec3& from, const fvec3& to, const fvec3& up = fvec3::fromUp())        noexcept;
        static inline       fquat_t     LookRotation            (const fvec3& forward, const fvec3& up)         noexcept;
        static inline       fquat_t     RandomUnitQuaternion    (void)                                          noexcept;

        // Static methods as members
        inline              float       Dot                     (const fquat_t& other)                  const   noexcept;
        inline              fquat_t     Lerp                    (const fquat_t& other, float t)         const   noexcept;
        inline              fquat_t     LerpFast                (const fquat_t& other, float t)         const   noexcept;
        inline              fquat_t     LerpUnclamped           (const fquat_t& other, float t)         const   noexcept;
        inline              fquat_t     Slerp                   (const fquat_t& other, float t)         const   noexcept;
        inline              fquat_t     SlerpAccurate           (const fquat_t& other, float t)         const   noexcept;
        inline              fquat_t     SlerpUnclamped          (const fquat_t& other, float t)         const   noexcept;
        inline              fquat_t     Squad                   (const fquat_t& a_tangent, const fquat_t& b, const fquat_t& b_tangent, float t) const noexcept;
        inline              radian      AngleBetween            (const fquat_t& other)                  const   noexcept;

        // Instance methods - Basic operations
        inline              float       Length                  (void)                                  const   noexcept;
        inline              float       LengthSq                (void)                                  const   noexcept;
        inline              fquat_t     NormalizeCopy           (void)                                  const   noexcept;
        inline              fquat_t&    Normalize               (void)                                          noexcept;
        inline              fquat_t     NormalizeSafeCopy       (void)                                  const   noexcept;
        inline              fquat_t&    NormalizeSafe           (void)                                          noexcept;
        inline              bool        isFinite                (void)                                  const   noexcept;
        inline              bool        isNormalized            (float tolerance = 1e-6f)               const   noexcept;
        inline              bool        isNearlyIdentity        (float tolerance = 1e-6f)               const   noexcept;
        inline              bool        isNearlyZero            (float tolerance = 1e-6f)               const   noexcept;
        inline              bool        Equals                  (const fquat_t& other, float tolerance = 1e-6f) const noexcept;

        // Instance methods - Quaternion specifics
        inline              fquat_t     ConjugateCopy           (void)                                  const   noexcept;
        inline              fquat_t&    Conjugate               (void)                                          noexcept;
        inline              fquat_t     InverseCopy             (void)                                  const   noexcept;
        inline              fquat_t&    Inverse                 (void)                                          noexcept;
        inline              fvec3       Axis                    (void)                                  const   noexcept;
        inline              radian      Angle                   (void)                                  const   noexcept;
        inline              radian3     ToEuler                 (void)                                  const   noexcept;
        inline              fvec3       Forward                 (void)                                  const   noexcept;
        inline              fvec3       Up                      (void)                                  const   noexcept;
        inline              fvec3       Right                   (void)                                  const   noexcept;
        inline std::pair<fvec3, radian> ToAxisAngle             (void)                                  const   noexcept;
        inline              fquat_t     Delta                   (const fquat_t& other)                  const   noexcept;
        inline              fquat_t     LogCopy                 (void)                                  const   noexcept;
        inline              fquat_t&    Log                     (void)                                          noexcept;
        inline              fquat_t     ExpCopy                 (void)                                  const   noexcept;
        inline              fquat_t&    Exp                     (void)                                          noexcept;

        // Instance methods - Rotation operations
        inline              fquat_t&    setupRotationX          (radian rx)                                     noexcept;
        inline              fquat_t&    setupRotationY          (radian ry)                                     noexcept;
        inline              fquat_t&    setupRotationZ          (radian rz)                                     noexcept;
        inline              fquat_t&    setupLookRotation       (const fvec3& forward, const fvec3& up)         noexcept;
        inline              fquat_t&    setupFromToRotation     (const fvec3& from, const fvec3& to, const fvec3& up) noexcept;
        inline              fquat_t     RotateXCopy             (radian rx)                             const   noexcept;
        inline              fquat_t     RotateYCopy             (radian ry)                             const   noexcept;
        inline              fquat_t     RotateZCopy             (radian rz)                             const   noexcept;
        inline              fquat_t&    RotateX                 (radian rx)                                     noexcept;
        inline              fquat_t&    RotateY                 (radian ry)                                     noexcept;
        inline              fquat_t&    RotateZ                 (radian rz)                                     noexcept;
        inline              fquat_t&    PreRotateX              (radian rx)                                     noexcept;
        inline              fquat_t&    PreRotateY              (radian ry)                                     noexcept;
        inline              fquat_t&    PreRotateZ              (radian rz)                                     noexcept;
        inline              fquat_t     PreRotateXCopy          (radian rx)                             const   noexcept;
        inline              fquat_t     PreRotateYCopy          (radian ry)                             const   noexcept;
        inline              fquat_t     PreRotateZCopy          (radian rz)                             const   noexcept;
        inline              fquat_t     RotateTowardsCopy       (const fquat_t& target, radian maxDelta) const  noexcept;
        inline              fquat_t&    RotateTowards           (const fquat_t& target, radian maxDelta)        noexcept;
        inline              fvec3       RotateVector            (const fvec3& v)                        const   noexcept;

        // Operator overloads
        inline              fquat_t     operator+               (const fquat_t& other)                  const   noexcept;
        inline              fquat_t     operator-               (const fquat_t& other)                  const   noexcept;
        inline              fquat_t     operator*               (float scalar)                          const   noexcept;
        inline              fquat_t     operator/               (float scalar)                          const   noexcept;
        inline              fquat_t     operator*               (const fquat_t& other)                  const   noexcept;
        inline              fquat_t&    operator+=              (const fquat_t& other)                          noexcept;
        inline              fquat_t&    operator-=              (const fquat_t& other)                          noexcept;
        inline              fquat_t&    operator*=              (float scalar)                                  noexcept;
        inline              fquat_t&    operator/=              (float scalar)                                  noexcept;
        inline              fquat_t&    operator*=              (const fquat_t& other)                          noexcept;
        inline              bool        operator==              (const fquat_t& other)                  const   noexcept;
        inline              bool        operator!=              (const fquat_t& other)                  const   noexcept;
        inline              float       operator[]              (std::int32_t index)                    const   noexcept;
        inline              float&      operator[]              (std::int32_t index)                            noexcept;

        // Friend operators
        template<bool> friend fquat_t   operator*               (float scalar, const fquat_t& q)                noexcept;
        template<bool> friend fquat_t   operator-               (const fquat_t& q)                              noexcept;
        template<bool> friend fvec3     operator*               (const fquat_t& q, const fvec3& v)              noexcept;
    };
}