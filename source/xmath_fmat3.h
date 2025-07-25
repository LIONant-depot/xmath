#pragma once
#ifndef XMATH_FLOAT_LINEAR_ALGEBRA_H
    #error "You should not include this header directly, just need to include xmath_flinear.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // details::f_mat3
    //------------------------------------------------------------------------------
    //
    // Internal data layouts for 3x3 float matrix (column-major storage).
    //
    // Notes:
    //  simd_data: Padded to 3x4 for SSE alignment/efficiency (floatx4 columns).
    //  cpu_data: Compact 3x3 without padding.
    //  Unions provide flexible access: columns, cells (col-major), flat elements, or named members (m_rowcol).
    //  Aligned to 16 bytes for SIMD variant.
    //
    namespace details::f_mat3
    {
        struct alignas(16) simd_data
        {
            union
            {
                floatx4                             m_Columns[3];
                std::array<std::array<float, 4>, 3> m_Cells;
                std::array<float, 12>               m_Elements;
                struct
                {
                    float m_00, m_10, m_20, m_30;
                    float m_01, m_11, m_21, m_31;
                    float m_02, m_12, m_22, m_32;
                };
            };
        };

        struct cpu_data
        {
            union
            {
                std::array<std::array<float, 3>, 3> m_Cells;
                std::array<float, 9>                m_Elements;
                struct
                {
                    float m_00, m_10, m_20;
                    float m_01, m_11, m_21;
                    float m_02, m_12, m_22;
                };
            };
        };
    }

    //------------------------------------------------------------------------------
    // fmat3_t
    //------------------------------------------------------------------------------
    //
    // Templated 3x3 float matrix class (column-major, optional SIMD via SSE).
    //
    // Notes:
    //  Template<bool t_use_simd_v="">: true for SIMD (padded columns), false for CPU (compact).
    //  Internal storage column-major; accessors feel row-major (op()(row,col)).
    //  Supports affine transforms (rotation, scale); no translation (use fmat4 for 3D homogeneous).
    //  Does not initialize if default constructed; use fromIdentity/setupIdentity.
    //  Prioritizes performance: mutable ops chainable, shorter names.
    //  Immutable ops suffixed Copy (e.g., RotateCopy).
    //  Asserts finite matrices in ops; isFinite for checks.
    //  Constexpr where possible; inline implementations.
    //  Assumes right-hand rule for rotations; extracts to quats/vecs.
    //
    template <bool T_USE_SIMD_V>
    struct fmat3_t : std::conditional_t<T_USE_SIMD_V, details::f_mat3::simd_data, details::f_mat3::cpu_data>
    {
        // Constructors
        constexpr                       fmat3_t                 (void)                                              noexcept = default;
        constexpr                       fmat3_t                 (float diagonal)                                    noexcept;
        constexpr                       fmat3_t                 (const std::array<float, 9>& arr)                   noexcept;
        constexpr                       fmat3_t                 (std::span<const float, 9> span)                    noexcept;
        inline                          fmat3_t                 (const fquat& q)                                    noexcept;
        inline                          fmat3_t                 (const radian3& Euler)                              noexcept;
        inline                          fmat3_t                 (const fquat& rotation, const fvec3& scale)         noexcept;
        constexpr explicit              fmat3_t                 (const fmat3_t<!T_USE_SIMD_V>& other)               noexcept;

        // Static constructors
        static constexpr fmat3_t        fromIdentity            (void)                                              noexcept;
        static constexpr fmat3_t        fromZero                (void)                                              noexcept;
        static inline    fmat3_t        fromRotation            (const fquat& q)                                    noexcept;
        static inline    fmat3_t        fromRotation            (const radian3& Euler)                              noexcept;
        static inline    fmat3_t        fromRotation            (const fvec3& axis, radian angle)                   noexcept;
        static inline    fmat3_t        fromScale               (const fvec3& s)                                    noexcept;
        static inline    fmat3_t        fromRotationX           ( radian Angle )                                    noexcept;
        static inline    fmat3_t        fromRotationY           ( radian Angle )                                    noexcept;
        static inline    fmat3_t        fromRotationZ           ( radian Angle )                                    noexcept;

        // Setup methods (mutable, replace)
        inline fmat3_t&                 setup                   (const fquat& rotation, const fvec3& scale)         noexcept;
        inline fmat3_t&                 setupIdentity           (void)                                              noexcept;
        inline fmat3_t&                 setupZero               (void)                                              noexcept;
        inline fmat3_t&                 setupRotation           (const fquat& q)                                    noexcept;
        inline fmat3_t&                 setupRotation           (const radian3& euler)                              noexcept;
        inline fmat3_t&                 setupScale              (const fvec3& s)                                    noexcept;
        inline fmat3_t&                 setupScale              (float s)                                           noexcept;

        // Accessors (transposed for row-major feel: operator()(row, col) accesses m_Cells[col][row])
        constexpr fvec3                 operator[]              (size_t row)                                const   noexcept;
        constexpr float&                operator()              (size_t row, size_t col)                            noexcept;
        constexpr const float&          operator()              (size_t row, size_t col)                    const   noexcept;
        constexpr                       operator std::span<const float,9> ()                                const   noexcept;
        inline                          operator fquat          ()                                          const   noexcept;

        // Operations
        inline fmat3_t                  operator+               (const fmat3_t& other)                      const   noexcept;
        inline fmat3_t                  operator-               (const fmat3_t& other)                      const   noexcept;
        inline fmat3_t                  operator*               (const fmat3_t& other)                      const   noexcept;
        inline fmat3_t&                 operator+=              (const fmat3_t& other)                              noexcept;
        inline fmat3_t&                 operator-=              (const fmat3_t& other)                              noexcept;
        inline fmat3_t&                 operator*=              (const fmat3_t& other)                              noexcept;
        inline fvec3                    operator*               (const fvec3& v)                            const   noexcept;
        inline bool                     Equals                  (const fmat3_t& other, float tolerance)     const   noexcept;

        // Math functions
        inline fmat3_t                  Transpose               (void)                                      const   noexcept;
        inline fmat3_t                  Inverse                 (void)                                      const   noexcept;
        inline float                    Determinant             (void)                                      const   noexcept;
        inline fmat3_t&                 Orthogonalize           (void)                                              noexcept;

        // Geometry helpers
        inline fquat                    ExtractRotation         (void)                                      const   noexcept;
        inline fvec3                    ExtractScale            (void)                                      const   noexcept;
        inline fvec3                    Forward                 (void)                                      const   noexcept;
        inline fvec3                    Back                    (void)                                      const   noexcept;
        inline fvec3                    Up                      (void)                                      const   noexcept;
        inline fvec3                    Down                    (void)                                      const   noexcept;
        inline fvec3                    Left                    (void)                                      const   noexcept;
        inline fvec3                    Right                   (void)                                      const   noexcept;
        inline fvec3                    RotateVector            (const fvec3& v)                            const   noexcept;
        inline fvec3                    InvRotateVector         (const fvec3& v)                            const   noexcept;
        inline fvec3                    TransformDirection      (const fvec3& d)                            const   noexcept;

        // Mutable chaining methods (post-multiply)
        inline fmat3_t&                 Rotate                  (const fquat& q)                                    noexcept;
        inline fmat3_t&                 Rotate                  (const fvec3& axis, radian angle)                   noexcept;
        inline fmat3_t&                 RotateX                 (radian angle)                                      noexcept;
        inline fmat3_t&                 RotateY                 (radian angle)                                      noexcept;
        inline fmat3_t&                 RotateZ                 (radian angle)                                      noexcept;
        inline fmat3_t&                 Scale                   (const fvec3& s)                                    noexcept;
        inline fmat3_t&                 Scale                   (float s)                                           noexcept;

        // Pre-multiply mutable chaining
        inline fmat3_t&                 PreRotate               (const fquat& q)                                    noexcept;
        inline fmat3_t&                 PreRotate               (const fvec3& axis, radian angle)                   noexcept;
        inline fmat3_t&                 PreRotateX              (radian angle)                                      noexcept;
        inline fmat3_t&                 PreRotateY              (radian angle)                                      noexcept;
        inline fmat3_t&                 PreRotateZ              (radian angle)                                      noexcept;
        inline fmat3_t&                 PreScale                (const fvec3& s)                                    noexcept;
        inline fmat3_t&                 PreScale                (float s)                                           noexcept;

        // Clear methods
        inline fmat3_t&                 ClearRotation           (void)                                              noexcept;
        inline fmat3_t&                 ClearScale              (void)                                              noexcept;

        // Immutable versions (Copy suffix)
        inline fmat3_t                  RotateCopy              (const fquat& q)                            const   noexcept;
        inline fmat3_t                  RotateCopy              (const fvec3& axis, radian angle)           const   noexcept;
        inline fmat3_t                  ScaleCopy               (const fvec3& s)                            const   noexcept;
        inline fmat3_t                  PreRotateCopy           (const fquat& q)                            const   noexcept;
        inline fmat3_t                  PreRotateCopy           (const fvec3& axis, radian angle)           const   noexcept;
        inline fmat3_t                  PreScaleCopy            (const fvec3& s)                            const   noexcept;

        // Safety and validation
        inline bool                     isFinite                (void)                                      const   noexcept;
        inline bool                     isIdentity              (void)                                      const   noexcept;
        inline void                     SanityCheck             (void)                                      const   noexcept;
    };
}