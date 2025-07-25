#pragma once
#ifndef XMATH_FLOAT_LINEAR_ALGEBRA_H
    #error "You should not include this header directly, just need to include xmath_flinear.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    // details::f_mat4
    //------------------------------------------------------------------------------
    //
    // Internal data layouts for 4x4 float matrix (column-major storage).
    //
    // Notes:
    //  simd_data: 4x4 with SSE alignment (floatx4 columns).
    //  cpu_data: Compact 4x4.
    //  Unions provide flexible access: columns, cells (col-major), flat elements, or named members (m_rowcol).
    //  Aligned to 16 bytes for SIMD variant.
    //
    namespace details::f_mat4
    {
        struct alignas(16) simd_data
        {
            union
            {
                floatx4                             m_Columns[4];
                std::array<std::array<float, 4>, 4> m_Cells;
                std::array<float, 16>               m_Elements;
                struct
                {
                    float m_00, m_10, m_20, m_30;
                    float m_01, m_11, m_21, m_31;
                    float m_02, m_12, m_22, m_32;
                    float m_03, m_13, m_23, m_33;
                };
            };
        };

        struct cpu_data
        {
            union
            {
                std::array<std::array<float, 4>, 4> m_Cells;
                std::array<float, 16>               m_Elements;
                struct
                {
                    float m_00, m_10, m_20, m_30;
                    float m_01, m_11, m_21, m_31;
                    float m_02, m_12, m_22, m_32;
                    float m_03, m_13, m_23, m_33;
                };
            };
        };
    }

    //------------------------------------------------------------------------------
    // fmat4_t
    //------------------------------------------------------------------------------
    //
    // Templated 4x4 float matrix class (column-major, optional SIMD via SSE).
    //
    // Notes:
    //  Template<bool t_use_simd_v="">: true for SIMD, false for CPU.
    //  Internal storage column-major; accessors feel row-major (op()(row,col)).
    //  Supports 3D transformations (translation, rotation, scale, projection); homogeneous coords.
    //  Does not initialize if default constructed; use fromIdentity/setupIdentity.
    //  Prioritizes performance: mutable ops chainable (&#x26; return), shorter names.
    //  Immutable ops suffixed Copy (e.g., RotateCopy).
    //  Asserts finite matrices in ops; isFinite for checks.
    //  Constexpr where possible; inline implementations.
    //  Assumes right-hand rule for rotations; extracts to quats/vecs.
    //  Projection matrices for graphics (perspective, ortho, look-at, billboard).
    //
    template <bool T_USE_SIMD_V>
    struct fmat4_t : std::conditional_t<T_USE_SIMD_V, details::f_mat4::simd_data, details::f_mat4::cpu_data>
    {
        // Constructors
        constexpr                       fmat4_t                 (void)                                              noexcept = default;
        constexpr                       fmat4_t                 (float diagonal)                                    noexcept;
        constexpr                       fmat4_t                 (const std::array<float, 16>& arr)                  noexcept;
        constexpr                       fmat4_t                 (std::span<const float, 16> span)                   noexcept;
        inline                          fmat4_t                 (const fquat& q)                                    noexcept;
        inline                          fmat4_t                 (const radian3& Euler)                              noexcept;
        inline                          fmat4_t                 (const fvec3& translation, const fquat& rotation, const fvec3& scale) noexcept;
        constexpr explicit              fmat4_t                 (const fmat4_t<!T_USE_SIMD_V>& other)               noexcept;

        // Static constructors
        static constexpr fmat4_t        fromIdentity            (void)                                              noexcept;
        static constexpr fmat4_t        fromZero                (void)                                              noexcept;
        static inline    fmat4_t        fromTranslation         (const fvec3& t)                                    noexcept;
        static inline    fmat4_t        fromRotation            (const fquat& q)                                    noexcept;
        static inline    fmat4_t        fromRotation            (const fvec3& axis, radian angle)                   noexcept;
        static inline    fmat4_t        fromRotation            (const radian3& Euler)                              noexcept;
        static inline    fmat4_t        fromScale               (const fvec3& s)                                    noexcept;
        static inline    fmat4_t        fromPerspective         (radian fov, float aspect, float near_plane, float far_plane) noexcept;
        static inline    fmat4_t        fromPerspective         (float left, float right, float bottom, float top, float near_plane, float far_plane) noexcept;
        static inline    fmat4_t        fromOrtho               (float left, float right, float bottom, float top, float near_plane, float far_plane) noexcept;
        static inline    fmat4_t        fromOrtho               (float width, float height, float near_plane, float far_plane) noexcept;
        static inline    fmat4_t        fromLookAt              (const fvec3& eye, const fvec3& target, const fvec3& up) noexcept;
        static inline    fmat4_t        fromBillboard           (const fvec3& from, const fvec3& to, const fvec3& up) noexcept;
        static inline    fmat4_t        fromRotationX           ( radian Angle )                                    noexcept;
        static inline    fmat4_t        fromRotationY           ( radian Angle )                                    noexcept;
        static inline    fmat4_t        fromRotationZ           ( radian Angle )                                    noexcept;

        // Setup methods (mutable, replace)
        inline fmat4_t&                 setup                   (const fvec3& translation, const fquat& rotation, const fvec3& scale) noexcept;
        inline fmat4_t&                 setupIdentity           (void)                                              noexcept;
        inline fmat4_t&                 setupZero               (void)                                              noexcept;
        inline fmat4_t&                 setupTranslation        (const fvec3& t)                                    noexcept;
        inline fmat4_t&                 setupRotation           (const fquat& q)                                    noexcept;
        inline fmat4_t&                 setupRotation           (const radian3& euler)                              noexcept;
        inline fmat4_t&                 setupScale              (const fvec3& s)                                    noexcept;
        inline fmat4_t&                 setupScale              (float s)                                           noexcept;

        // Accessors (transposed for row-major feel: operator()(row, col) accesses m_Cells[col][row])
        constexpr fvec4                 operator[]              (size_t row)                                const   noexcept;
        constexpr float&                operator()              (size_t row, size_t col)                            noexcept;
        constexpr const float&          operator()              (size_t row, size_t col)                    const   noexcept;
        constexpr                       operator std::span<const float,16> ()                               const   noexcept;
        inline                          operator fquat          ()                                          const   noexcept;

        // Operations
        inline fmat4_t                  operator+               (const fmat4_t& other)                      const   noexcept;
        inline fmat4_t                  operator-               (const fmat4_t& other)                      const   noexcept;
        inline fmat4_t                  operator*               (const fmat4_t& other)                      const   noexcept;
        inline fmat4_t&                 operator+=              (const fmat4_t& other)                              noexcept;
        inline fmat4_t&                 operator-=              (const fmat4_t& other)                              noexcept;
        inline fmat4_t&                 operator*=              (const fmat4_t& other)                              noexcept;
        inline fvec4                    operator*               (const fvec4& v)                            const   noexcept;
        inline fvec3                    operator*               (const fvec3& v)                            const   noexcept;
        inline bool                     Equals                  (const fmat4_t& other, float tolerance)     const   noexcept;

        // Math functions
        inline fmat4_t                  Transpose               (void)                                      const   noexcept;
        inline fmat4_t                  Inverse                 (void)                                      const   noexcept;
        inline fmat4_t                  InverseSRT              (void)                                      const   noexcept;
        inline fmat4_t                  InverseRT               (void)                                      const   noexcept;
        inline float                    Determinant             (void)                                      const   noexcept;
        inline fmat4_t&                 Orthogonalize           (void)                                              noexcept;

        // Geometry helpers
        inline fvec3                    ExtractPosition         (void)                                      const   noexcept;
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
        inline fvec3                    TransformPosition       (const fvec3& p)                            const   noexcept;
        inline fvec3                    TransformDirection      (const fvec3& d)                            const   noexcept;

        // Mutable chaining methods (post-multiply)
        inline fmat4_t&                 Translate               (const fvec3& t)                                    noexcept;
        inline fmat4_t&                 Rotate                  (const fquat& q)                                    noexcept;
        inline fmat4_t&                 Rotate                  (const fvec3& axis, radian angle)                   noexcept;
        inline fmat4_t&                 RotateX                 (radian angle)                                      noexcept;
        inline fmat4_t&                 RotateY                 (radian angle)                                      noexcept;
        inline fmat4_t&                 RotateZ                 (radian angle)                                      noexcept;
        inline fmat4_t&                 Scale                   (const fvec3& s)                                    noexcept;
        inline fmat4_t&                 Scale                   (float s)                                           noexcept;

        // Pre-multiply mutable chaining
        inline fmat4_t&                 PreTranslate            (const fvec3& t)                                    noexcept;
        inline fmat4_t&                 PreRotate               (const fquat& q)                                    noexcept;
        inline fmat4_t&                 PreRotate               (const fvec3& axis, radian angle)                   noexcept;
        inline fmat4_t&                 PreRotateX              (radian angle)                                      noexcept;
        inline fmat4_t&                 PreRotateY              (radian angle)                                      noexcept;
        inline fmat4_t&                 PreRotateZ              (radian angle)                                      noexcept;
        inline fmat4_t&                 PreScale                (const fvec3& s)                                    noexcept;
        inline fmat4_t&                 PreScale                (float s)                                           noexcept;

        // Clear methods
        inline fmat4_t&                 ClearTranslation        (void)                                              noexcept;
        inline fmat4_t&                 ClearRotation           (void)                                              noexcept;
        inline fmat4_t&                 ClearScale              (void)                                              noexcept;

        // Immutable versions (Copy suffix)
        inline fmat4_t                  TranslateCopy           (const fvec3& t)                            const   noexcept;
        inline fmat4_t                  RotateCopy              (const fquat& q)                            const   noexcept;
        inline fmat4_t                  RotateCopy              (const fvec3& axis, radian angle)           const   noexcept;
        inline fmat4_t                  ScaleCopy               (const fvec3& s)                            const   noexcept;
        inline fmat4_t                  PreTranslateCopy        (const fvec3& t)                            const   noexcept;
        inline fmat4_t                  PreRotateCopy           (const fquat& q)                            const   noexcept;
        inline fmat4_t                  PreRotateCopy           (const fvec3& axis, radian angle)           const   noexcept;
        inline fmat4_t                  PreScaleCopy            (const fvec3& s)                            const   noexcept;

        // Safety and validation
        inline bool                     isFinite                (void)                                      const   noexcept;
        inline bool                     isIdentity              (void)                                      const   noexcept;
        inline void                     SanityCheck             (void)                                      const   noexcept;
    };
}