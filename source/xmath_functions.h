#ifndef XMATH_FUNCTIONS_H
#define XMATH_FUNCTIONS_H
#pragma once

namespace xmath
{
    constexpr static auto       e_v         = std::numbers::e_v<float>;
    constexpr static auto       log2e_v     = std::numbers::log2e_v<float>;
    constexpr static auto       log10e_v    = std::numbers::log10e_v<float>;
    constexpr static auto       ln2_v       = std::numbers::ln2_v<float>;
    constexpr static auto       ln10_v      = std::numbers::ln10_v<float>;
    constexpr static auto       sqrt2_v     = std::numbers::sqrt2_v<float>;
    constexpr static auto       flt_tol_v   = 0.001f;

    template<std::floating_point T> inline      T               Exp             ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Pow             ( T a, T b )                                        noexcept;
    template<std::floating_point T> inline      T               FMod            ( T x, T y )                                        noexcept;
    template<std::floating_point T> inline      T               ModFX           ( T x, T& y )                                       noexcept;
    template<std::floating_point T> inline      T               Log             ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Log2            ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Log10           ( T x )                                             noexcept;
    template<std::floating_point T> constexpr   T               i2f             ( T i )                                             noexcept;
    template<std::floating_point T> constexpr   std::int32_t    f2i             ( T f )                                             noexcept;
    template <typename T>           constexpr   T               FSel            ( T a, T b, T c )                                   noexcept;
    template <typename T>           constexpr   T               Sqr             ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Sqrt            ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               InvSqrt         ( T x )                                             noexcept;
    template <typename T1, typename T2> constexpr auto          Min             ( T1 a, T2 b )                                      noexcept -> decltype(a + b);
    template <typename T1, typename T2> constexpr auto          Max             ( T1 a, T2 b )                                      noexcept -> decltype(a + b);
    template<std::floating_point T> inline      bool            FEqual          ( T f0, T f1, T tol = flt_tol_v )                   noexcept;
    template<std::floating_point T> constexpr   bool            FLess           ( T f0, T f1, T tol = flt_tol_v )                   noexcept;
    template<std::floating_point T> constexpr   bool            FGreater        ( T f0, T f1, T tol = flt_tol_v )                   noexcept;
    template <typename T> constexpr             bool            Sign            ( T x )                                             noexcept;
    template<std::floating_point T> inline      std::int32_t    LRound          ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Round           ( T a, T b )                                        noexcept;
    template <std::floating_point T> constexpr  T               Round           ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Ceil            ( T x )                                             noexcept;
    template<std::floating_point T> inline      T               Floor           ( T x )                                             noexcept;
    template <typename T>           constexpr   bool            isInRange       ( T x, T min, T max )                               noexcept;
    template <typename T>           constexpr   T               Range           ( T x, T min, T max )                               noexcept;
    template <typename T>           constexpr   T               Abs             ( T x )                                             noexcept;
    template <typename T>           constexpr   T               Clamp           ( const T& value, const T& low, const T& high )     noexcept;
    template <typename T>           constexpr   T               Lerp            ( float t, T a, T b )                               noexcept;
    template<std::floating_point T> inline      bool            isValid         ( T x )                                             noexcept;
    template<std::floating_point T> constexpr   T               Trunc           ( T x )                                             noexcept;
    template <std::floating_point T>constexpr   bool            isFinite        ( T x )                                             noexcept;
    template <std::floating_point T>constexpr   T               CopySign        ( T x, T y )                                        noexcept;

    inline bool                                                 SolvedQuadraticRoots(float& root1, float& root2, float a, float b, float c) noexcept;
}
#endif

