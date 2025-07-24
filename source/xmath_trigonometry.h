#ifndef XMATH_TRIGONOMETRY_H
#define XMATH_TRIGONOMETRY_H
#pragma once
//==============================================================================
//  TRIGONOMETRIC MATH FUNCTIONS
//==============================================================================
//
// Types:
//
//  degree - Typed wrapper for degree values.
//  radian - Typed wrapper for radian values (single precision).
//  radian64 - Typed wrapper for radian values (double precision).
//
// Functions:
//  
//  float    sin   ( radian angle   )        - Sine.
//  float    cos   ( radian angle   )        - Cosine.
//  float    tan   ( radian angle   )        - Tangent.
//  radian   asin  ( float  sine    )        - Arc sine.
//  radian   acos  ( float  cosine  )        - Arc cosine.
//  radian   atan  ( float  tangent )        - Arc tangent.
//  radian   atan2 ( float m_Y, float m_X    )   - Standard "atan2" arc tangent where m_Y can equal 0.
//                
// Additional functions:
//
//  void     sinCos       ( radian angle, float& sin, float& cos )      - Sine and cosine in one function call.
//  radian   modAngle     ( radian angle )                              - Provide equivalent angle in [    0, 360 ) degrees.
//  radian   modAngle2    ( radian angle )                              - Provide equivalent angle in [ -180, 180 ) degrees.
//  radian   minAngleDiff ( radian angle1, radian angle2 )              - Provide smallest angle between two given angles.
//
// Constants:
//
//  Use std::numbers for PI, E, etc., but provide aliases for convenience.
//
// Conversions:
//
//  degToRad / radToDeg templates for angle conversions.
//  User-defined literals: e.g., 90.0_xdeg.
//
//==============================================================================
namespace xmath
{
    //==============================================================================
    // ANGLE TYPES
    //==============================================================================
    template <std::floating_point T> constexpr T DegToRad(T Deg) noexcept { return Deg * std::numbers::pi_v<T> / T{ 180 }; }
    template <std::floating_point T> constexpr T RadToDeg(T Rad) noexcept { return Rad * T{ 180 } / std::numbers::pi_v<T>; }

    //------------------------------------------------------------------------------
    // degree
    //------------------------------------------------------------------------------
    template<std::floating_point T>
    struct radian_t;

    template<std::floating_point T>
    struct degree_t : strong_typing_numerics_t<degree_t<T>, T>
    {
        using parent_t = strong_typing_numerics_t<degree_t<T>, T>; using parent_t::parent_t;
        using parent_t::operator +;  using parent_t::operator -;  using parent_t::operator /;  using parent_t::operator *;
        using parent_t::operator +=; using parent_t::operator -=; using parent_t::operator /=; using parent_t::operator *=;
        constexpr auto  operator<=> (const degree_t& rhs) const noexcept = default;

        template<std::floating_point G>
        constexpr               degree_t(radian_t<G> deg)             noexcept : parent_t{ (static_cast<T>(RadToDeg(deg.m_Value))) } {}
        constexpr radian_t<T>   getRadians(void)              const   noexcept { return { DegToRad(this->m_Value) }; }
    };

    //------------------------------------------------------------------------------
    // radian
    //------------------------------------------------------------------------------
    template<std::floating_point T>
    struct radian_t : strong_typing_numerics_t<radian_t<T>, T>
    {
        using parent_t = strong_typing_numerics_t<radian_t<T>, T>; using parent_t::parent_t;
        using parent_t::operator +;  using parent_t::operator -;  using parent_t::operator /;  using parent_t::operator *;
        using parent_t::operator +=; using parent_t::operator -=; using parent_t::operator /=; using parent_t::operator *=;
        constexpr auto  operator<=> (const radian_t& rhs) const noexcept = default;

        template<std::floating_point G>
        constexpr               radian_t    ( degree_t<G> deg )             noexcept : parent_t{(static_cast<T>(DegToRad(deg.m_Value)))} {}
        constexpr degree_t<T>   getDegrees  ( void )                const   noexcept { return { RadToDeg(this->m_Value) }; }
    };

    //------------------------------------------------------------------------------
    // radian degree radian64 degree64
    //------------------------------------------------------------------------------
    using radian    = radian_t<float>;
    using degree    = degree_t<float>;
    using dradian   = radian_t<double>;
    using ddegree   = degree_t<double>;


    // User-defined literals
    constexpr degree   operator"" _xdeg(long double deg)          noexcept { return degree{ static_cast<float>(deg) }; }
    constexpr degree   operator"" _xdeg(unsigned long long deg)   noexcept { return degree{ static_cast<float>(deg) }; }

    constexpr radian   operator"" _xrad(long double rad)          noexcept { return radian{ static_cast<float>(rad) }; }
    constexpr radian   operator"" _xrad(unsigned long long rad)   noexcept { return radian{ static_cast<float>(rad) }; }

    constexpr dradian operator"" _xrad64(long double rad)        noexcept { return dradian{ static_cast<double>(rad) }; }
    constexpr dradian operator"" _xrad64(unsigned long long rad) noexcept { return dradian{ static_cast<double>(rad) }; }

    //==============================================================================
    // CONSTANTS (Aliases from std::numbers for convenience)
    //==============================================================================
    constexpr static radian     pi_v        = radian{ std::numbers::pi_v<float> };
    constexpr static dradian    pi64_v      = dradian{ std::numbers::pi_v<double> };
    constexpr static radian     pi2_v       = radian{ std::numbers::pi_v<float> *2 };
    constexpr static radian     pi_over2_v  = radian{ static_cast<float>(std::numbers::pi_v<double> / 2) };
    constexpr static radian     pi_over4_v  = radian{ std::numbers::pi_v<float> / 4 };

    //==============================================================================
    // BASIC FUNCTIONS
    //==============================================================================
    template<std::floating_point T> inline      radian_t<T>     ModAngle        ( radian_t<T> Angle )                               noexcept;
    template<std::floating_point T> inline      radian_t<T>     ModAngle2       ( radian_t<T> Angle )                               noexcept;
    template<std::floating_point T> inline      radian_t<T>     MinAngleDiff    ( radian_t<T> Angle1, radian_t<T> Angle2 )          noexcept;
    template<std::floating_point T> inline      radian_t<T>     LerpAngle       ( T t, radian_t<T> Angle1, radian_t<T> Angle2 )     noexcept;
    template<std::floating_point T> inline      void            SinCos          ( radian_t<T> Angle, T& S, T& C )                   noexcept;
    template<std::floating_point T> inline      T               Sin             ( radian_t<T> x )                                   noexcept;
    template<std::floating_point T> inline      T               Cos             ( radian_t<T> x )                                   noexcept;
    template<std::floating_point T> inline      T               Tan             ( radian_t<T> x )                                   noexcept;
    template<std::floating_point T> inline      radian_t<T>     ATan2           ( T y, T x )                                        noexcept;
    template<std::floating_point T> inline      radian_t<T>     ATan            ( T x )                                             noexcept;
    template<std::floating_point T> inline      radian_t<T>     Asin            ( T x )                                             noexcept;
    template<std::floating_point T> inline      radian_t<T>     Acos            ( T x )                                             noexcept;
}





#endif