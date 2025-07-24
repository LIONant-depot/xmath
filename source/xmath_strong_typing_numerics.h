#ifndef XMATH_STRONG_TYPE_NUMBERICS_H
#define XMATH_STRONG_TYPE_NUMBERICS_H
#pragma once
namespace xmath
{
    //------------------------------------------------------------------------------
    // strong_typing_numerics_t (CRTP base for numeric strong types)
    //------------------------------------------------------------------------------
    template <typename T_DERIVED, typename T>
    struct strong_typing_numerics_t
    {
        T m_Value;

        constexpr               strong_typing_numerics_t(void)                              noexcept = default;
        constexpr explicit      strong_typing_numerics_t(T value)                           noexcept : m_Value(value) {}

        static consteval T_DERIVED  fromZero    (void)                                      noexcept { return T_DERIVED{ T{0} }; }

        // Arithmetic operators
        constexpr T_DERIVED         operator+   (const T_DERIVED rhs)               const   noexcept { return T_DERIVED{ m_Value + rhs.m_Value }; }
        constexpr T_DERIVED         operator-   (const T_DERIVED rhs)               const   noexcept { return T_DERIVED{ m_Value - rhs.m_Value }; }
        constexpr T_DERIVED         operator*   (T scalar)                          const   noexcept { return T_DERIVED{ m_Value * scalar }; }
        constexpr T_DERIVED         operator/   (T scalar)                          const   noexcept { return T_DERIVED{ m_Value / scalar }; }
        friend constexpr T_DERIVED  operator*   (T scalar, T_DERIVED rhs)                   noexcept { return rhs * scalar; }

        // Compound assignment
        constexpr T_DERIVED&        operator+=  (const T_DERIVED rhs)                       noexcept { m_Value += rhs.m_Value; return static_cast<T_DERIVED&>(*this); }
        constexpr T_DERIVED&        operator-=  (const T_DERIVED rhs)                       noexcept { m_Value -= rhs.m_Value; return static_cast<T_DERIVED&>(*this); }
        constexpr T_DERIVED&        operator*=  (T scalar)                                  noexcept { m_Value *= scalar; return static_cast<T_DERIVED&>(*this); }
        constexpr T_DERIVED&        operator/=  (T scalar)                                  noexcept { m_Value /= scalar; return static_cast<T_DERIVED&>(*this); }

        // Unary operators
        constexpr T_DERIVED         operator-   ()                                  const   noexcept requires (std::is_unsigned_v<T> == false) { return T_DERIVED{ -m_Value }; }

        // Comparison operators via spaceship
        constexpr auto              operator<=> (const strong_typing_numerics_t& rhs) const noexcept = default;
    };
}
#endif