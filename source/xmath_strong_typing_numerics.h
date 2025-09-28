#ifndef XMATH_STRONG_TYPE_NUMBERICS_H
#define XMATH_STRONG_TYPE_NUMBERICS_H
#pragma once
#include <type_traits>
#include <compare>
namespace xmath
{
    //------------------------------------------------------------------------------
    // strong_typing_numerics_t (CRTP base for numeric strong types)
    // This version requires the user to make this class a parent class
    // Example:
    //      struct degree : strong_typing_numerics_t<degree, T>
    //      {
    //              // Here the user is able to override or add any functions 
    //      }
    //
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
        constexpr T_DERIVED         operator*   (const T_DERIVED rhs)               const   noexcept { return T_DERIVED{ m_Value * rhs.m_Value }; }
        constexpr T_DERIVED         operator/   (const T_DERIVED rhs)               const   noexcept { return T_DERIVED{ m_Value / rhs.m_Value }; }
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

    //------------------------------------------------------------------------------
    // strong_typing_direct_numerics_t (CRTP base for numeric strong types)
    // This version the user can not override anything...
    // Example:
    //      using degree = strong_typing_direct_numerics_t<degree, struct degree_tag>;
    //
    //------------------------------------------------------------------------------
    template <typename T, typename T_TAG >
    struct strong_typing_direct_numerics_t
    {
        T m_Value;
        using class_t = strong_typing_direct_numerics_t<T, T_TAG>;

        constexpr                   strong_typing_direct_numerics_t(void)                   noexcept = default;
        constexpr explicit          strong_typing_direct_numerics_t(T value)                noexcept : m_Value(value) {}

        static consteval class_t    fromZero    (void)                                      noexcept { return class_t{ T{0} }; }

        // Arithmetic operators
        constexpr class_t           operator+   (const class_t rhs)                 const   noexcept { return class_t{ m_Value + rhs.m_Value }; }
        constexpr class_t           operator-   (const class_t rhs)                 const   noexcept { return class_t{ m_Value - rhs.m_Value }; }
        constexpr class_t           operator*   (const class_t rhs)                 const   noexcept { return class_t{ m_Value * rhs.m_Value }; }
        constexpr class_t           operator/   (const class_t rhs)                 const   noexcept { return class_t{ m_Value / rhs.m_Value }; }
        constexpr class_t           operator*   (T scalar)                          const   noexcept { return class_t{ m_Value * scalar }; }
        constexpr class_t           operator/   (T scalar)                          const   noexcept { return class_t{ m_Value / scalar }; }
        friend constexpr class_t    operator*   (T scalar, class_t rhs)                     noexcept { return rhs * scalar; }

        // Compound assignment
        constexpr class_t&          operator+=  (const class_t rhs)                         noexcept { m_Value += rhs.m_Value; return static_cast<class_t&>(*this); }
        constexpr class_t&          operator-=  (const class_t rhs)                         noexcept { m_Value -= rhs.m_Value; return static_cast<class_t&>(*this); }
        constexpr class_t&          operator*=  (T scalar)                                  noexcept { m_Value *= scalar; return static_cast<class_t&>(*this); }
        constexpr class_t&          operator/=  (T scalar)                                  noexcept { m_Value /= scalar; return static_cast<class_t&>(*this); }

        // Unary operators
        constexpr class_t           operator-   ()                                  const   noexcept requires (std::is_unsigned_v<T> == false) { return class_t{ -m_Value }; }

        // Comparison operators via spaceship
        constexpr auto              operator<=> (const class_t& rhs)                const   noexcept = default;
    };
}

//
// Make our strong types compatible with std::maps...
//
namespace std
{
    template <typename T_DERIVED, typename T>
    struct hash<xmath::strong_typing_numerics_t<T_DERIVED, T>>
    {
        inline
        std::size_t operator()(const xmath::strong_typing_numerics_t<T_DERIVED, T>& k) const noexcept
        {
            return std::hash<T>{}(k.m_Value);
        }
    };

    template <typename T, typename T_TAG >
    struct hash<xmath::strong_typing_direct_numerics_t<T, T_TAG>>
    {
        inline
        std::size_t operator()(const xmath::strong_typing_direct_numerics_t<T, T_TAG>& k) const noexcept
        {
            return std::hash<T>{}(k.m_Value);
        }
    };
}

#endif