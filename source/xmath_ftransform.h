#pragma once
#ifndef XMATH_FLOAT_LINEAR_ALGEBRA_H
    #error "You should not include this header directly, just need to include xmath_flinear.h"
#endif
namespace xmath
{
    //------------------------------------------------------------------------------
    struct transform3
    {
        static inline   transform3      Blend           (const transform3& From, const transform3& To, const float T)   noexcept;
        inline          transform3&     Blend           (const transform3& To, const float T)                           noexcept;
        inline          void            setupIdentity   (void)                                                          noexcept;
        inline          fmat4           toMatrix        (void)                                                  const   noexcept;
        inline                          operator fmat4  (void)                                                  const   noexcept;

        fvec3           m_Scale     {fvec3::fromOne()};
        fquat           m_Rotation  {fquat::fromIdentity()};
        fvec3           m_Position  {fvec3::fromZero()};
    };

    //------------------------------------------------------------------------------
    struct transform2
    {
        static inline   transform2      Blend           (const transform2& From, const transform2& To, const float T)   noexcept;
        inline          transform2&     Blend           (const transform2& To, const float T)                           noexcept;
        inline          void            setupIdentity   (void)                                                          noexcept;

        fvec2           m_Scale     {fvec2::fromOne()};
        radian          m_Rotation  {radian::fromZero()};
        fvec2           m_Position  {fvec2::fromZero()};
    };
}
