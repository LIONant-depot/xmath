#ifndef _XMATH_RADIAN3_H
#define _XMATH_RADIAN3_H
#pragma once

#ifndef XMATH_H
    #include "xmath.h"
#endif

namespace xmath
{
    //------------------------------------------------------------------------------
    // Description:
    //     This class represents a 3D rotation using the classical Roll, Pitch, Yaw.
    //     The order of rotation there for is (Z,X,Y)
    //     It is very useful when trying to minimize the memory footprint needed to 
    //     represent a rotation or if you just want to deal with only angles.
    //
    //<P>  The class is not hardware accelerated and works only with floats. There
    //     is not need for special alignments except for the atomic types.
    // See Also:
    //     radian xquaternion matrix4
    //------------------------------------------------------------------------------
    struct radian3 
    {
        constexpr                               radian3                 ( void )                                                                        noexcept = default; 
        constexpr                               radian3                 ( const radian Angles )                                                         noexcept;
        constexpr                               radian3                 ( radian Pitch, radian Yaw, radian Roll )                                       noexcept;

        static consteval    radian3             fromZero                ( void )                                                                        noexcept;

        inline              void                Zero                    ( void )                                                                        noexcept;
        inline              radian3&            ModAngle                ( void )                                                                        noexcept;
        inline              radian3             ModAngleCopy            ( void )                                                                const   noexcept;
        inline              radian3&            ModAngle2               ( void )                                                                        noexcept;
        inline              radian3             ModAngle2Copy           ( void )                                                                const   noexcept;
        inline              radian3             MinAngleDiff            ( const radian3& R )                                                    const   noexcept;

        constexpr           radian              Difference              ( const radian3& R )                                                    const   noexcept;
        constexpr           bool                isInrange               ( radian Min, radian Max )                                              const   noexcept;
        inline              bool                isValid                 ( void )                                                                const   noexcept;
        inline              bool                hasGimbalLock           ( void )                                                                const   noexcept;
        inline              radian3&            Normalize               ( void )                                                                        noexcept;
        inline              radian3             NormalizeCopy           ( void )                                                                const   noexcept;
        inline              radian3&            Clamp                   ( radian Min, radian Max )                                                      noexcept;
        inline              radian3             ClampCopy               ( radian Min, radian Max )                                              const   noexcept;
        inline              radian3             Lerp                    ( float t, const radian3& Other )                                       const   noexcept;

        constexpr           bool                operator ==             ( const radian3& R )                                                    const   noexcept;
        inline              const radian3&      operator +=             ( const radian3& R )                                                            noexcept;
        inline              const radian3&      operator -=             ( const radian3& R )                                                            noexcept;
        inline              const radian3&      operator *=             ( float Scalar )                                                                noexcept;
        inline              const radian3&      operator /=             ( float Scalar )                                                                noexcept;

        inline              radian3             operator +              ( const radian3& R )                                                    const   noexcept;
        inline              radian3             operator -              ( const radian3& R )                                                    const   noexcept;
        inline              radian3             operator -              ( void )                                                                const   noexcept;
        inline              radian3             operator /              ( float S )                                                             const   noexcept;
        inline              radian3             operator *              ( float S )                                                             const   noexcept;
        inline friend       radian3             operator *              ( float S, const radian3& R )                                                   noexcept;


        radian     m_Pitch;
        radian     m_Yaw;
        radian     m_Roll;
    };
}

#include "implementation/xmath_radian3_inline.h"
#endif