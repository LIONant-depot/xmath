#ifndef XMATH_TO_XPROPERTY_H
#define XMATH_TO_XPROPERTY_H
#pragma once

namespace xmath
{
    struct vec2_friend : xmath::fvec2
    {
        XPROPERTY_DEF
        ( "vector2", xmath::fvec2, xproperty::settings::vector2_group
        , obj_member<"X", &xmath::fvec2::m_X >
        , obj_member<"Y", &xmath::fvec2::m_Y >
        )
    };
    XPROPERTY_REG(vec2_friend)

    //------------------------------------------------------------------------------------------------

    struct vec3_friend : xmath::fvec3
    {
        XPROPERTY_DEF
        ( "vector3", xmath::fvec3, xproperty::settings::vector3_group
        , obj_member<"X", &xmath::fvec3::m_X >
        , obj_member<"Y", &xmath::fvec3::m_Y >
        , obj_member<"Z", &xmath::fvec3::m_Z >
        )
    };
    XPROPERTY_REG(vec3_friend)

    //------------------------------------------------------------------------------------------------

    struct fbbox_friend : xmath::fbbox
    {
        XPROPERTY_DEF
        ( "fbbox", xmath::fbbox
        , obj_member<"Min", &xmath::fbbox::m_Min >
        , obj_member<"Max", &xmath::fbbox::m_Max >
        )
    };
    XPROPERTY_REG(fbbox_friend)

}

#endif