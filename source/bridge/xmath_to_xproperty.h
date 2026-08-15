#ifndef XMATH_TO_XPROPERTY_H
#define XMATH_TO_XPROPERTY_H
#pragma once

namespace xmath
{
    // This header is meant to be includable on its own (from my_properties.h itself, or from any
    // descriptor file, regardless of what that caller has already brought into scope) - so it can't
    // rely on an ambient "using namespace xproperty" the way a normal descriptor file happens to have.
    // Scoped to this namespace block rather than file scope to keep the pollution contained.
    using namespace xproperty;
    using namespace xproperty::settings;

    struct vec2_friend;
    struct vec3_friend;
    struct fbbox_friend;
}

// These specializations must come before ANY friend struct body below is defined/registered.
// fbbox_friend's own Min/Max members are typed xmath::fvec3 - registering fbbox_friend forces the
// property system to resolve reflected_type<fvec3> right then, at fbbox_friend's own XPROPERTY_REG.
// If that specialization isn't visible yet, the compiler implicitly instantiates the primary (identity)
// reflected_type<fvec3> template to answer that, and a later explicit specialization for the same type
// becomes ill-formed ("explicit specialization after instantiation") - this bit the fvec3 case exactly
// when these specializations lived after the namespace closed, past fbbox_friend's definition.
template<> struct xproperty::settings::reflected_type<xmath::fvec2>  { using type = xmath::vec2_friend; };
template<> struct xproperty::settings::reflected_type<xmath::fvec3>  { using type = xmath::vec3_friend; };
template<> struct xproperty::settings::reflected_type<xmath::fbbox>  { using type = xmath::fbbox_friend; };

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