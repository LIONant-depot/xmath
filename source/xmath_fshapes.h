#ifndef XMATH_FLOAT_SHAPES_H
#define XMATH_FLOAT_SHAPES_H
#pragma once

#include "xmath_flinear.h"

//------------------------------------------------------------------------------
// pre-definitions
//------------------------------------------------------------------------------
namespace xmath
{
    template< bool T_USE_SIMD_V>
    struct fplane_t;
    using fplane  = fplane_t<true>;
    using fplaned = fplane_t<false>;


    template< bool T_USE_SIMD_V>
    struct fbbox_t;
    using fbbox  = fbbox_t<true>;
    using fbboxd = fbbox_t<false>;

}

#include "xmath_fplane.h"
#include "xmath_fbbox.h"




#include "implementation/xmath_fplane_inline.h"
#include "implementation/xmath_fbbox_inline.h"

#endif