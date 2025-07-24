#ifndef XMATH_FLOAT_LINEAR_ALGEBRA_H
#define XMATH_FLOAT_LINEAR_ALGEBRA_H
#pragma once

#include "xmath_fvector.h"

//------------------------------------------------------------------------------
// pre-definitions
//------------------------------------------------------------------------------
namespace xmath
{
    template< bool T_USE_SIMD_V>
    struct fquat_t;
    using fquat     = fquat_t<true>;
    using fquatd    = fquat_t<false>;

    template< bool T_USE_SIMD_V>
    struct fmat4x3_t;
    using fmat4x3  = fmat4x3_t<true>;
    using fmat4x3d = fmat4x3_t<false>;

    template< bool T_USE_SIMD_V>
    struct fmat4_t;
    using fmat4  = fmat4_t<true>;
    using fmat4d = fmat4_t<false>;

    template< bool T_USE_SIMD_V>
    struct fmat3_t;
    using fmat3  = fmat3_t<true>;
    using fmat3d = fmat3_t<false>;
}

//------------------------------------------------------------------------------
// classes
//------------------------------------------------------------------------------
#include "xmath_fquat.h"
#include "xmath_fmat4.h"
#include "xmath_fmat3.h"
#include "xmath_ftransform.h"

//------------------------------------------------------------------------------
// inline functions
//------------------------------------------------------------------------------
#include "implementation/xmath_fquat_inline.h"
#include "implementation/xmath_fmat4_inline.h"
#include "implementation/xmath_fmat3_inline.h"
#include "implementation/xmath_ftransform_inline.h"

#endif