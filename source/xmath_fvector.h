#ifndef XMATH_FLOAT_VECTOR_H
#define XMATH_FLOAT_VECTOR_H
#pragma once

#include <cstdint>      // for fixed-width integer types
#include <utility>      // for std::move, std::swap if used
#include <array>        // for std::array
#include <span>         // for std::span
#include <random>       // For std::mt19937, std::uniform_real_distribution, std::random_device
#include <cmath>        // For std::sqrt, std::cos, std::sin, std::fabs (used in your vector methods)
#include <format>       // Debug output
#include <ostream>


#ifndef XMATH_BASICS_H
    #include "xmath_basics.h"
#endif

#ifndef XMATH_RADIAN3_H
    #include "xmath_radian3.h"
#endif

//------------------------------------------------------------------------------
// pre-definitions
//------------------------------------------------------------------------------
namespace xmath
{
    template <bool T_USE_SIMD_V >
    struct  fvec3_t;
    using   fvec3     = fvec3_t<true>;      // * simd_data - is design for simd instructions which makes it very fast but it wastes 1 float, and requires a big alignment
    using   fvec3d    = fvec3_t<false>;     // * cpu_data  - is your standard 3 floats... but because is just standard it does not run as fast...
    using   floatx4   = __m128;             // xmath own alias for simd data
    struct  fvec4;                          // This class is simd optimize so alignment requirements are more aggressive than your standard floats
    struct  fvec2;                          // Pretty much your standard vector2

    namespace simde
    {
        enum : int { X = 0, Y = 1, Z = 2, W = 3 };
    }
}

//------------------------------------------------------------------------------
// classes
//------------------------------------------------------------------------------
#include "xmath_fvec4.h"
#include "xmath_fvec3.h"
#include "xmath_fvec2.h"

//------------------------------------------------------------------------------
// inline functions
//------------------------------------------------------------------------------
#include "implementation/xmath_fvec3_inline.h"
#include "implementation/xmath_fvec4_inline.h"
#include "implementation/xmath_fvec2_inline.h"

#endif
