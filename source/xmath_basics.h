#ifndef XMATH_BASICS_H
#define XMATH_BASICS_H
#pragma once

#include <concepts>
#include <numbers>
#include <type_traits>
#include <compare>
#include <smmintrin.h>
#include <cmath>
#include <bit>
#include <limits>
#include <cassert>
#include <span>

namespace xmath
{
    using   floatx4 = __m128;             // xmath own alias for simd data
}

// Portable construction / lane access for floatx4. MSVC's __m128 is a union (m128_f32[]); GCC/Clang's
// is a vector extension type that takes a plain brace list and supports operator[] directly.
#if defined(_MSC_VER) && !defined(__clang__)
    #define XMATH_FLOATX4(...)       ::xmath::floatx4{ .m128_f32{ __VA_ARGS__ } }
    #define XMATH_FLOATX4_LANE(V, I) (V).m128_f32[I]
#else
    #define XMATH_FLOATX4(...)       ::xmath::floatx4{ __VA_ARGS__ }
    #define XMATH_FLOATX4_LANE(V, I) (V)[I]
#endif

namespace xmath
{
}

#include "xmath_strong_typing_numerics.h"
#include "xmath_functions.h"
#include "xmath_trigonometry.h"

#include "implementation/xmath_functions_inline.h"
#include "implementation/xmath_trigonometry_inline.h"

#endif