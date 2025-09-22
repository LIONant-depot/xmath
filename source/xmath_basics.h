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

#include "xmath_strong_typing_numerics.h"
#include "xmath_functions.h"
#include "xmath_trigonometry.h"

#include "implementation/xmath_functions_inline.h"
#include "implementation/xmath_trigonometry_inline.h"

#endif