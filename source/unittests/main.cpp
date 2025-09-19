#include "../../source/xmath.h"

#include "../../source/unittests/xmath_unittest.h"
#include "../../source/unittests/xmath_fvec4_unittest.h"
#include "../../source/unittests/xmath_fvec3_unittest.h"
#include "../../source/unittests/xmath_fvec2_unittest.h"
#include "../../source/unittests/xmath_radian3_unittest.h"
#include "../../source/unittests/xmath_fquat_unittest.h"
#include "../../source/unittests/xmath_fmat4_unittest.h"
#include "../../source/unittests/xmath_fmat3_unittest.h"
#include "../../source/unittests/xmath_fbbox_unittest.h"
#include "../../source/unittests/xmath_fplane_unittest.h"
#include "../../source/unittests/xmath_irect_unittest.h"
#include "../../source/unittests/xmath_ivec2_unittest.h"

int main()
{
    if (true) xmath::unit_test::_math::RunTests();
    if (true) xmath::unit_test::_fvec2::RunTests();
    if (true) xmath::unit_test::_fvec3::RunTests();
    if (true) xmath::unit_test::_fvec4::RunTests();
    if (true) xmath::unit_test::_radian3::RunTests();
    if (true) xmath::unit_test::_fquat::RunTests();
    if (true) xmath::unit_test::_fmat4::RunTests();
    if (true) xmath::unit_test::_fmat3::RunTests();
    if (true) xmath::unit_test::_fbbox::RunTests();
    if (true) xmath::unit_test::_fplane::RunTests();
    if (true) xmath::unit_test::_irect::RunTests();
    if (true) xmath::unit_test::_ivec2::RunTests();

    return 0;
}