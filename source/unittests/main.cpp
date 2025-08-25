#include "../../source/xmath_flinear.h"

#include "../../source/unittests/xmath_unittest.h"
#include "../../source/unittests/xmath_fvec4_unittest.h"
#include "../../source/unittests/xmath_fvec3_unittest.h"
#include "../../source/unittests/xmath_fvec2_unittest.h"
#include "../../source/unittests/xmath_radian3_unittest.h"
#include "../../source/unittests/xmath_fquat_unittest.h"
#include "../../source/unittests/xmath_fmat4_unittest.h"
#include "../../source/unittests/xmath_fmat3_unittest.h"

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

    return 0;
}