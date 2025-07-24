#ifndef _XMATH_RADIAN3_UNIT_TEST_H
#define _XMATH_RADIAN3_UNIT_TEST_H
#pragma once

namespace  xmath::unit_test::_radian3
{
    constexpr float TEST_TOL = 0.001f;  // Tolerance for floating-point comparisons

    // Helper to check if two radians are approximately equal
    inline bool ApproxEqual(radian a, radian b, float tol = TEST_TOL)
    {
        return std::abs(a.m_Value - b.m_Value) < tol;
    }

    // Helper to check if two radian3 are approximately equal
    inline bool ApproxEqual(const radian3& a, const radian3& b, float tol = TEST_TOL)
    {
        return ApproxEqual(a.m_Pitch, b.m_Pitch, tol) &&
               ApproxEqual(a.m_Yaw,   b.m_Yaw,   tol) &&
               ApproxEqual(a.m_Roll,  b.m_Roll,  tol);
    }

    // Test constructors
    inline void TestConstructors()
    {
        // Default constructor
        radian3 r0;
        assert(ApproxEqual(r0.m_Pitch, radian{0.0f}));
        assert(ApproxEqual(r0.m_Yaw, radian{0.0f}));
        assert(ApproxEqual(r0.m_Roll, radian{0.0f}));

        // Single angle constructor
        radian3 r1(radian{1.0f});
        assert(ApproxEqual(r1.m_Pitch, radian{1.0f}));
        assert(ApproxEqual(r1.m_Yaw, radian{1.0f}));
        assert(ApproxEqual(r1.m_Roll, radian{1.0f}));

        // Three angles constructor
        radian3 r2(radian{0.5f}, radian{1.5f}, radian{2.5f});
        assert(ApproxEqual(r2.m_Pitch, radian{0.5f}));
        assert(ApproxEqual(r2.m_Yaw, radian{1.5f}));
        assert(ApproxEqual(r2.m_Roll, radian{2.5f}));
    }

    // Test fromZero
    inline void TestZero()
    {
        radian3 r(radian{1.0f}, radian{2.0f}, radian{3.0f});
        r.Zero();
        assert(ApproxEqual(r.m_Pitch, radian{0.0f}));
        assert(ApproxEqual(r.m_Yaw, radian{0.0f}));
        assert(ApproxEqual(r.m_Roll, radian{0.0f}));
    }

    // Test ModAngle and ModAngleCopy
    inline void TestModAngle()
    {
        radian3 r(radian{pi_v.m_Value * 3.0f}, radian{-pi_v.m_Value * 2.5f}, radian{pi_v.m_Value * 0.5f});
        radian3 copy = r.ModAngleCopy();
        assert(ApproxEqual(copy.m_Pitch, radian{pi_v.m_Value}));
        assert(ApproxEqual(copy.m_Yaw, radian{pi_v.m_Value * 1.5f}));
        assert(ApproxEqual(copy.m_Roll, radian{pi_v.m_Value * 0.5f}));

        r.ModAngle();
        assert(ApproxEqual(r, copy));
    }

    // Test ModAngle2 and ModAngle2Copy
    inline void TestModAngle2()
    {
        radian3 r(radian{pi_v.m_Value * 3.0f}, radian{-pi_v.m_Value * 2.5f}, radian{pi_v.m_Value * 0.5f});
        radian3 copy = r.ModAngle2Copy();
        assert(ApproxEqual(copy.m_Pitch, radian{-pi_v.m_Value}));
        assert(ApproxEqual(copy.m_Yaw,   radian{-pi_v.m_Value * 0.5f}));
        assert(ApproxEqual(copy.m_Roll,  radian{pi_v.m_Value * 0.5f}));

        r.ModAngle2();
        assert(ApproxEqual(r, copy));
    }

    // Test MinAngleDiff
    inline void TestMinAngleDiff()
    {
        radian3 r1(radian{0.0f}, radian{pi_v.m_Value * 1.8f}, radian{-pi_v.m_Value * 0.2f});
        radian3 r2(radian{pi_v.m_Value * 2.0f}, radian{0.0f}, radian{pi_v.m_Value * 1.9f});
        radian3 diff = r1.MinAngleDiff(r2);
        assert(ApproxEqual(diff.m_Pitch, radian{0.0f}));
        assert(ApproxEqual(diff.m_Yaw,  radian{-pi_v.m_Value * 0.2f}));
        assert(ApproxEqual(diff.m_Roll, radian{-pi_v.m_Value * 0.1f}));
    }

    // Test Difference
    inline void TestDifference()
    {
        radian3 r1(radian{1.0f}, radian{2.0f}, radian{3.0f});
        radian3 r2(radian{4.0f}, radian{5.0f}, radian{6.0f});
        radian diff = r1.Difference(r2);
        assert(ApproxEqual(diff, radian{std::sqrt(3.0f*3.0f + 3.0f*3.0f + 3.0f*3.0f)}));
    }

    // Test isInrange
    inline void TestIsInrange()
    {
        radian3 r(radian{0.5f}, radian{1.5f}, radian{2.5f});
        assert(r.isInrange(radian{0.0f}, radian{3.0f}));
        assert(!r.isInrange(radian{1.0f}, radian{2.0f}));
    }

    // Test isValid
    inline void TestIsValid()
    {
        radian3 r;
        assert(r.isValid());

        r.m_Pitch.m_Value = std::numeric_limits<float>::quiet_NaN();
        assert(!r.isValid());

        r.m_Pitch.m_Value = std::numeric_limits<float>::infinity();
        assert(!r.isValid());
    }

    // Test hasGimbalLock
    inline void TestHasGimbalLock()
    {
        radian3 r;
        assert(!r.hasGimbalLock());

        r.m_Pitch = pi_over2_v;
        assert(r.hasGimbalLock());

        r.m_Pitch = radian{pi_over2_v.m_Value + 0.1f};
        assert(!r.hasGimbalLock());
    }

    // Test Normalize and NormalizeCopy
    inline void TestNormalize()
    {
        radian3 r(radian{pi_v.m_Value * 3.0f}, radian{-pi_v.m_Value * 2.5f}, radian{pi_v.m_Value * 0.5f});
        radian3 copy = r.NormalizeCopy();
        assert(ApproxEqual(copy.m_Pitch, radian{-pi_v.m_Value}));
        assert(ApproxEqual(copy.m_Yaw,   radian{-pi_v.m_Value * 0.5f}));
        assert(ApproxEqual(copy.m_Roll,  radian{pi_v.m_Value * 0.5f}));

        r.Normalize();
        assert(ApproxEqual(r, copy));
    }

    // Test Clamp and ClampCopy
    inline void TestClamp()
    {
        radian3 r(radian{-1.0f}, radian{1.5f}, radian{3.0f});
        radian3 copy = r.ClampCopy(radian{0.0f}, radian{2.0f});
        assert(ApproxEqual(copy.m_Pitch, radian{0.0f}));
        assert(ApproxEqual(copy.m_Yaw, radian{1.5f}));
        assert(ApproxEqual(copy.m_Roll, radian{2.0f}));

        r.Clamp(radian{0.0f}, radian{2.0f});
        assert(ApproxEqual(r, copy));
    }

    // Test Lerp
    inline void TestLerp()
    {
        radian3 r1(radian{0.0f}, radian{pi_v.m_Value * 2.0f}, radian{-pi_v.m_Value});
        radian3 r2(radian{pi_v.m_Value}, radian{0.0f}, radian{pi_v.m_Value});
        radian3 lerped = r1.Lerp(0.5f, r2);
        assert(ApproxEqual(lerped.m_Pitch, radian{pi_v.m_Value * 0.5f}));
        assert(ApproxEqual(lerped.m_Yaw, radian{0.0f}));
        assert(ApproxEqual(lerped.m_Roll, radian{ pi_v.m_Value }));
    }

    // Test operators
    inline void TestOperators()
    {
        radian3 r1(radian{1.0f}, radian{2.0f}, radian{3.0f});
        radian3 r2(radian{4.0f}, radian{5.0f}, radian{6.0f});

        // ==
        assert(!(r1 == r2));
        assert(r1 == r1);

        // += -= *= /=
        radian3 r = r1;
        r += r2;
        assert(ApproxEqual(r, radian3(radian{5.0f}, radian{7.0f}, radian{9.0f})));
        r -= r2;
        assert(ApproxEqual(r, r1));
        r *= 2.0f;
        assert(ApproxEqual(r, radian3(radian{2.0f}, radian{4.0f}, radian{6.0f})));
        r /= 2.0f;
        assert(ApproxEqual(r, r1));

        // + - * /
        assert(ApproxEqual(r1 + r2, radian3(radian{5.0f}, radian{7.0f}, radian{9.0f})));
        assert(ApproxEqual(r1 - r2, radian3(radian{-3.0f}, radian{-3.0f}, radian{-3.0f})));
        assert(ApproxEqual(-r1, radian3(radian{-1.0f}, radian{-2.0f}, radian{-3.0f})));
        assert(ApproxEqual(r1 / 2.0f, radian3(radian{0.5f}, radian{1.0f}, radian{1.5f})));
        assert(ApproxEqual(r1 * 2.0f, radian3(radian{2.0f}, radian{4.0f}, radian{6.0f})));
        assert(ApproxEqual(2.0f * r1, radian3(radian{2.0f}, radian{4.0f}, radian{6.0f})));
    }

    // Run all tests
    inline void RunTests()
    {
        TestConstructors();
        TestZero();
        TestModAngle();
        TestModAngle2();
        TestMinAngleDiff();
        TestDifference();
        TestIsInrange();
        TestIsValid();
        TestHasGimbalLock();
        TestNormalize();
        TestClamp();
        TestLerp();
        TestOperators();
    }

} // namespace xmath::test

#endif