#ifndef _XMATH_FVEC2_UNITTEST_H
#define _XMATH_FVEC2_UNITTEST_H
#pragma once

#include <cassert>
#include <cmath>

namespace xmath::unit_test::_fvec2
{
    //------------------------------------------------------------------------------
    // Floating point comparison macro
    //------------------------------------------------------------------------------

#define ASSERT_FLOAT_EQ(expected, actual) assert(std::fabs((expected) - (actual)) < 1e-6f)
#define ASSERT_FLOAT2_EQ(expected, actual) \
    ASSERT_FLOAT_EQ((expected).m_X, (actual).m_X); \
    ASSERT_FLOAT_EQ((expected).m_Y, (actual).m_Y)

    //------------------------------------------------------------------------------
    // RunFloat2Tests
    //------------------------------------------------------------------------------
    //
    // Runs all unit tests for fvec2 class using assert.
    //
    // Notes:
    //  Throws assertion failure on error.
    //  Assumes fvec3 and fvec4 have similar component access for swizzle tests.
    //
    void RunTests(void)
    {
        //------------------------------------------------------------------------------
        // Test Constructors
        //------------------------------------------------------------------------------

        {
            fvec2 v(1.f, 2.f);
            ASSERT_FLOAT_EQ(1.f, v.m_X);
            ASSERT_FLOAT_EQ(2.f, v.m_Y);
        }

        {
            fvec2 v(3.f);
            ASSERT_FLOAT_EQ(3.f, v.m_X);
            ASSERT_FLOAT_EQ(3.f, v.m_Y);
        }

        {
            std::array<float, 2> arr = { 4.f, 5.f };
            std::span<float> sp(arr);
            fvec2 v(sp);
            ASSERT_FLOAT_EQ(4.f, v.m_X);
            ASSERT_FLOAT_EQ(5.f, v.m_Y);
        }

        //------------------------------------------------------------------------------
        // Test Static Properties
        //------------------------------------------------------------------------------

        {
            auto v = fvec2::fromZero();
            ASSERT_FLOAT_EQ(0.f, v.m_X);
            ASSERT_FLOAT_EQ(0.f, v.m_Y);
        }

        {
            auto v = fvec2::fromOne();
            ASSERT_FLOAT_EQ(1.f, v.m_X);
            ASSERT_FLOAT_EQ(1.f, v.m_Y);
        }

        {
            auto v = fvec2::fromUnitX();
            ASSERT_FLOAT_EQ(1.f, v.m_X);
            ASSERT_FLOAT_EQ(0.f, v.m_Y);
        }

        {
            auto v = fvec2::fromUnitY();
            ASSERT_FLOAT_EQ(0.f, v.m_X);
            ASSERT_FLOAT_EQ(1.f, v.m_Y);
        }

        {
            auto v = fvec2::fromUp();
            ASSERT_FLOAT_EQ(0.f, v.m_X);
            ASSERT_FLOAT_EQ(1.f, v.m_Y);
        }

        {
            auto v = fvec2::fromDown();
            ASSERT_FLOAT_EQ(0.f, v.m_X);
            ASSERT_FLOAT_EQ(-1.f, v.m_Y);
        }

        {
            auto v = fvec2::fromLeft();
            ASSERT_FLOAT_EQ(-1.f, v.m_X);
            ASSERT_FLOAT_EQ(0.f, v.m_Y);
        }

        {
            auto v = fvec2::fromRight();
            ASSERT_FLOAT_EQ(1.f, v.m_X);
            ASSERT_FLOAT_EQ(0.f, v.m_Y);
        }

        //------------------------------------------------------------------------------
        // Test Static Methods
        //------------------------------------------------------------------------------

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            float dot = fvec2::Dot(a, b);
            ASSERT_FLOAT_EQ(11.f, dot);
        }

        {
            fvec2 a(1.f, 4.f);
            fvec2 b(3.f, 2.f);
            auto min = fvec2::Min(a, b);
            ASSERT_FLOAT_EQ(1.f, min.m_X);
            ASSERT_FLOAT_EQ(2.f, min.m_Y);
        }

        {
            fvec2 a(1.f, 4.f);
            fvec2 b(3.f, 2.f);
            auto max = fvec2::Max(a, b);
            ASSERT_FLOAT_EQ(3.f, max.m_X);
            ASSERT_FLOAT_EQ(4.f, max.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            auto lerp = fvec2::Lerp(a, b, 0.5f);
            ASSERT_FLOAT_EQ(2.f, lerp.m_X);
            ASSERT_FLOAT_EQ(3.f, lerp.m_Y);
        }

        {
            fvec2 a(0.f, 0.f);
            fvec2 b(3.f, 4.f);
            float dist = fvec2::Distance(a, b);
            ASSERT_FLOAT_EQ(5.f, dist);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            float cross = fvec2::Cross(a, b);
            ASSERT_FLOAT_EQ(-2.f, cross);
        }

        //------------------------------------------------------------------------------
        // Test Static Methods as Members
        //------------------------------------------------------------------------------

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            float dot = a.Dot(b);
            ASSERT_FLOAT_EQ(11.f, dot);
        }

        {
            fvec2 a(1.f, 4.f);
            fvec2 b(3.f, 2.f);
            auto min = a.Min(b);
            ASSERT_FLOAT_EQ(1.f, min.m_X);
            ASSERT_FLOAT_EQ(2.f, min.m_Y);
        }

        {
            fvec2 a(1.f, 4.f);
            fvec2 b(3.f, 2.f);
            auto max = a.Max(b);
            ASSERT_FLOAT_EQ(3.f, max.m_X);
            ASSERT_FLOAT_EQ(4.f, max.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            auto lerp = a.Lerp(b, 0.5f);
            ASSERT_FLOAT_EQ(2.f, lerp.m_X);
            ASSERT_FLOAT_EQ(3.f, lerp.m_Y);
        }

        {
            fvec2 a(0.f, 0.f);
            fvec2 b(3.f, 4.f);
            float dist = a.Distance(b);
            ASSERT_FLOAT_EQ(5.f, dist);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            float cross = a.Cross(b);
            ASSERT_FLOAT_EQ(-2.f, cross);
        }

        //------------------------------------------------------------------------------
        // Test Instance Methods
        //------------------------------------------------------------------------------

        {
            fvec2 v(3.f, 4.f);
            float len = v.Length();
            ASSERT_FLOAT_EQ(5.f, len);
        }

        {
            fvec2 v(3.f, 4.f);
            float lenSq = v.LengthSq();
            ASSERT_FLOAT_EQ(25.f, lenSq);
        }

        {
            fvec2 v(3.f, 4.f);
            auto norm = v.NormalizeCopy();
            ASSERT_FLOAT_EQ(0.6f, norm.m_X);
            ASSERT_FLOAT_EQ(0.8f, norm.m_Y);
        }

        {
            fvec2 v(3.f, 4.f);
            v.Normalize();
            ASSERT_FLOAT_EQ(0.6f, v.m_X);
            ASSERT_FLOAT_EQ(0.8f, v.m_Y);
        }

        {
            fvec2 v(3.f, 4.f);
            auto norm = v.NormalizeSafeCopy();
            ASSERT_FLOAT_EQ(0.6f, norm.m_X);
            ASSERT_FLOAT_EQ(0.8f, norm.m_Y);
        }

        {
            fvec2 v(3.f, 4.f);
            v.NormalizeSafe();
            ASSERT_FLOAT_EQ(0.6f, v.m_X);
            ASSERT_FLOAT_EQ(0.8f, v.m_Y);
        }

        {
            fvec2 v(1.f, 2.f);
            bool finite = v.isFinite();
            assert(finite);
        }

        {
            fvec2 v(1.f, 2.f);
            bool inRange = v.isInRange(0.f, 3.f);
            assert(inRange);
        }

        {
            fvec2 v(2.f, 4.f);
            auto oneOver = v.OneOverCopy();
            ASSERT_FLOAT_EQ(0.5f, oneOver.m_X);
            ASSERT_FLOAT_EQ(0.25f, oneOver.m_Y);
        }

        {
            fvec2 v(2.f, 4.f);
            v.OneOver();
            ASSERT_FLOAT_EQ(0.5f, v.m_X);
            ASSERT_FLOAT_EQ(0.25f, v.m_Y);
        }

        {
            fvec2 v(-1.f, 2.f);
            auto abs = v.AbsCopy();
            ASSERT_FLOAT_EQ(1.f, abs.m_X);
            ASSERT_FLOAT_EQ(2.f, abs.m_Y);
        }

        {
            fvec2 v(-1.f, 2.f);
            v.Abs();
            ASSERT_FLOAT_EQ(1.f, v.m_X);
            ASSERT_FLOAT_EQ(2.f, v.m_Y);
        }

        {
            fvec2 v(1.f, 1.f);
            fvec2 normal(0.f, 1.f);
            auto refl = v.Reflection(normal);
            ASSERT_FLOAT_EQ(1.f, refl.m_X);
            ASSERT_FLOAT_EQ(-1.f, refl.m_Y);
        }

        {
            fvec2 a(0.f, 0.f);
            fvec2 b(3.f, 4.f);
            float distSq = a.DistanceSquare(b);
            ASSERT_FLOAT_EQ(25.f, distSq);
        }

        {
            fvec2 a(1.f, 0.f);
            fvec2 b(0.f, 1.f);
            radian angle = a.AngleBetween(b);
            ASSERT_FLOAT_EQ(1.570796f, angle.m_Value); // pi/2
        }

        {
            fvec2 v(1.3f, 2.7f);
            v.GridSnap(1.f, 1.f);
            ASSERT_FLOAT_EQ(1.f, v.m_X);
            ASSERT_FLOAT_EQ(3.f, v.m_Y);
        }

        {
            fvec2 v(1.f, 2.f);
            auto perp = v.Perp();
            ASSERT_FLOAT_EQ(-2.f, perp.m_X);
            ASSERT_FLOAT_EQ(1.f, perp.m_Y);
        }

        //------------------------------------------------------------------------------
        // Test Swizzle for float
        //------------------------------------------------------------------------------

        {
            fvec2 v(1.f, 2.f);
            ASSERT_FLOAT_EQ(1.f, v.x());
            ASSERT_FLOAT_EQ(2.f, v.y());
        }

        //------------------------------------------------------------------------------
        // Test Swizzle for fvec2
        //------------------------------------------------------------------------------

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
        }

        //------------------------------------------------------------------------------
        // Test Swizzle for fvec3 (assuming fvec3 has m_X, m_Y, m_Z)
        //------------------------------------------------------------------------------

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xxx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xxy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xyx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xyy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yxx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yxy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yyx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yyy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
        }

        //------------------------------------------------------------------------------
        // Test Swizzle for fvec4 (assuming fvec4 has m_X, m_Y, m_Z, m_W)
        //------------------------------------------------------------------------------

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xxxx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xxxy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xxyx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xxyy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xyxx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xyxy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xyyx();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.xyyy();
            ASSERT_FLOAT_EQ(1.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yxxx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yxxy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yxyx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yxyy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(1.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yyxx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yyxy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(1.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yyyx();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(1.f, sw.m_W);
        }

        {
            fvec2 v(1.f, 2.f);
            auto sw = v.yyyy();
            ASSERT_FLOAT_EQ(2.f, sw.m_X);
            ASSERT_FLOAT_EQ(2.f, sw.m_Y);
            ASSERT_FLOAT_EQ(2.f, sw.m_Z);
            ASSERT_FLOAT_EQ(2.f, sw.m_W);
        }

        //------------------------------------------------------------------------------
        // Test Operators
        //------------------------------------------------------------------------------

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            auto sum = a + b;
            ASSERT_FLOAT_EQ(4.f, sum.m_X);
            ASSERT_FLOAT_EQ(6.f, sum.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            auto diff = a - b;
            ASSERT_FLOAT_EQ(-2.f, diff.m_X);
            ASSERT_FLOAT_EQ(-2.f, diff.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            auto scaled = a * 3.f;
            ASSERT_FLOAT_EQ(3.f, scaled.m_X);
            ASSERT_FLOAT_EQ(6.f, scaled.m_Y);
        }

        {
            fvec2 a(3.f, 6.f);
            auto divided = a / 3.f;
            ASSERT_FLOAT_EQ(1.f, divided.m_X);
            ASSERT_FLOAT_EQ(2.f, divided.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            a += fvec2(3.f, 4.f);
            ASSERT_FLOAT_EQ(4.f, a.m_X);
            ASSERT_FLOAT_EQ(6.f, a.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            a -= fvec2(3.f, 4.f);
            ASSERT_FLOAT_EQ(-2.f, a.m_X);
            ASSERT_FLOAT_EQ(-2.f, a.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            a *= 3.f;
            ASSERT_FLOAT_EQ(3.f, a.m_X);
            ASSERT_FLOAT_EQ(6.f, a.m_Y);
        }

        {
            fvec2 a(3.f, 6.f);
            a /= 3.f;
            ASSERT_FLOAT_EQ(1.f, a.m_X);
            ASSERT_FLOAT_EQ(2.f, a.m_Y);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(1.f, 2.f);
            assert(a == b);
        }

        {
            fvec2 a(1.f, 2.f);
            fvec2 b(3.f, 4.f);
            assert(a != b);
        }

        {
            fvec2 v(1.f, 2.f);
            ASSERT_FLOAT_EQ(1.f, v[0]);
            ASSERT_FLOAT_EQ(2.f, v[1]);
        }

        {
            fvec2 v(1.f, 2.f);
            v[0] = 3.f;
            v[1] = 4.f;
            ASSERT_FLOAT_EQ(3.f, v.m_X);
            ASSERT_FLOAT_EQ(4.f, v.m_Y);
        }

        //------------------------------------------------------------------------------
        // Test Friend Operators
        //------------------------------------------------------------------------------

        {
            fvec2 v(1.f, 2.f);
            auto scaled = 3.f * v;
            ASSERT_FLOAT_EQ(3.f, scaled.m_X);
            ASSERT_FLOAT_EQ(6.f, scaled.m_Y);
        }

        {
            fvec2 v(1.f, 2.f);
            auto neg = -v;
            ASSERT_FLOAT_EQ(-1.f, neg.m_X);
            ASSERT_FLOAT_EQ(-2.f, neg.m_Y);
        }

        //------------------------------------------------------------------------------
        // Test Edge Cases
        //------------------------------------------------------------------------------

        {
            fvec2 zero(0.f, 0.f);
            auto norm = zero.NormalizeCopy();
            ASSERT_FLOAT_EQ(0.f, norm.m_X);
            ASSERT_FLOAT_EQ(0.f, norm.m_Y);
        }

        {
            fvec2 inf(std::numeric_limits<float>::infinity(), 0.f);
            assert(!inf.isFinite());
        }

        {
            fvec2 v(0.f, 0.f);
            radian angle = v.AngleBetween(fvec2(1.f, 0.f));
            ASSERT_FLOAT_EQ(0.f, angle.m_Value);
        }
    }

#undef ASSERT_FLOAT_EQ
#undef ASSERT_FLOAT2_EQ

} // namespace xmath::test

#endif