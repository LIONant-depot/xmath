#ifndef XMATH_FVEC3_UNIT_TEST_H
#define XMATH_FVEC3_UNIT_TEST_H
#pragma once

// Unit tests for fvec3_t using standard C++ asserts.
// No external libraries; plain asserts from <cassert>.
// Tests cover key functionality for both SIMD and CPU variants.
// Floating-point comparisons use a small epsilon.

#include <cassert>         // For assert
#include <cmath>           // For std::abs, std::sqrt, etc.
#include <limits>          // For numeric_limits

namespace xmath::unit_test::_fvec3
{
    constexpr float EPSILON = std::numeric_limits<float>::epsilon() * 100.0f;  // Adjusted for typical float precision

    //------------------------------------------------------------------------------
    // Helper function for approximate equality between fvec3_t instances
    template <bool UseSIMD>
    bool ApproxEqual(const fvec3_t<UseSIMD>& a, const fvec3_t<UseSIMD>& b, float eps = EPSILON)
    {
        return std::abs(a.m_X - b.m_X) < eps &&
               std::abs(a.m_Y - b.m_Y) < eps &&
               std::abs(a.m_Z - b.m_Z) < eps;
    }

    //------------------------------------------------------------------------------

    bool ApproxEqual(const float a, const float b, float eps = EPSILON)
    {
        return std::abs(a - b) < eps;
    }

    //------------------------------------------------------------------------------
    // Templated test functions for both SIMD and CPU
    template <bool UseSIMD>
    void TestConstructors()
    {
        assert(fvec3(3.0f, 2.0f, 1.0f).m_X == 3.0f );
        assert(fvec3(3.0f, 2.0f, 1.0f).m_Y == 2.0f);
        assert(fvec3(3.0f, 2.0f, 1.0f).m_Z == 1.0f);

        // Default Constructor (uninitialized) - just check size and alignment
        fvec3_t<UseSIMD> v_default;
        assert(UseSIMD && sizeof(v_default) == 16 || UseSIMD==false && sizeof(v_default) == 4*3 );

        // Component Constructor
        fvec3_t<UseSIMD> v(1.0f, 2.0f, 3.0f);
        assert(v.m_X == 1.0f);
        assert(v.m_Y == 2.0f);
        assert(v.m_Z == 3.0f);

        // Single Value Constructor
        fvec3_t<UseSIMD> v_single(5.0f);
        assert(v_single.m_X == 5.0f);
        assert(v_single.m_Y == 5.0f);
        assert(v_single.m_Z == 5.0f);

        // Pitch Yaw Constructor (assuming radian is float or similar; test forward direction)
        fvec3_t<UseSIMD> v_pitch_yaw(radian(0.0f), radian(0.0f));
        assert(ApproxEqual(v_pitch_yaw, fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));

        // Conversion Between Variants (test SIMD to CPU or vice versa)
        fvec3_t<!UseSIMD> v_other(1.0f, 2.0f, 3.0f);
        fvec3_t<UseSIMD> v_convert(v_other);
        assert(ApproxEqual(v_convert, fvec3_t<UseSIMD>(1.0f, 2.0f, 3.0f)));
    }

    //------------------------------------------------------------------------------

    template <bool UseSIMD>
    void TestStaticConstants()
{
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromZero(),       fvec3_t<UseSIMD>(0.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromOne(),        fvec3_t<UseSIMD>(1.0f, 1.0f, 1.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromUp(),         fvec3_t<UseSIMD>(0.0f, 1.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromDown(),       fvec3_t<UseSIMD>(0.0f, -1.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromLeft(),       fvec3_t<UseSIMD>(-1.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromRight(),      fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromForward(),    fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromBack(),       fvec3_t<UseSIMD>(0.0f, 0.0f, -1.0f)));
    }

    //------------------------------------------------------------------------------

    template <bool UseSIMD>
    void TestOperators()
    {
        fvec3_t<UseSIMD> a(1.0f, 2.0f, 3.0f);
        fvec3_t<UseSIMD> b(4.0f, 5.0f, 6.0f);

        // Addition
        fvec3_t<UseSIMD> c_add = a + b;
        assert(ApproxEqual(c_add, fvec3_t<UseSIMD>(5.0f, 7.0f, 9.0f)));

        // Subtraction
        fvec3_t<UseSIMD> c_sub = a - b;
        assert(ApproxEqual(c_sub, fvec3_t<UseSIMD>(-3.0f, -3.0f, -3.0f)));

        // Scalar Multiplication
        fvec3_t<UseSIMD> c_mul = a * 2.0f;
        assert(ApproxEqual(c_mul, fvec3_t<UseSIMD>(2.0f, 4.0f, 6.0f)));
        fvec3_t<UseSIMD> d_mul = 3.0f * a;
        assert(ApproxEqual(d_mul, fvec3_t<UseSIMD>(3.0f, 6.0f, 9.0f)));

        // Scalar Division
        fvec3_t<UseSIMD> c_div = a / 2.0f;
        assert(ApproxEqual(c_div, fvec3_t<UseSIMD>(0.5f, 1.0f, 1.5f)));

        // Compound Assignment
        fvec3_t<UseSIMD> c_comp = a;
        c_comp += b;
        assert(ApproxEqual(c_comp, fvec3_t<UseSIMD>(5.0f, 7.0f, 9.0f)));
        c_comp -= b;
        assert(ApproxEqual(c_comp, a));
        c_comp *= 2.0f;
        assert(ApproxEqual(c_comp, fvec3_t<UseSIMD>(2.0f, 4.0f, 6.0f)));
        c_comp /= 2.0f;
        assert(ApproxEqual(c_comp, a));

        // Equality
        assert(a == a);
        assert(a != b);

        // Index Operator
        assert(a[0] == 1.0f);
        assert(a[1] == 2.0f);
        assert(a[2] == 3.0f);
        fvec3_t<UseSIMD> idx_test = a;
        idx_test[0] = 10.0f;
        assert(idx_test[0] == 10.0f);
    }

    //------------------------------------------------------------------------------

    template <bool UseSIMD>
    void TestVectorOperations()
    {
        fvec3_t<UseSIMD> a(1.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> b(0.0f, 1.0f, 0.0f);

        // Dot Product
        {
            using TestCase = std::tuple<fvec3_t<UseSIMD>, fvec3_t<UseSIMD>, float>;
            std::vector<TestCase> testCases = {
                { {1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, 0.0f },   // Orthogonal
                { {1.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, 1.0f },   // Parallel
                { {1.0f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, -1.0f }, // Anti-parallel
                { {1.0f, 2.0f, 3.0f}, {4.0f, -5.0f, 6.0f}, 12.0f }, // Arbitrary
                { {0.0f, 0.0f, 0.0f}, {1.0f, 2.0f, 3.0f}, 0.0f },   // fromZero vector
                { {1.0f, 2.0f, 3.0f}, {1.0f, 2.0f, 3.0f}, 14.0f },  // Self dot
                { {-1.0f, -2.0f, -3.0f}, {1.0f, 2.0f, 3.0f}, -14.0f } // Negative components
            };

            for (const auto& [v1, v2, expected] : testCases)
            {
                float resultStatic = fvec3_t<UseSIMD>::Dot(v1, v2);
                float resultMember = v1.Dot(v2);

                assert(resultStatic == expected);
                assert(resultMember == expected);
            }
        }

        // Cross Product
        fvec3_t<UseSIMD> cross_ab = fvec3_t<UseSIMD>::Cross(a, b);
        assert(ApproxEqual(cross_ab, fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));
        assert(ApproxEqual(a.Cross(b), fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));

        // Length and Normalize
        fvec3_t<UseSIMD> v_len(3.0f, 4.0f, 0.0f);
        assert(v_len.Length() == 5.0f);
        assert(v_len.LengthSq() == 25.0f);
        fvec3_t<UseSIMD> v_norm = v_len.NormalizeCopy();
        assert(ApproxEqual(v_norm, fvec3_t<UseSIMD>(0.6f, 0.8f, 0.0f)));
        v_len.Normalize();
        assert(ApproxEqual(v_len, fvec3_t<UseSIMD>(0.6f, 0.8f, 0.0f)));

        // Safe Normalize (test zero vector)
        fvec3_t<UseSIMD> zero = fvec3_t<UseSIMD>::fromZero();
        fvec3_t<UseSIMD> safe_norm = zero.NormalizeSafeCopy();
        assert(ApproxEqual(safe_norm.Length(), 1.0f));  // Assuming safe normalize returns a normalized (unit length) vector even for zero input

        // Lerp
        fvec3_t<UseSIMD> lerp_ab = fvec3_t<UseSIMD>::Lerp(a, b, 0.5f);
        assert(ApproxEqual(lerp_ab, fvec3_t<UseSIMD>(0.5f, 0.5f, 0.0f)));
        assert(ApproxEqual(a.Lerp(b, 0.5f), fvec3_t<UseSIMD>(0.5f, 0.5f, 0.0f)));

        // Distance
        assert(fvec3_t<UseSIMD>::Distance(a, b) == std::sqrt(2.0f));
        assert(a.Distance(b) == std::sqrt(2.0f));

        // Min/Max
        assert(ApproxEqual(fvec3_t<UseSIMD>::Min(a, b), fvec3_t<UseSIMD>(0.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::Max(a, b), fvec3_t<UseSIMD>(1.0f, 1.0f, 0.0f)));
    }

    //------------------------------------------------------------------------------

    template <bool UseSIMD>
    void TestAdvancedOperations()
    {
        fvec3_t<UseSIMD> v(1.0f, 2.0f, 3.0f);

        // Abs
        fvec3_t<UseSIMD> v_abs = (-v).AbsCopy();
        assert(ApproxEqual(v_abs, v));

        // OneOver
        fvec3_t<UseSIMD> v_one_over = v.OneOverCopy();
        assert(ApproxEqual(v_one_over, fvec3_t<UseSIMD>(1.0f / 1.0f, 1.0f / 2.0f, 1.0f / 3.0f), 0.01f));

        // Reflection (simple case)
        fvec3_t<UseSIMD> incident(1.0f, -1.0f, 0.0f);
        fvec3_t<UseSIMD> normal(0.0f, 1.0f, 0.0f);
        fvec3_t<UseSIMD> reflected = incident.Reflection(normal);
        assert(ApproxEqual(reflected, fvec3_t<UseSIMD>(1.0f, 1.0f, 0.0f)));

        // isFinite
        assert(v.isFinite());
        fvec3_t<UseSIMD> inf_v(1.0f, std::numeric_limits<float>::infinity(), 3.0f);
        assert(!inf_v.isFinite());

        // AngleBetween
        fvec3_t<UseSIMD> a(1.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> b(0.0f, 1.0f, 0.0f);
        assert(xmath::Abs(a.AngleBetween(b).m_Value - (xmath::pi_v.m_Value / 2.0f)) < EPSILON);
    }

    //------------------------------------------------------------------------------
    // Run all tests for a specific variant
    template <bool UseSIMD>
    void RunVariantTests()
    {
        TestConstructors<UseSIMD>();
        TestStaticConstants<UseSIMD>();
        TestOperators<UseSIMD>();
        TestVectorOperations<UseSIMD>();
        TestAdvancedOperations<UseSIMD>();
    }

    //------------------------------------------------------------------------------
    // Main test runner - call this to run all tests
    void RunTests()
    {
        RunVariantTests<true>();   // SIMD
        RunVariantTests<false>();  // CPU
    }

}

#endif