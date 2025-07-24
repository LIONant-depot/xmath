#ifndef _XMATH_UNITTEST_H
#define _XMATH_UNITTEST_H
#pragma once

#include <cassert>
#include <cmath>

namespace xmath::unit_test::_math
{
    //------------------------------------------------------------------------------
    // RunUnitTests
    //------------------------------------------------------------------------------
    //
    // Runs unit tests for all xmath functions.
    //
    // Notes:
    //  Uses assert for failure; tests nominal/edge cases.
    //  Floating-point comparisons use FEqual with default tol.
    //  Covers trig, basic ops, angles, solvers; skips invalid inputs (assumed handled by func asserts).
    //
    inline void RunTests() noexcept
    {
        // Test angle functions
        {
            radian angle{ 3.0f * pi_v.m_Value };
            radian mod = ModAngle(angle);
            assert(FEqual(mod.m_Value, pi_v.m_Value));

            radian mod2 = ModAngle2(radian{ 1.5f * pi_v.m_Value });
            assert(FEqual(mod2.m_Value, -0.5f * pi_v.m_Value));

            radian diff = MinAngleDiff(radian{ pi_v.m_Value }, radian{ 0.0f });
            assert(FEqual(diff.m_Value, -pi_v.m_Value));

            radian lerped = LerpAngle(0.5f, radian{ 0.0f }, radian{ pi_v.m_Value });
            assert(FEqual(lerped.m_Value, 0.5f * pi_v.m_Value));
        }

        // Test SinCos
        {
            float s, c;
            SinCos(radian{ pi_over2_v.m_Value }, s, c);
            assert(FEqual(s, 1.0f) && FEqual(c, 0.0f));
        }

        // Test Sin, Cos, Tan
        {
            assert(FEqual(Sin(radian{ 0.0f }), 0.0f));
            assert(FEqual(Cos(radian{ pi_v.m_Value }), -1.0f));
            assert(FEqual(Tan(radian{ pi_over4_v.m_Value }), 1.0f));
        }

        // Test Asin, Acos
        {
            assert(FEqual(Asin(0.0f).m_Value, 0.0f));
            assert(FEqual(Acos(1.0f).m_Value, 0.0f));
        }

        // Test Exp, Pow
        {
            assert(FEqual(Exp(0.0f), 1.0f));
            assert(FEqual(Pow(2.0f, 3.0f), 8.0f));
        }

        // Test FMod, ModFX
        {
            assert(FEqual(FMod(5.5f, 2.0f), 1.5f));
            float intPart;
            assert(FEqual(ModFX(3.7f, intPart), 0.7f) && FEqual(intPart, 3.0f));
        }

        // Test Log, Log2, Log10
        {
            assert(FEqual(Log(1.0f), 0.0f));
            assert(FEqual(Log2(8.0f), 3.0f));
            assert(FEqual(Log10(100.0f), 2.0f));
        }

        // Test i2f, f2i
        {
            assert(FEqual(i2f(5.0f), 5.0f));
            assert(f2i(3.9f) == 3);
        }

        // Test FSel
        {
            assert(FSel(1.0f, 2.0f, 3.0f) == 2.0f);
            assert(FSel(-1.0f, 2.0f, 3.0f) == 3.0f);
        }

        // Test Sqr, Sqrt, InvSqrt
        {
            assert(FEqual(Sqr(3.0f), 9.0f));
            assert(FEqual(Sqrt(16.0f), 4.0f));
            assert(FEqual(InvSqrt(4.0f), 0.5f));
        }

        // Test Min, Max
        {
            assert(Min(1.0f, 2.0f) == 1.0f);
            assert(Max(1.0f, 2.0f) == 2.0f);
        }

        // Test FEqual, FLess, FGreater
        {
            assert(FEqual(1.0f, 1.0001f));
            assert(FLess(1.0f, 1.002f));
            assert(FGreater(1.002f, 1.0f));
        }

        // Test Sign
        {
            assert(Sign(0.0f));
            assert(!Sign(-1.0f));
        }

        // Test LRound, Round, Ceil, Floor
        {
            assert(LRound(2.6f) == 3);
            assert(FEqual(Round(5.3f, 2.0f), 6.0f));
            assert(FEqual(Ceil(2.1f), 3.0f));
            assert(FEqual(Floor(2.9f), 2.0f));
        }

        // Test isInRange, Range
        {
            assert(isInRange(5.0f, 1.0f, 10.0f));
            assert(FEqual(Range(15.0f, 1.0f, 10.0f), 10.0f));
        }

        // Test Abs, Lerp
        {
            assert(FEqual(Abs(-3.0f), 3.0f));
            assert(FEqual(Lerp(0.5f, 0.0f, 10.0f), 5.0f));
        }

        // Test isValid
        {
            assert(isValid(1.0f));
            assert(!isValid(std::numeric_limits<float>::infinity()));
        }

        // Test SolvedQuadraticRoots
        {
            float r1, r2;
            assert(SolvedQuadraticRoots(r1, r2, 1.0f, -3.0f, 2.0f));
            assert(FEqual(r1, 1.0f) || FEqual(r1, 2.0f));
            assert(!SolvedQuadraticRoots(r1, r2, 1.0f, 0.0f, 1.0f)); // No real roots
        }
    }
}
#endif