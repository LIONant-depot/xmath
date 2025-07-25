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
        assert(fvec3(3.0f, 2.0f, 1.0f).m_X == 3.0f);
        assert(fvec3(3.0f, 2.0f, 1.0f).m_Y == 2.0f);
        assert(fvec3(3.0f, 2.0f, 1.0f).m_Z == 1.0f);

        // Default Constructor (uninitialized) - just check size and alignment
        fvec3_t<UseSIMD> v_default;
        assert(UseSIMD && sizeof(v_default) == 16 || UseSIMD == false && sizeof(v_default) == 4 * 3);

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

        // fvec2 + float Constructor
        fvec2 xy(4.0f, 5.0f);
        fvec3_t<UseSIMD> v_f2_float(xy, 6.0f);
        assert(v_f2_float.m_X == 4.0f && v_f2_float.m_Y == 5.0f && v_f2_float.m_Z == 6.0f);

        // float + fvec2 Constructor
        fvec3_t<UseSIMD> v_float_f2(7.0f, xy);
        assert(v_float_f2.m_X == 7.0f && v_float_f2.m_Y == 4.0f && v_float_f2.m_Z == 5.0f);

        // Span Constructor
        std::array<float, 3> span_arr = { 8.0f, 9.0f, 10.0f };
        std::span<float> span(span_arr);
        fvec3_t<UseSIMD> v_span(span);
        assert(v_span.m_X == 8.0f && v_span.m_Y == 9.0f && v_span.m_Z == 10.0f);

        // Array<double> Constructor
        std::array<double, 3> double_arr = { 11.0, 12.0, 13.0 };
        fvec3_t<UseSIMD> v_double(double_arr);
        assert(v_double.m_X == 11.0f && v_double.m_Y == 12.0f && v_double.m_Z == 13.0f);
    }

    //------------------------------------------------------------------------------

    template <bool UseSIMD>
    void TestStaticConstants()
    {
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromZero(), fvec3_t<UseSIMD>(0.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromOne(), fvec3_t<UseSIMD>(1.0f, 1.0f, 1.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromUp(), fvec3_t<UseSIMD>(0.0f, 1.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromDown(), fvec3_t<UseSIMD>(0.0f, -1.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromLeft(), fvec3_t<UseSIMD>(-1.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromRight(), fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromForward(), fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));
        assert(ApproxEqual(fvec3_t<UseSIMD>::fromBack(), fvec3_t<UseSIMD>(0.0f, 0.0f, -1.0f)));

        // Test fromRandomUnitVector - check if length is approximately 1
        fvec3_t<UseSIMD> random_unit = fvec3_t<UseSIMD>::fromRandomUnitVector();
        assert(ApproxEqual(random_unit.Length(), 1.0f, 1e-4f));
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

        // Component-wise Multiplication
        fvec3_t<UseSIMD> e_mul = a * b;
        assert(ApproxEqual(e_mul, fvec3_t<UseSIMD>(4.0f, 10.0f, 18.0f)));

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

                assert(ApproxEqual(resultStatic, expected));
                assert(ApproxEqual(resultMember, expected));
            }
        }

        // Cross Product
        fvec3_t<UseSIMD> cross_ab = fvec3_t<UseSIMD>::Cross(a, b);
        assert(ApproxEqual(cross_ab, fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));
        assert(ApproxEqual(a.Cross(b), fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));

        // Self cross should be zero
        assert(ApproxEqual(a.Cross(a), fvec3_t<UseSIMD>(0.0f, 0.0f, 0.0f)));

        // Anti-commutative
        assert(ApproxEqual(a.Cross(b), -b.Cross(a)));

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

        // Pitch and Yaw
        fvec3_t<UseSIMD> dir(1.0f, 0.0f, 1.0f);
        assert(ApproxEqual(dir.Pitch().m_Value, 0.0f));
        assert(ApproxEqual(dir.Yaw().m_Value, xmath::pi_over4_v.m_Value));

        auto [pitch, yaw] = dir.PitchYaw();
        assert(ApproxEqual(pitch.m_Value, 0.0f));
        assert(ApproxEqual(yaw.m_Value, xmath::pi_over4_v.m_Value));

        // RotationTowards
        fvec3_t<UseSIMD> from(1.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> to(0.0f, 1.0f, 0.0f);
        auto [axis, angle] = from.RotationTowards(to);
        assert(ApproxEqual(axis, fvec3_t<UseSIMD>(0.0f, 0.0f, 1.0f)));
        assert(ApproxEqual(angle.m_Value, xmath::pi_over2_v.m_Value));

        // SignedAngleBetween
        assert(ApproxEqual(from.SignedAngleBetween(to).m_Value, xmath::pi_over2_v.m_Value));

        // VectorToLineSegment, SquareDistToLineSeg, ClosestPointInLineSegment
        fvec3_t<UseSIMD> point(0.0f, 1.0f, 0.0f);
        fvec3_t<UseSIMD> start(0.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> end(2.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> vec_to_seg = point.VectorToLineSegment(start, end);
        assert(ApproxEqual(vec_to_seg, fvec3_t<UseSIMD>(0.0f, -1.0f, 0.0f)));
        assert(point.SquareDistToLineSeg(start, end) == 1.0f);
        assert(ApproxEqual(point.ClosestPointInLineSegment(start, end), fvec3_t<UseSIMD>(0.0f, 0.0f, 0.0f)));

        // ClosestPointToRectangle
        fvec3_t<UseSIMD> p0(0.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> e0(1.0f, 0.0f, 0.0f);
        fvec3_t<UseSIMD> e1(0.0f, 1.0f, 0.0f);
        fvec3_t<UseSIMD> rect_point(1.5f, 1.5f, 0.0f);
        fvec3_t<UseSIMD> closest;
        float dist_sq = rect_point.ClosestPointToRectangle(p0, e0, e1, closest);
        assert(ApproxEqual(closest, fvec3_t<UseSIMD>(1.0f, 1.0f, 0.0f)));
        assert(dist_sq == 0.5f);

        // RotateXCopy, RotateX
        fvec3_t<UseSIMD> rot_x = fvec3_t<UseSIMD>(0.0f, 1.0f, 0.0f).RotateXCopy(xmath::pi_over2_v);
        assert(ApproxEqual(rot_x, fvec3_t<UseSIMD>(0.0f, 0.0f, -1.0f)));
        rot_x.RotateX(radian(-xmath::pi_over2_v.m_Value));
        assert(ApproxEqual(rot_x, fvec3_t<UseSIMD>(0.0f, 1.0f, 0.0f)));

        // RotateYCopy, RotateY
        fvec3_t<UseSIMD> rot_y = fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f).RotateYCopy(xmath::pi_over2_v);
        assert(ApproxEqual(rot_y, fvec3_t<UseSIMD>(0.0f, 0.0f, -1.0f)));
        rot_y.RotateY(radian(-xmath::pi_over2_v.m_Value));
        assert(ApproxEqual(rot_y, fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f)));

        // RotateZCopy, RotateZ
        fvec3_t<UseSIMD> rot_z = fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f).RotateZCopy(xmath::pi_over2_v);
        assert(ApproxEqual(rot_z, fvec3_t<UseSIMD>(0.0f, 1.0f, 0.0f)));
        rot_z.RotateZ(radian(-xmath::pi_over2_v.m_Value));
        assert(ApproxEqual(rot_z, fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f)));

        // RotateCopy, Rotate
        radian3 euler(xmath::pi_over2_v, xmath::pi_over2_v, xmath::pi_over2_v);
        fvec3_t<UseSIMD> rot_euler = fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f).RotateCopy(euler);
        assert(ApproxEqual(rot_euler, fvec3_t<UseSIMD>(-1.0f, 0.0f, 0.0f), 1e-4f));

        // Test RotateCopy, Rotate
        euler = radian3(radian(1.57079637f), radian(1.57079637f), radian(1.57079637f));
        rot_euler = fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f).RotateCopy(euler);
        assert(ApproxEqual(rot_euler, fvec3_t<UseSIMD>(-1.0f, 0.0f, 0.0f), 1e-4f));
        fvec3_t<UseSIMD> back = rot_euler.RotateInverseCopy(euler);
        assert(ApproxEqual(back, fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f), 1e-4f));

        // Test in-place Rotate and RotateInverse
        fvec3_t<UseSIMD> rot_mut = fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f);
        rot_mut.Rotate(euler);
        assert(ApproxEqual(rot_mut, fvec3_t<UseSIMD>(-1.0f, 0.0f, 0.0f), 1e-4f));
        rot_mut.RotateInverse(euler);
        assert(ApproxEqual(rot_mut, fvec3_t<UseSIMD>(1.0f, 0.0f, 0.0f), 1e-4f));

        // GridSnap
        fvec3_t<UseSIMD> gs(1.2f, 2.6f, 3.4f);
        gs.GridSnap(1.f, 1.f, 1.f);
        assert(ApproxEqual(gs, fvec3_t<UseSIMD>(1.f, 3.f, 3.f)));

        // isRightHanded
        p0 = fvec3_t<UseSIMD>::fromZero();
        fvec3_t<UseSIMD> p1(1.f, 0.f, 0.f);
        fvec3_t<UseSIMD> p2(0.f, 1.f, 0.f);
        assert(p0.isRightHanded(p1, p2) == false); // Left-handed in 2D projection
        assert(p0.isRightHanded(p2, p1) == true); // Right-handed

        // ProjectCopy, Project
        fvec3_t<UseSIMD> proj_on = fvec3_t<UseSIMD>(3.f, 0.f, 0.f);
        fvec3_t<UseSIMD> proj_vec(1.f, 1.f, 1.f);
        fvec3_t<UseSIMD> projected = proj_vec.ProjectCopy(proj_on);
        assert(ApproxEqual(projected, fvec3_t<UseSIMD>(1.f, 0.f, 0.f)));
        proj_vec.Project(proj_on);
        assert(ApproxEqual(proj_vec, fvec3_t<UseSIMD>(1.f, 0.f, 0.f)));

        // Test Perpendicular
        proj_vec = fvec3_t<UseSIMD>(1.0f, 1.0f, 1.0f);
        normal = fvec3_t<UseSIMD>(0.0f, 1.0f, 0.0f);
        fvec3_t<UseSIMD> perp = proj_vec.Perpendicular(normal);
        assert(ApproxEqual(perp.Dot(normal), 0.0f, 1e-4f)); // Perpendicular to normal
        assert(ApproxEqual(perp.Dot(proj_vec), 0.0f, 1e-4f)); // Perpendicular to input vector
        assert(perp.Length() > 0.0f); // Non-zero length

        // ProjectOntoPlane
        fvec3_t<UseSIMD> plane_normal(0.f, 1.f, 0.f);
        fvec3_t<UseSIMD> proj_plane = proj_vec.ProjectOntoPlane(plane_normal);
        assert(ApproxEqual(proj_plane, fvec3_t<UseSIMD>(1.f, 0.f, 1.f)));

        // isNearlyZero
        assert(fvec3_t<UseSIMD>(0.f, 0.f, 0.f).isNearlyZero());
        assert(!fvec3_t<UseSIMD>(0.1f, 0.f, 0.f).isNearlyZero());

        // isNormalized
        fvec3_t<UseSIMD> unit(1.f, 0.f, 0.f);
        assert(unit.isNormalized());
        fvec3_t<UseSIMD> non_unit(2.f, 0.f, 0.f);
        assert(!non_unit.isNormalized());

        // MoveTowardsCopy, MoveTowards
        from = fvec3_t<UseSIMD>::fromZero();
        to = fvec3_t<UseSIMD>(0.f, 0.f, 10.f);
        fvec3_t<UseSIMD> moved = from.MoveTowardsCopy(to, 5.f);
        assert(ApproxEqual(moved, fvec3_t<UseSIMD>(0.f, 0.f, 5.f)));
        from.MoveTowards(to, 5.f);
        assert(ApproxEqual(from, fvec3_t<UseSIMD>(0.f, 0.f, 5.f)));

        // Test swizzles (a few examples)
        assert(ApproxEqual(v.x(), 1.0f));
        assert(ApproxEqual(v.xy().x(), fvec2(1.0f, 2.0f).x()));
        assert(ApproxEqual(v.xy().y(), fvec2(1.0f, 2.0f).y()));
        assert(ApproxEqual(v.xyz(), v));
        assert(ApproxEqual(v.xxxx().x(), fvec4(1.0f, 1.0f, 1.0f, 1.0f).x()));
        assert(ApproxEqual(v.xxxx().y(), fvec4(1.0f, 1.0f, 1.0f, 1.0f).y()));
        assert(ApproxEqual(v.xxxx().z(), fvec4(1.0f, 1.0f, 1.0f, 1.0f).z()));
        assert(ApproxEqual(v.xxxx().w(), fvec4(1.0f, 1.0f, 1.0f, 1.0f).w()));
    }

    //------------------------------------------------------------------------------

    template <bool UseSIMD>
    void TestVariantTests()
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
        TestVariantTests<true>();   // SIMD
        TestVariantTests<false>();  // CPU
    }

}

#endif




