#ifndef _XMATH_FQUAT_UNITTEST_H
#define _XMATH_FQUAT_UNITTEST_H
#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <array>
#include <string>
#include <utility>
#include <sstream>

namespace xmath::unit_test::_fquat
    {
    // Tolerance for floating point comparisons
    constexpr float TOL = 1e-5f;

    // Helper function to check if two floats are approximately equal
    bool approx_equal(float a, float b, float tol = TOL) {
        return std::fabs(a - b) < tol;
    }

    // Helper function to check if two fvec3 are approximately equal
    bool vec3_approx_equal(const fvec3& a, const fvec3& b, float tol = TOL) {
        return approx_equal(a.m_X, b.m_X, tol) &&
            approx_equal(a.m_Y, b.m_Y, tol) &&
            approx_equal(a.m_Z, b.m_Z, tol);
    }

    // Helper function to check if two fquat_t<V> are approximately equal
    template <bool V>
    bool quat_approx_equal(const fquat_t<V>& a, const fquat_t<V>& b, float tol = TOL) {
        return approx_equal(a.m_X, b.m_X, tol) &&
            approx_equal(a.m_Y, b.m_Y, tol) &&
            approx_equal(a.m_Z, b.m_Z, tol) &&
            approx_equal(a.m_W, b.m_W, tol);
    }

    // Helper function to check if radian3 are approximately equal
    bool radian3_approx_equal(const radian3& a, const radian3& b, float tol = TOL) {
        return approx_equal(a.m_Pitch.m_Value, b.m_Pitch.m_Value, tol) &&
            approx_equal(a.m_Yaw.m_Value, b.m_Yaw.m_Value, tol) &&
            approx_equal(a.m_Roll.m_Value, b.m_Roll.m_Value, tol);
    }

    // Helper function to check if pair<fvec3, radian> are approximately equal
    bool axis_angle_approx_equal(const std::pair<fvec3, radian>& a, const std::pair<fvec3, radian>& b, float tol = TOL) {
        return vec3_approx_equal(a.first, b.first, tol) && approx_equal(a.second.m_Value, b.second.m_Value, tol);
    }

    // Test functions for SIMD (true) and non-SIMD (false)
    template <bool V>
    void test_constructors()
    {
        assert( fquat_t<V>::fromIdentity() == fquat_t<V>(0,0,0,1) );

        // Default constructor - uninitialized, but check for zero in practice
        fquat_t<V> q_default;
        // No assert as uninitialized, but for test, set to zero and check
        q_default = fquat_t<V>::fromZero();
        assert(quat_approx_equal(q_default, fquat_t<V>::fromZero()));

        // Component constructor
        fquat_t<V> q_comp(0.1f, 0.2f, 0.3f, 0.4f);
        assert(approx_equal(q_comp.m_X, 0.1f));
        assert(approx_equal(q_comp.m_Y, 0.2f));
        assert(approx_equal(q_comp.m_Z, 0.3f));
        assert(approx_equal(q_comp.m_W, 0.4f));

        // Axis-angle constructor
        fvec3 axis(0.0f, 1.0f, 0.0f);
        radian angle(xmath::pi_v / 2.0f); // 90 degrees
        fquat_t<V> q_axis_angle(axis, angle);
        assert(approx_equal(q_axis_angle.m_X, 0.0f));
        assert(approx_equal(q_axis_angle.m_Y, 0.70710678118f, TOL));
        assert(approx_equal(q_axis_angle.m_Z, 0.0f));
        assert(approx_equal(q_axis_angle.m_W, 0.70710678118f, TOL));

        // Euler constructor
        radian3 euler(radian(0.0f), radian(xmath::pi_v / 2.0f), radian(0.0f));
        fquat_t<V> q_euler(euler);
        assert(quat_approx_equal(q_euler, fquat_t<V>(0.0f, 0.70710678118f, 0.0f, 0.70710678118f), TOL));

        // FromToRotation constructor with default up
        fvec3 from(1.0f, 0.0f, 0.0f);
        fvec3 to(0.0f, 1.0f, 0.0f);
        fquat_t<V> q_from_to(from, to, fvec3::fromUp() );
        assert(quat_approx_equal(q_from_to, fquat_t<V>(0.0f, 0.0f, 0.70710678118f, 0.70710678118f), TOL));

        // FromToRotation with custom up
        fvec3 up_custom(0.0f, 0.0f, 1.0f);
        fquat_t<V> q_from_to_up(from, to, up_custom);
        assert(quat_approx_equal(q_from_to_up, fquat_t<V>(0.0f, 0.0f, 0.70710678118f, 0.70710678118f), TOL));

        // LookRotation constructor
        fvec3 forward(0.0f, 0.0f, 1.0f);
        fvec3 up(0.0f, 1.0f, 0.0f);
        fquat_t<V> q_look(forward, up);
        assert(quat_approx_equal(q_look, fquat_t<V>::fromIdentity()));

        // SIMD register constructor (only for V=true)
        if constexpr (V) {
            floatx4 reg = _mm_set_ps(0.4f, 0.3f, 0.2f, 0.1f);
            fquat_t<true> q_reg(reg);
            assert(approx_equal(q_reg.m_X, 0.1f));
            assert(approx_equal(q_reg.m_Y, 0.2f));
            assert(approx_equal(q_reg.m_Z, 0.3f));
            assert(approx_equal(q_reg.m_W, 0.4f));
        }

        // Copy from other variant
        fquat_t<!V> q_other(0.1f, 0.2f, 0.3f, 0.4f);
        fquat_t<V> q_copy(q_other);
        assert(quat_approx_equal(q_copy, fquat_t<V>(0.1f, 0.2f, 0.3f, 0.4f)));

        // From double array
        std::array<double, 4> arr = { 0.1, 0.2, 0.3, 0.4 };
        fquat_t<V> q_arr(arr);
        assert(quat_approx_equal(q_arr, fquat_t<V>(0.1f, 0.2f, 0.3f, 0.4f)));
    }

    template <bool V>
    void test_assignment_conversion()
    {
        fquat_t<V> q(0.1f, 0.2f, 0.3f, 0.4f);

        // To double array
        std::array<double, 4> arr = q;
        assert(approx_equal(static_cast<float>(arr[0]), 0.1f));
        assert(approx_equal(static_cast<float>(arr[1]), 0.2f));
        assert(approx_equal(static_cast<float>(arr[2]), 0.3f));
        assert(approx_equal(static_cast<float>(arr[3]), 0.4f));

        // To string
        std::string str = q;
        // Check format "(0.100000, 0.200000, 0.300000, 0.400000)" - approximate match
        assert(str.find("0.1") != std::string::npos && str.find("0.2") != std::string::npos && str.find("0.3") != std::string::npos && str.find("0.4") != std::string::npos);

        // ToString
        std::string to_str = q.ToString();
        assert(to_str == str);

        // operator<<
        std::ostringstream os;
        os << q;
        assert(os.str() == str);
    }

    template <bool V>
    void test_static_properties()
    {
        // fromIdentity
        fquat_t<V> id = fquat_t<V>::fromIdentity();
        assert(quat_approx_equal(id, fquat_t<V>(0.0f, 0.0f, 0.0f, 1.0f)));

        // fromZero
        fquat_t<V> zero = fquat_t<V>::fromZero();
        assert(quat_approx_equal(zero, fquat_t<V>(0.0f, 0.0f, 0.0f, 0.0f)));
    }

    template <bool V>
    void test_static_methods()
    {
        fquat_t<V> a = fquat_t<V>::fromIdentity();
        fquat_t<V> b(0.0f, 0.70710678118f, 0.0f, 0.70710678118f);

        // Dot
        float dot = fquat_t<V>::Dot(a, b);
        assert(approx_equal(dot, 0.70710678118f));

        // Lerp
        fquat_t<V> lerp = fquat_t<V>::Lerp(a, b, 0.5f);
        assert(lerp.isNormalized());

        // Slerp
        fquat_t<V> slerp = fquat_t<V>::Slerp(a, b, 0.5f);
        assert(slerp.isNormalized());

        // SlerpAccurate
        fquat_t<V> slerp_acc = fquat_t<V>::SlerpAccurate(a, b, 0.5f);
        assert(slerp_acc.isNormalized());
        assert(quat_approx_equal(slerp, slerp_acc, 1e-4f));

        // AngleBetween
        radian ang = fquat_t<V>::AngleBetween(a, b);
        assert(approx_equal(ang.m_Value, xmath::pi_over2_v.m_Value, TOL));

        // FromToRotation
        fvec3 from(1.0f, 0.0f, 0.0f);
        fvec3 to(0.0f, 1.0f, 0.0f);
        fvec3 up(0.0f, 0.0f, 1.0f);
        fquat_t<V> from_to = fquat_t<V>::FromToRotation(from, to, up);
        assert(from_to.isNormalized());

        // LookRotation
        fvec3 forward(0.0f, 0.0f, 1.0f);
        fquat_t<V> look = fquat_t<V>::LookRotation(forward, fvec3::fromUp());
        assert(quat_approx_equal(look, fquat_t<V>::fromIdentity()));

        // RandomUnitQuaternion
        fquat_t<V> rand = fquat_t<V>::RandomUnitQuaternion();
        assert(rand.isNormalized());
    }

    template <bool V>
    void test_instance_basic()
    {
        fquat_t<V> q(0.5f, 0.5f, 0.5f, 0.5f);

        // Length
        float len = q.Length();
        assert(approx_equal(len, 1.0f));

        // LengthSq
        float len_sq = q.LengthSq();
        assert(approx_equal(len_sq, 1.0f));

        // NormalizeCopy
        fquat_t<V> norm_copy = q.NormalizeCopy();
        assert(quat_approx_equal(norm_copy, q)); // Already unit

        // Normalize
        fquat_t<V> norm = q;
        norm.Normalize();
        assert(quat_approx_equal(norm, q));

        // NormalizeSafeCopy
        fquat_t<V> zero = fquat_t<V>::fromZero();
        fquat_t<V> safe_copy = zero.NormalizeSafeCopy();
        assert(quat_approx_equal(safe_copy, fquat_t<V>::fromIdentity()));

        // NormalizeSafe
        fquat_t<V> safe = zero;
        safe.NormalizeSafe();
        assert(quat_approx_equal(safe, fquat_t<V>::fromIdentity()));

        // isFinite
        assert(q.isFinite());
        fquat_t<V> inf(std::numeric_limits<float>::infinity(), 0.0f, 0.0f, 0.0f);
        assert(!inf.isFinite());

        // isNormalized
        assert(q.isNormalized());

        // isNearlyIdentity
        assert(fquat_t<V>::fromIdentity().isNearlyIdentity());
        assert(!q.isNearlyIdentity());

        // Equals
        assert(q.Equals(q));
        assert(!q.Equals(fquat_t<V>::fromIdentity()));
    }

    template <bool V>
    void test_quat_specific()
    {
        fquat_t<V> q(0.0f, 0.70710678118f, 0.0f, 0.70710678118f); // 90 deg Y

        // ConjugateCopy
        fquat_t<V> conj_copy = q.ConjugateCopy();
        assert(approx_equal(conj_copy.m_Y, -0.70710678118f, TOL));
        assert(approx_equal(conj_copy.m_W, 0.70710678118f, TOL));

        // Conjugate
        fquat_t<V> conj = q;
        conj.Conjugate();
        assert(quat_approx_equal(conj, conj_copy));

        // InverseCopy
        fquat_t<V> inv_copy = q.InverseCopy();
        assert(quat_approx_equal(inv_copy, conj_copy)); // Unit quat

        // Inverse
        fquat_t<V> inv = q;
        inv.Inverse();
        assert(quat_approx_equal(inv, inv_copy));

        // Axis
        fvec3 axis = q.Axis();
        assert(vec3_approx_equal(axis, fvec3(0.0f, 1.0f, 0.0f)));

        // Angle
        radian ang = q.Angle();
        assert(approx_equal(ang.m_Value, xmath::pi_over2_v.m_Value, TOL));

        // ToEuler
        radian3 euler = q.ToEuler();
        assert(radian3_approx_equal(euler, radian3(radian(0.0f), radian(xmath::pi_v / 2.0f), radian(0.0f)), TOL));

        // Forward
        fvec3 fwd = q.Forward();
        assert(vec3_approx_equal(fwd, fvec3(1.0f, 0.0f, 0.0f), TOL));

        // Up
        fvec3 up_vec = q.Up();
        assert(vec3_approx_equal(up_vec, fvec3(0.0f, 1.0f, 0.0f)));

        // Right
        fvec3 right = q.Right();
        assert(vec3_approx_equal(right, fvec3(0.0f, 0.0f, -1.0f), TOL));

        // ToAxisAngle
        auto [axis_aa, ang_aa] = q.ToAxisAngle();
        assert(vec3_approx_equal(axis_aa, fvec3(0.0f, 1.0f, 0.0f)));
        assert(approx_equal(ang_aa.m_Value, xmath::pi_over2_v.m_Value, TOL));

        // Delta
        fquat_t<V> delta = q.Delta(fquat_t<V>::fromIdentity());
        assert(quat_approx_equal(delta, inv_copy));

        // LogCopy
        fquat_t<V> log_copy = q.LogCopy();
        assert(approx_equal(log_copy.m_Y, xmath::pi_over4_v.m_Value, TOL));
        assert(approx_equal(log_copy.m_W, 0.0f));

        // Log
        fquat_t<V> log = q;
        log.Log();
        assert(quat_approx_equal(log, log_copy));

        // ExpCopy
        fquat_t<V> exp_copy = log_copy.ExpCopy();
        assert(quat_approx_equal(exp_copy, q, TOL));

        // Exp
        fquat_t<V> exp = log_copy;
        exp.Exp();
        assert(quat_approx_equal(exp, exp_copy));
    }

    template <bool V>
    void test_rotation_operations()
{
        fquat_t<V> q = fquat_t<V>::fromIdentity();

        // setupRotationX
        fquat_t<V> rx = q;
        rx.setupRotationX(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rx, fquat_t<V>(0.70710678118f, 0.0f, 0.0f, 0.70710678118f), TOL));

        // setupRotationY
        fquat_t<V> ry = q;
        ry.setupRotationY(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(ry, fquat_t<V>(0.0f, 0.70710678118f, 0.0f, 0.70710678118f), TOL));

        // setupRotationZ
        fquat_t<V> rz = q;
        rz.setupRotationZ(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rz, fquat_t<V>(0.0f, 0.0f, 0.70710678118f, 0.70710678118f), TOL));

        // RotateXCopy
        fquat_t<V> rot_x_copy = q.RotateXCopy(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rot_x_copy, rx));

        // RotateYCopy
        fquat_t<V> rot_y_copy = q.RotateYCopy(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rot_y_copy, ry));

        // RotateZCopy
        fquat_t<V> rot_z_copy = q.RotateZCopy(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rot_z_copy, rz));

        // RotateX
        fquat_t<V> rot_x = q;
        rot_x.RotateX(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rot_x, rot_x_copy));

        // RotateY
        fquat_t<V> rot_y = q;
        rot_y.RotateY(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rot_y, rot_y_copy));

        // RotateZ
        fquat_t<V> rot_z = q;
        rot_z.RotateZ(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(rot_z, rot_z_copy));

        // PreRotateX
        fquat_t<V> pre_rot_x = q;
        pre_rot_x.PreRotateX(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(pre_rot_x, fquat_t<V>(0.70710678118f, 0.0f, 0.0f, 0.70710678118f), TOL));

        // PreRotateY
        fquat_t<V> pre_rot_y = q;
        pre_rot_y.PreRotateY(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(pre_rot_y, ry));

        // PreRotateZ
        fquat_t<V> pre_rot_z = q;
        pre_rot_z.PreRotateZ(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(pre_rot_z, rz));

        // PreRotateXCopy
        fquat_t<V> pre_rot_x_copy = q.PreRotateXCopy(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(pre_rot_x_copy, rx));

        // PreRotateYCopy
        fquat_t<V> pre_rot_y_copy = q.PreRotateYCopy(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(pre_rot_y_copy, ry));

        // PreRotateZCopy
        fquat_t<V> pre_rot_z_copy = q.PreRotateZCopy(radian(xmath::pi_v / 2.0f));
        assert(quat_approx_equal(pre_rot_z_copy, rz));

        // RotateTowardsCopy
        fquat_t<V> target = fquat_t<V>(0.0f, 0.0f, 0.0f, 1.0f);
        fquat_t<V> towards_copy = rx.RotateTowardsCopy(target, radian(xmath::pi_v / 4.0f));
        assert(towards_copy.isNormalized());

        // RotateTowards
        fquat_t<V> towards = rx;
        towards.RotateTowards(target, radian(xmath::pi_v / 4.0f));
        assert(quat_approx_equal(towards, towards_copy));

        // RotateVector
        fvec3 vec(1.0f, 0.0f, 0.0f);
        fvec3 rotated = rx.RotateVector(vec);
        assert(vec3_approx_equal(rotated, fvec3(1.0f, 0.0f, 0.0f), TOL)); // Adjust for actual rotation

        // In test_static_methods<V>
        fquat_t<V> a_tan(0.0f, 0.5f, 0.0f, 0.86602540378f); // Arbitrary tangent
        fquat_t<V> b_tan(0.0f, 0.0f, 0.5f, 0.86602540378f);
        fquat_t<V> a = fquat_t<V>::fromIdentity(); // a is the identity quaternion (0,0,0,1)
        fquat_t<V> b(0.0f, 0.70710678f, 0.0f, 0.70710678f); // b is a quaternion representing a 90-degree rotation around the Y-axis (approximately)
        fquat_t<V> squad = fquat_t<V>::Squad(a, a_tan, b, b_tan, 0.5f);
        assert(squad.isNormalized());

        // In test_instance_basic<V>
        fquat_t<V> q_from_to;
        q_from_to.setupFromToRotation(fvec3(1.0f, 0.0f, 0.0f), fvec3(0.0f, 1.0f, 0.0f), fvec3(0.0f, 0.0f, 1.0f));
        assert(quat_approx_equal(q_from_to, fquat_t<V>(0.0f, 0.0f, 0.70710678118f, 0.70710678118f), TOL));

        fquat_t<V> q_look;
        q_look.setupLookRotation(fvec3(0.0f, 0.0f, 1.0f), fvec3(0.0f, 1.0f, 0.0f));
        assert(quat_approx_equal(q_look, fquat_t<V>::fromIdentity()));

        assert(!q.isNearlyZero());
        fquat_t<V> zero = fquat_t<V>::fromZero();
        assert(zero.isNearlyZero());
    }

    template <bool V>
    void test_operators()
    {
        fquat_t<V> q1(0.1f, 0.2f, 0.3f, 0.4f);
        fquat_t<V> q2(0.4f, 0.3f, 0.2f, 0.1f);

        // +
        fquat_t<V> sum = q1 + q2;
        assert(quat_approx_equal(sum, fquat_t<V>(0.5f, 0.5f, 0.5f, 0.5f)));

        // -
        fquat_t<V> diff = q1 - q2;
        assert(quat_approx_equal(diff, fquat_t<V>(-0.3f, -0.1f, 0.1f, 0.3f)));

        // * scalar
        fquat_t<V> scaled = q1 * 2.0f;
        assert(quat_approx_equal(scaled, fquat_t<V>(0.2f, 0.4f, 0.6f, 0.8f)));

        // / scalar
        fquat_t<V> divided = q1 / 2.0f;
        assert(quat_approx_equal(divided, fquat_t<V>(0.05f, 0.1f, 0.15f, 0.2f)));

        // * quat
        fquat_t<V> prod = q1 * q2;
        assert(approx_equal(prod.m_Y, 0.24f, TOL));
        assert(approx_equal(prod.m_Z, 0.06f, TOL));
        assert(approx_equal(prod.m_W, -0.12f, TOL));

        // +=
        fquat_t<V> op_add = q1;
        op_add += q2;
        assert(quat_approx_equal(op_add, sum));

        // -=
        fquat_t<V> op_sub = q1;
        op_sub -= q2;
        assert(quat_approx_equal(op_sub, diff));

        // *= scalar
        fquat_t<V> op_mul_s = q1;
        op_mul_s *= 2.0f;
        assert(quat_approx_equal(op_mul_s, scaled));

        // /= scalar
        fquat_t<V> op_div_s = q1;
        op_div_s /= 2.0f;
        assert(quat_approx_equal(op_div_s, divided));

        // *= quat
        fquat_t<V> op_mul_q = q1;
        op_mul_q *= q2;
        assert(quat_approx_equal(op_mul_q, prod));

        // ==
        assert(q1 == q1);
        assert(!(q1 == q2));

        // !=
        assert(q1 != q2);
        assert(!(q1 != q1));

        // []
        assert(approx_equal(q1[0], 0.1f));
        assert(approx_equal(q1[1], 0.2f));
        assert(approx_equal(q1[2], 0.3f));
        assert(approx_equal(q1[3], 0.4f));

        // Mutable []
        fquat_t<V> mutable_q = q1;
        mutable_q[0] = 0.5f;
        assert(approx_equal(mutable_q.m_X, 0.5f));

        // Friend * scalar left
        fquat_t<V> scalar_left = 2.0f * q1;
        assert(quat_approx_equal(scalar_left, scaled));

        // Friend unary -
        fquat_t<V> neg = -q1;
        assert(quat_approx_equal(neg, fquat_t<V>(-0.1f, -0.2f, -0.3f, -0.4f)));

        // Friend * quat vec
        fvec3 vec(1.0f, 0.0f, 0.0f);
        fvec3 rotated = q1 * vec;
        assert(vec3_approx_equal(rotated, q1.RotateVector(vec)));
    }

    template <bool V>
    void test_quat_euler_stability()
    {
        // Helper to check if two quaternions represent the same rotation (q or -q)
        auto same_rotation = [](const fquat_t<V>& a, const fquat_t<V>& b, float tol = TOL) {
            return quat_approx_equal(a, b, tol) || quat_approx_equal(a, -b, tol);
            };

        // Test 1: fromIdentity quaternion
        fquat_t<V> q_id = fquat_t<V>::fromIdentity();
        radian3 euler_id = q_id.ToEuler();
        fquat_t<V> q_id_back = fquat_t<V>(euler_id);
        assert(same_rotation(q_id, q_id_back));

        // Test 2: 90 deg yaw (YXZ order)
        radian3 euler_yaw(radian(0.0f), radian(xmath::pi_v / 2.0f), radian(0.0f));
        fquat_t<V> q_yaw(euler_yaw);
        radian3 euler_yaw_back = q_yaw.ToEuler();
        fquat_t<V> q_yaw_recon(euler_yaw_back);
        assert(radian3_approx_equal(euler_yaw, euler_yaw_back));
        assert(same_rotation(q_yaw, q_yaw_recon));

        // Test 3: 90 deg pitch (gimbal lock case)
        radian3 euler_pitch(radian(xmath::pi_v / 2.0f), radian(0.0f), radian(0.0f));
        fquat_t<V> q_pitch(euler_pitch);
        radian3 euler_pitch_back = q_pitch.ToEuler();
        fquat_t<V> q_pitch_recon(euler_pitch_back);
        // Gimbal lock: roll/yaw ambiguous, but rotation should match
        assert(same_rotation(q_pitch, q_pitch_recon, 1e-3f)); // Increased tolerance for precision loss

        // Test 4: -90 deg pitch (gimbal lock)
        radian3 euler_neg_pitch(radian(-xmath::pi_v / 2.0f), radian(0.0f), radian(0.0f));
        fquat_t<V> q_neg_pitch(euler_neg_pitch);
        radian3 euler_neg_pitch_back = q_neg_pitch.ToEuler();
        fquat_t<V> q_neg_pitch_recon(euler_neg_pitch_back);
        assert(same_rotation(q_neg_pitch, q_neg_pitch_recon, 1e-3f)); // Increased tolerance for precision loss

        // Test 5: Arbitrary rotation
        radian3 euler_arb(radian(xmath::pi_v / 4.0f), radian(xmath::pi_v / 3.0f), radian(xmath::pi_v / 6.0f));
        fquat_t<V> q_arb(euler_arb);
        radian3 euler_arb_back = q_arb.ToEuler();
        fquat_t<V> q_arb_recon(euler_arb_back);
        assert(radian3_approx_equal(euler_arb, euler_arb_back, 1e-4f));
        assert(same_rotation(q_arb, q_arb_recon));

        // Test 6: Random quaternion (normalized)
        fquat_t<V> q_rand = fquat_t<V>::RandomUnitQuaternion();
        radian3 euler_rand = q_rand.ToEuler();
        fquat_t<V> q_rand_recon(euler_rand);
        assert(same_rotation(q_rand, q_rand_recon));

        // Test 7: Near gimbal lock with small roll
        radian3 euler_near_lock(radian(xmath::pi_over2_v.m_Value - 0.01f), radian(0.5f), radian(0.3f));
        fquat_t<V> q_near(euler_near_lock);
        radian3 euler_near_back = q_near.ToEuler();
        fquat_t<V> q_near_recon(euler_near_back);
        assert(same_rotation(q_near, q_near_recon, 1e-4f)); // Increased tolerance

        // Test 8: 180 deg rotation (axis ambiguous)
        fquat_t<V> q_180(0.0f, 0.0f, 1.0f, 0.0f); // 180 deg around Z
        radian3 euler_180 = q_180.ToEuler();
        fquat_t<V> q_180_recon(euler_180);
        assert(same_rotation(q_180, q_180_recon));
    }

    template <bool V>
    void run_all_tests()
    {
        test_constructors<V>();
        test_assignment_conversion<V>();
        test_static_properties<V>();
        test_static_methods<V>();
        test_instance_basic<V>();
        test_quat_specific<V>();
        test_rotation_operations<V>();
        test_operators<V>();
        test_quat_euler_stability<V>();
    }

    void RunTests(void)
    {
        run_all_tests<true>();
        run_all_tests<false>();
    }
}
#endif