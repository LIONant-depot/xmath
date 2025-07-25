#ifndef XMATH_FVEC4_UNITTEST_H
#define XMATH_FVEC4_UNITTEST_H

#include <cassert>
#include <cmath>
#include <limits>

namespace xmath::unit_test::_fvec4
{
    bool ApproxEqual(float a, float b, float eps = 1e-5f)
    {
        return std::abs(a - b) < eps;
    }

    bool ApproxEqual(const xmath::fvec4& a, const xmath::fvec4& b, float eps = 1e-5f)
    {
        return ApproxEqual(a.m_X, b.m_X, eps) &&
            ApproxEqual(a.m_Y, b.m_Y, eps) &&
            ApproxEqual(a.m_Z, b.m_Z, eps) &&
            ApproxEqual(a.m_W, b.m_W, eps);
    }

    void RunTests()
    {
        // Constructors
        xmath::fvec4 v1(1.f, 2.f, 3.f, 4.f);
        assert(v1.m_X == 1.f && v1.m_Y == 2.f && v1.m_Z == 3.f && v1.m_W == 4.f);

        xmath::fvec4 v_all(5.f);
        assert(v_all.m_X == 5.f && v_all.m_Y == 5.f && v_all.m_Z == 5.f && v_all.m_W == 5.f);

        // Assuming fvec3 is available with similar structure
        fvec3 f3(6.f, 7.f, 8.f);
        xmath::fvec4 v_from_f3(f3, 9.f);
        assert(v_from_f3.m_X == 6.f && v_from_f3.m_Y == 7.f && v_from_f3.m_Z == 8.f && v_from_f3.m_W == 9.f);

        xmath::fvec4 v_from_f3_front(10.f, f3);
        assert(v_from_f3_front.m_X == 10.f && v_from_f3_front.m_Y == 6.f && v_from_f3_front.m_Z == 7.f && v_from_f3_front.m_W == 8.f);

        floatx4 reg = _mm_set_ps(13.f, 12.f, 11.f, 10.f);
        fvec4 v_from_reg(reg);
        assert(v_from_reg.m_X == 10.f && v_from_reg.m_Y == 11.f && v_from_reg.m_Z == 12.f && v_from_reg.m_W == 13.f);

        // Assuming fvec2 is available
        fvec2 xy(14.f, 15.f);
        fvec2 zw(16.f, 17.f);
        fvec4 v_from_f2(xy, zw);
        assert(v_from_f2.m_X == 14.f && v_from_f2.m_Y == 15.f && v_from_f2.m_Z == 16.f && v_from_f2.m_W == 17.f);

        std::array<float, 4> arr = { 18.f, 19.f, 20.f, 21.f };
        std::span<float> span(arr);
        fvec4 v_from_span(span);
        assert(v_from_span.m_X == 18.f && v_from_span.m_Y == 19.f && v_from_span.m_Z == 20.f && v_from_span.m_W == 21.f);

        // Static properties
        fvec4 zero = fvec4::fromZero();
        assert(zero.m_X == 0.f && zero.m_Y == 0.f && zero.m_Z == 0.f && zero.m_W == 0.f);

        fvec4 one = fvec4::fromOne();
        assert(one.m_X == 1.f && one.m_Y == 1.f && one.m_Z == 1.f && one.m_W == 1.f);

        fvec4 unit_x = fvec4::fromUnitX();
        assert(unit_x.m_X == 1.f && unit_x.m_Y == 0.f && unit_x.m_Z == 0.f && unit_x.m_W == 0.f);

        fvec4 unit_y = fvec4::fromUnitY();
        assert(unit_y.m_X == 0.f && unit_y.m_Y == 1.f && unit_y.m_Z == 0.f && unit_y.m_W == 0.f);

        fvec4 unit_z = fvec4::fromUnitZ();
        assert(unit_z.m_X == 0.f && unit_z.m_Y == 0.f && unit_z.m_Z == 1.f && unit_z.m_W == 0.f);

        fvec4 unit_w = fvec4::fromUnitW();
        assert(unit_w.m_X == 0.f && unit_w.m_Y == 0.f && unit_w.m_Z == 0.f && unit_w.m_W == 1.f);

        fvec4 up = fvec4::fromUp();
        assert(up.m_X == 0.f && up.m_Y == 1.f && up.m_Z == 0.f && up.m_W == 0.f);

        fvec4 down = fvec4::fromDown();
        assert(down.m_X == 0.f && down.m_Y == -1.f && down.m_Z == 0.f && down.m_W == 0.f);

        fvec4 left = fvec4::fromLeft();
        assert(left.m_X == -1.f && left.m_Y == 0.f && left.m_Z == 0.f && left.m_W == 0.f);

        fvec4 right = fvec4::fromRight();
        assert(right.m_X == 1.f && right.m_Y == 0.f && right.m_Z == 0.f && right.m_W == 0.f);

        fvec4 forward = fvec4::fromForward();
        assert(forward.m_X == 0.f && forward.m_Y == 0.f && forward.m_Z == 1.f && forward.m_W == 0.f);

        fvec4 back = fvec4::fromBack();
        assert(back.m_X == 0.f && back.m_Y == 0.f && back.m_Z == -1.f && back.m_W == 0.f);

        // Test fromRandomUnitVector - check if length is approximately 1
        fvec4 random_unit = fvec4::fromRandomUnitVector();
        assert(ApproxEqual(random_unit.Length(), 1.f, 1e-4f));

        // Static methods
        float dot = fvec4::Dot(v1, v1);
        assert(dot == 30.f);

        fvec4 min_v = fvec4::Min(v1, fvec4(0.f, 3.f, 2.f, 5.f));
        assert(min_v.m_X == 0.f && min_v.m_Y == 2.f && min_v.m_Z == 2.f && min_v.m_W == 4.f);

        fvec4 max_v = fvec4::Max(v1, fvec4(0.f, 3.f, 2.f, 5.f));
        assert(max_v.m_X == 1.f && max_v.m_Y == 3.f && max_v.m_Z == 3.f && max_v.m_W == 5.f);

        fvec4 lerp_v = fvec4::Lerp(v1, zero, 0.5f);
        assert(lerp_v.m_X == 0.5f && lerp_v.m_Y == 1.f && lerp_v.m_Z == 1.5f && lerp_v.m_W == 2.f);

        float dist = fvec4::Distance(v1, zero);
        assert(ApproxEqual(dist, std::sqrt(30.f)));

        // Static methods as members
        assert(v1.Dot(v1) == 30.f);
        assert(ApproxEqual(v1.Min(fvec4(0.f, 3.f, 2.f, 5.f)), fvec4(0.f, 2.f, 2.f, 4.f)));
        assert(ApproxEqual(v1.Max(fvec4(0.f, 3.f, 2.f, 5.f)), fvec4(1.f, 3.f, 3.f, 5.f)));
        assert(ApproxEqual(v1.Lerp(zero, 0.5f), fvec4(0.5f, 1.f, 1.5f, 2.f)));
        assert(ApproxEqual(v1.Distance(zero), std::sqrt(30.f)));

        // Instance methods
        assert(ApproxEqual(v1.Length(), std::sqrt(30.f)));
        assert(v1.LengthSq() == 30.f);

        fvec4 norm_copy = v1.NormalizeCopy();
        assert(ApproxEqual(norm_copy.Length(), 1.f));
        assert(ApproxEqual(norm_copy, v1 / v1.Length()));

        fvec4 norm_mut = v1;
        norm_mut.Normalize();
        assert(ApproxEqual(norm_mut.Length(), 1.f));

        fvec4 zero_norm_safe_copy = zero.NormalizeSafeCopy();
        assert(zero_norm_safe_copy == zero);

        fvec4 zero_norm_safe_mut = zero;
        zero_norm_safe_mut.NormalizeSafe();
        assert(zero_norm_safe_mut == zero);

        fvec4 hom_copy = fvec4(2.f, 4.f, 6.f, 2.f).HomogeneousCopy();
        assert(hom_copy == fvec4(1.f, 2.f, 3.f, 1.f));

        fvec4 hom_mut = fvec4(2.f, 4.f, 6.f, 2.f);
        hom_mut.Homogenize();
        assert(hom_mut == fvec4(1.f, 2.f, 3.f, 1.f));

        assert(v1.isFinite());

        fvec4 inf_v(1.f, 2.f, std::numeric_limits<float>::infinity(), 4.f);
        assert(!inf_v.isFinite());

        assert(v1.isInRange(1.f, 4.f));
        assert(!v1.isInRange(2.f, 3.f));

        fvec4 one_over_copy = v1.OneOverCopy();
        assert(one_over_copy.m_X == 1.f && one_over_copy.m_Y == 0.5f && one_over_copy.m_Z == 1.f / 3.f && one_over_copy.m_W == 0.25f);

        fvec4 one_over_mut = v1;
        one_over_mut.OneOver();
        assert(ApproxEqual(one_over_mut, one_over_copy));

        fvec4 neg_v(-1.f, -2.f, -3.f, -4.f);
        fvec4 abs_copy = neg_v.AbsCopy();
        assert(abs_copy == fvec4(1.f, 2.f, 3.f, 4.f));

        fvec4 abs_mut = neg_v;
        abs_mut.Abs();
        assert(abs_mut == fvec4(1.f, 2.f, 3.f, 4.f));

        fvec4 incident(1.f, -1.f, 0.f, 0.f);
        fvec4 normal(0.f, 1.f, 0.f, 0.f);
        fvec4 refl = incident.Reflection(normal);
        assert(refl == fvec4(1.f, 1.f, 0.f, 0.f));

        assert(v1.DistanceSquare(zero) == 30.f);

        fvec4 _unit_x(1.f, 0.f, 0.f, 0.f);
        fvec4 _unit_y(0.f, 1.f, 0.f, 0.f);
        assert(ApproxEqual(_unit_x.AngleBetween(_unit_y).m_Value, xmath::pi_over2_v.m_Value));

        fvec4 gs(1.2f, 2.6f, 3.4f, 4.5f);
        gs.GridSnap(1.f, 1.f, 1.f, 1.f);
        assert(gs.m_X == 1.f && gs.m_Y == 3.f && gs.m_Z == 3.f && gs.m_W == 5.f);

        // Swizzle methods (test a few representatives)
        assert(v1.x() == 1.f);
        assert(v1.w() == 4.f);

        auto _xy = v1.xy();
        assert(_xy.m_X == 1.f && _xy.m_Y == 2.f);

        auto _zw = v1.zw();
        assert(_zw.m_X == 3.f && _zw.m_Y == 4.f);

        auto xxx = v1.xxx();
        assert(xxx.m_X == 1.f && xxx.m_Y == 1.f && xxx.m_Z == 1.f);

        auto yzw = v1.yzw();
        assert(yzw.m_X == 2.f && yzw.m_Y == 3.f && yzw.m_Z == 4.f);

        auto xxxx = v1.xxxx();
        assert(xxxx.m_X == 1.f && xxxx.m_Y == 1.f && xxxx.m_Z == 1.f && xxxx.m_W == 1.f);

        auto wzyx = v1.wzyx();
        assert(wzyx.m_X == 4.f && wzyx.m_Y == 3.f && wzyx.m_Z == 2.f && wzyx.m_W == 1.f);

        // Operator overloads
        assert(v1 + v1 == fvec4(2.f, 4.f, 6.f, 8.f));
        assert(v1 - v1 == zero);
        assert(v1 * 2.f == fvec4(2.f, 4.f, 6.f, 8.f));
        assert(v1 / 2.f == fvec4(0.5f, 1.f, 1.5f, 2.f));

        fvec4 mut_op = v1;
        mut_op += v1;
        assert(mut_op == fvec4(2.f, 4.f, 6.f, 8.f));

        mut_op -= v1;
        assert(mut_op == v1);

        mut_op *= 2.f;
        assert(mut_op == fvec4(2.f, 4.f, 6.f, 8.f));

        mut_op /= 2.f;
        assert(mut_op == v1);

        assert(v1 == fvec4(1.f, 2.f, 3.f, 4.f));
        assert(v1 != zero);

        assert(v1[0] == 1.f);
        assert(v1[1] == 2.f);
        assert(v1[2] == 3.f);
        assert(v1[3] == 4.f);

        // Friend operators
        assert(2.f * v1 == fvec4(2.f, 4.f, 6.f, 8.f));
        assert(-v1 == fvec4(-1.f, -2.f, -3.f, -4.f));

        // Component-wise math
        fvec4 pos_v(1.f, 2.f, 3.f, 4.f);
        neg_v = (-1.f, -2.f, -3.f, -4.f);
        fvec4 mixed_v(-1.f, 2.f, -3.f, 4.f);

        fvec4 abs_c = mixed_v.AbsCopy();
        assert(ApproxEqual(abs_c, pos_v));

        fvec4 abs_m = mixed_v;
        abs_m.Abs();
        assert(ApproxEqual(abs_m, pos_v));

        fvec4 one_over_c = pos_v.OneOverCopy();
        assert(ApproxEqual(one_over_c, fvec4(1.f, 0.5f, 1.f / 3.f, 0.25f)));

        fvec4 one_over_m = pos_v;
        one_over_m.OneOver();
        assert(ApproxEqual(one_over_m, fvec4(1.f, 0.5f, 1.f / 3.f, 0.25f)));

        fvec4 sqrt_c = fvec4(4.f, 9.f, 16.f, 25.f).SqrtCopy();
        assert(ApproxEqual(sqrt_c, fvec4(2.f, 3.f, 4.f, 5.f)));

        fvec4 sqrt_m = fvec4(4.f, 9.f, 16.f, 25.f);
        sqrt_m.Sqrt();
        assert(ApproxEqual(sqrt_m, fvec4(2.f, 3.f, 4.f, 5.f)));

        // Test SmoothStep
        fvec4 smoothstep = fvec4(0.f, 0.25f, 0.5f, 0.75f).SmoothStep(0.f, 1.f);
        assert(ApproxEqual(smoothstep.m_X, 0.f));
        assert(ApproxEqual(smoothstep.m_Y, 0.15625f, 1e-4f));
        assert(ApproxEqual(smoothstep.m_Z, 0.5f));
        assert(ApproxEqual(smoothstep.m_W, 0.84375f, 1e-4f));

        // Test InvSqrtCopy (accurate)
        fvec4 inv_sqrt_c = fvec4(4.f, 9.f, 16.f, 25.f).InvSqrtCopy();
        assert(ApproxEqual(inv_sqrt_c, fvec4(0.5f, 1.f / 3.f, 0.25f, 0.2f), 1e-5f));

        // Test InvSqrtCopyFast (fast, lower precision)
        fvec4 inv_sqrt_fast_c = fvec4(4.f, 9.f, 16.f, 25.f).InvSqrtFastCopy();
        assert(ApproxEqual(inv_sqrt_fast_c, fvec4(0.5f, 1.f / 3.f, 0.25f, 0.2f), 1e-3f));

        // Test in-place InvSqrt (accurate)
        fvec4 inv_sqrt_m = fvec4(4.f, 9.f, 16.f, 25.f);
        inv_sqrt_m.InvSqrt();

        assert(ApproxEqual(inv_sqrt_m, fvec4(0.5f, 1.f / 3.f, 0.25f, 0.2f), 1e-5f));
        inv_sqrt_m = fvec4(4.f, 9.f, 16.f, 25.f);

        inv_sqrt_m.InvSqrt();
        assert(ApproxEqual(inv_sqrt_m, fvec4(0.5f, 1.f / 3.f, 0.25f, 0.2f)));

        fvec4 sign_c = mixed_v.SignCopy();
        assert(ApproxEqual(sign_c, fvec4(-1.f, 1.f, -1.f, 1.f)));

        fvec4 sign_m = mixed_v;
        sign_m.Sign();
        assert(ApproxEqual(sign_m, fvec4(-1.f, 1.f, -1.f, 1.f)));

        fvec4 floor_c = fvec4(1.8f, 2.2f, -3.7f, 4.5f).FloorCopy();
        assert(ApproxEqual(floor_c, fvec4(1.f, 2.f, -4.f, 4.f)));

        fvec4 floor_m = fvec4(1.8f, 2.2f, -3.7f, 4.5f);
        floor_m.Floor();
        assert(ApproxEqual(floor_m, fvec4(1.f, 2.f, -4.f, 4.f)));

        fvec4 ceil_c = fvec4(1.2f, 2.8f, -3.1f, 4.5f).CeilCopy();
        assert(ApproxEqual(ceil_c, fvec4(2.f, 3.f, -3.f, 5.f)));

        fvec4 ceil_m = fvec4(1.2f, 2.8f, -3.1f, 4.5f);
        ceil_m.Ceil();
        assert(ApproxEqual(ceil_m, fvec4(2.f, 3.f, -3.f, 5.f)));

        fvec4 fract_c = fvec4(1.2f, 2.8f, -3.1f, 4.5f).FractCopy();
        assert(ApproxEqual(fract_c, fvec4(0.2f, 0.8f, -0.1f, 0.5f), 1e-4f));

        fvec4 fract_m = fvec4(1.2f, 2.8f, -3.1f, 4.5f);
        fract_m.Fract();
        assert(ApproxEqual(fract_m, fvec4(0.2f, 0.8f, -0.1f, 0.5f), 1e-4f));

        fvec4 round_c = fvec4(1.2f, 2.8f, -3.1f, 4.5f).RoundCopy();
        assert(ApproxEqual(round_c, fvec4(1.f, 3.f, -3.f, 4.f)));

        fvec4 round_m = fvec4(1.2f, 2.8f, -3.1f, 4.5f);
        round_m.Round();
        assert(ApproxEqual(round_m, fvec4(1.f, 3.f, -3.f, 4.f)));

        fvec4 trunc_c = fvec4(1.2f, 2.8f, -3.1f, 4.5f).TruncCopy();
        assert(ApproxEqual(trunc_c, fvec4(1.f, 2.f, -3.f, 4.f)));

        fvec4 trunc_m = fvec4(1.2f, 2.8f, -3.1f, 4.5f);
        trunc_m.Trunc();
        assert(ApproxEqual(trunc_m, fvec4(1.f, 2.f, -3.f, 4.f)));

        fvec4 mod_c = fvec4(5.f, 8.f, -3.f, 4.f).ModCopy(3.f);
        assert(ApproxEqual(mod_c, fvec4(2.f, 2.f, 0.f, 1.f)));

        fvec4 mod_m = fvec4(5.f, 8.f, -3.f, 4.f);
        mod_m.Mod(3.f);
        assert(ApproxEqual(mod_m, fvec4(2.f, 2.f, 0.f, 1.f)));

        fvec4 clamp_s_c = fvec4(-1.f, 0.5f, 2.f, 3.f).ClampCopy(0.f, 2.f);
        assert(ApproxEqual(clamp_s_c, fvec4(0.f, 0.5f, 2.f, 2.f)));

        fvec4 clamp_s_m = fvec4(-1.f, 0.5f, 2.f, 3.f);
        clamp_s_m.Clamp(0.f, 2.f);
        assert(ApproxEqual(clamp_s_m, fvec4(0.f, 0.5f, 2.f, 2.f)));

        fvec4 clamp_v_c = fvec4(-1.f, 0.5f, 2.f, 3.f).ClampCopy(fvec4(0.f, 0.f, 1.f, 2.f), fvec4(2.f, 2.f, 3.f, 4.f));
        assert(ApproxEqual(clamp_v_c, fvec4(0.f, 0.5f, 2.f, 3.f)));

        fvec4 clamp_v_m = fvec4(-1.f, 0.5f, 2.f, 3.f);
        clamp_v_m.Clamp(fvec4(0.f, 0.f, 1.f, 2.f), fvec4(2.f, 2.f, 3.f, 4.f));
        assert(ApproxEqual(clamp_v_m, fvec4(0.f, 0.5f, 2.f, 3.f)));

        fvec4 step = fvec4(-1.f, 0.5f, 1.f, 2.f).Step(1.f);
        assert(ApproxEqual(step, fvec4(0.f, 0.f, 1.f, 1.f)));

        // Test SmoothStep
        smoothstep = fvec4(0.f, 0.25f, 0.5f, 0.75f).SmoothStep(0.f, 1.f);
        assert(ApproxEqual(smoothstep.m_X, 0.f));
        assert(ApproxEqual(smoothstep.m_Y, 0.15625f, 1e-4f));
        assert(ApproxEqual(smoothstep.m_Z, 0.5f));
        assert(ApproxEqual(smoothstep.m_W, 0.84375f, 1e-4f));
        // Test with fvec4(0.f, 0.25f, 0.5f, 0.75f)
        // Expected: 0, 0.15625, 0.5, 0.84375

        fvec4 log_c = fvec4(1.f, xmath::e_v, xmath::e_v * xmath::e_v, 8.f).LogCopy();
        assert(ApproxEqual(log_c.m_X, 0.f));
        assert(ApproxEqual(log_c.m_Y, 1.f));
        assert(ApproxEqual(log_c.m_Z, 2.f));
        assert(ApproxEqual(log_c.m_W, std::log(8.f)));

        fvec4 log_m = fvec4(1.f, xmath::e_v, xmath::e_v * xmath::e_v, 8.f);
        log_m.Log();
        assert(ApproxEqual(log_m, log_c));

        fvec4 log2_c = fvec4(1.f, 2.f, 4.f, 8.f).Log2Copy();
        assert(ApproxEqual(log2_c.m_X, 0.f));
        assert(ApproxEqual(log2_c.m_Y, 1.f));
        assert(ApproxEqual(log2_c.m_Z, 2.f));
        assert(ApproxEqual(log2_c.m_W, 3.f));

        fvec4 log2_m = fvec4(1.f, 2.f, 4.f, 8.f);
        log2_m.Log2();
        assert(ApproxEqual(log2_m, log2_c));

        fvec4 pow_c = fvec4(2.f, 3.f, 4.f, 5.f).PowCopy(2.f);
        assert(ApproxEqual(pow_c, fvec4(4.f, 9.f, 16.f, 25.f)));

        fvec4 pow_m = fvec4(2.f, 3.f, 4.f, 5.f);
        pow_m.Pow(2.f);
        assert(ApproxEqual(pow_m, fvec4(4.f, 9.f, 16.f, 25.f)));

        fvec4 sin_c = fvec4(0.f, xmath::pi_v.m_Value / 2.f, xmath::pi_v.m_Value, 3 * xmath::pi_v.m_Value / 2.f).SinCopy();
        assert(ApproxEqual(sin_c.m_X, 0.f));
        assert(ApproxEqual(sin_c.m_Y, 1.f));
        assert(ApproxEqual(sin_c.m_Z, 0.f, 1e-4f));
        assert(ApproxEqual(sin_c.m_W, -1.f, 1e-4f));

        fvec4 sin_m = fvec4(0.f, xmath::pi_v.m_Value / 2.f, xmath::pi_v.m_Value, 3 * xmath::pi_v.m_Value / 2.f);
        sin_m.Sin();
        assert(ApproxEqual(sin_m, sin_c));

        fvec4 cos_c = fvec4(0.f, xmath::pi_v.m_Value / 2.f, xmath::pi_v.m_Value, 3 * xmath::pi_v.m_Value / 2.f).CosCopy();
        assert(ApproxEqual(cos_c.m_X, 1.f));
        assert(ApproxEqual(cos_c.m_Y, 0.f, 1e-4f));
        assert(ApproxEqual(cos_c.m_Z, -1.f, 1e-4f));
        assert(ApproxEqual(cos_c.m_W, 0.f, 1e-4f));

        fvec4 cos_m = fvec4(0.f, xmath::pi_v.m_Value / 2.f, xmath::pi_v.m_Value, 3 * xmath::pi_v.m_Value / 2.f);
        cos_m.Cos();
        assert(ApproxEqual(cos_m, cos_c));

        fvec4 tan_c = fvec4(0.f, xmath::pi_v.m_Value / 4.f, xmath::pi_v.m_Value / 2.f - 0.01f, -xmath::pi_v.m_Value / 4.f).TanCopy();
        assert(ApproxEqual(tan_c.m_X, 0.f));
        assert(ApproxEqual(tan_c.m_Y, 1.f));
        assert(ApproxEqual(tan_c.m_Z, xmath::Tan(xmath::pi_v / 2.f - radian(0.01f)), 1e-2f));
        assert(ApproxEqual(tan_c.m_W, -1.f));

        fvec4 tan_m = fvec4(0.f, xmath::pi_v.m_Value / 4.f, xmath::pi_v.m_Value / 2.f - 0.01f, -xmath::pi_v.m_Value / 4.f);
        tan_m.Tan();
        assert(ApproxEqual(tan_m, tan_c));

        fvec4 asin_c = fvec4(0.f, 0.5f, 1.f, -0.5f).AsinCopy();
        assert(ApproxEqual(asin_c.m_X, 0.f));
        assert(ApproxEqual(asin_c.m_Y, xmath::pi_v.m_Value / 6.f, 1e-4f));
        assert(ApproxEqual(asin_c.m_Z, xmath::pi_v.m_Value / 2.f, 1e-4f));
        assert(ApproxEqual(asin_c.m_W, -xmath::pi_v.m_Value / 6.f, 1e-4f));

        fvec4 asin_m = fvec4(0.f, 0.5f, 1.f, -0.5f);
        asin_m.Asin();
        assert(ApproxEqual(asin_m, asin_c));

        fvec4 acos_c = fvec4(1.f, 0.5f, 0.f, -0.5f).AcosCopy();
        assert(ApproxEqual(acos_c.m_X, 0.f));
        assert(ApproxEqual(acos_c.m_Y, xmath::pi_v.m_Value / 3.f, 1e-4f));
        assert(ApproxEqual(acos_c.m_Z, xmath::pi_v.m_Value / 2.f, 1e-4f));
        assert(ApproxEqual(acos_c.m_W, 2 * xmath::pi_v.m_Value / 3.f, 1e-4f));

        fvec4 acos_m = fvec4(1.f, 0.5f, 0.f, -0.5f);
        acos_m.Acos();
        assert(ApproxEqual(acos_m, acos_c));

        fvec4 atan_c = fvec4(0.f, 1.f, std::sqrt(3.f), -1.f).AtanCopy();
        assert(ApproxEqual(atan_c.m_X, 0.f));
        assert(ApproxEqual(atan_c.m_Y, xmath::pi_v.m_Value / 4.f, 1e-4f));
        assert(ApproxEqual(atan_c.m_Z, xmath::pi_v.m_Value / 3.f, 1e-4f));
        assert(ApproxEqual(atan_c.m_W, -xmath::pi_v.m_Value / 4.f, 1e-4f));

        fvec4 atan_m = fvec4(0.f, 1.f, std::sqrt(3.f), -1.f);
        atan_m.Atan();
        assert(ApproxEqual(atan_m, atan_c));

        fvec4 y(3.f, 4.f, 0.f, -3.f);
        fvec4 x(4.f, 3.f, 1.f, 4.f);
        fvec4 atan2_c = y.Atan2Copy(x);
        assert(ApproxEqual(atan2_c.m_X, xmath::Atan2(3.f, 4.f).m_Value));
        assert(ApproxEqual(atan2_c.m_Y, xmath::Atan2(4.f, 3.f).m_Value));
        assert(ApproxEqual(atan2_c.m_Z, xmath::Atan2(0.f, 1.f).m_Value));
        assert(ApproxEqual(atan2_c.m_W, xmath::Atan2(-3.f, 4.f).m_Value));

        fvec4 atan2_m = y;
        atan2_m.Atan2(x);
        assert(ApproxEqual(atan2_m, atan2_c));
    }
}

#endif
