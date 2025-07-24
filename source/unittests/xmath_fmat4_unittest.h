namespace xmath::unit_test::_fmat4
{
    // Helper functions
    constexpr float EPSILON = 1e-5f;
    constexpr float SMALL_EPSILON = 1e-6f; // For more precise checks where needed

    bool approx_equal(float a, float b, float eps = EPSILON)
    {
        return std::fabs(a - b) < eps;
    }

    bool vec3_approx_equal(const fvec3& a, const fvec3& b, float eps = EPSILON)
    {
        return approx_equal(a.m_X, b.m_X, eps) && approx_equal(a.m_Y, b.m_Y, eps) && approx_equal(a.m_Z, b.m_Z, eps);
    }

    bool quat_approx_equal(const fquat& a, const fquat& b, float eps = EPSILON)
    {
        // Quaternions can be equivalent if negated, so check both
        if (approx_equal(a.m_W, b.m_W, eps) && vec3_approx_equal(fvec3(a.m_X, a.m_Y, a.m_Z), fvec3(b.m_X, b.m_Y, b.m_Z), eps)) return true;
        return approx_equal(a.m_W, -b.m_W, eps) && vec3_approx_equal(fvec3(a.m_X, a.m_Y, a.m_Z), fvec3(-b.m_X, -b.m_Y, -b.m_Z), eps);
    }

    bool mat4_approx_equal(const std::span<const float, 16>&& a, const std::span<const float, 16>&& b, float eps = EPSILON)
    {
        for (size_t i = 0; i < 16; ++i)
        {
            if (!approx_equal(a[i], b[i], eps)) return false;
        }
        return true;
    }

    bool vec4_approx_equal(const fvec4& a, const fvec4& b, float eps = EPSILON)
    {
        return approx_equal(a.m_X, b.m_X, eps) && approx_equal(a.m_Y, b.m_Y, eps) && approx_equal(a.m_Z, b.m_Z, eps) && approx_equal(a.m_W, b.m_W, eps);
    }

    template <bool V>
    void test_constructors()
    {
        // Diagonal
        fmat4_t<V> diag(2.0f);
        assert(approx_equal(diag(0, 0), 2.0f) && approx_equal(diag(1, 1), 2.0f) && approx_equal(diag(2, 2), 2.0f) && approx_equal(diag(3, 3), 2.0f));
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) if (i != j) assert(approx_equal(diag(i, j), 0.0f));

        // Array
        std::array<float, 16> arr = { 1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16 };
        fmat4_t<V> from_arr(arr);
        for (size_t i = 0; i < 16; ++i) assert(approx_equal(from_arr.m_Elements[i], static_cast<float>(i + 1)));

        // Span
        std::array<float, 16> span_arr = arr;
        fmat4_t<V> from_span(span_arr);

        assert(mat4_approx_equal(from_span, from_arr));

        // Quaternion
        fquat q_id = fquat::fromIdentity();
        fmat4_t<V> from_q(q_id);
        assert(mat4_approx_equal(from_q, fmat4_t<V>::fromIdentity()));

        fquat q_rot = fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v);
        fmat4_t<V> from_q_rot(q_rot);
        assert(quat_approx_equal(static_cast<fquat>(from_q_rot), q_rot));

        // Euler (assuming radian3 is fvec3 with radians)
        radian3 euler_zero = radian3::fromZero();
        fmat4_t<V> from_euler_zero(euler_zero);
        assert(mat4_approx_equal(from_euler_zero, fmat4_t<V>::fromIdentity()));

        radian3 euler_pi2(radian{ 0.f }, xmath::pi_over2_v, radian{ 0.f });
        fmat4_t<V> from_euler_pi2(euler_pi2);
        assert(vec3_approx_equal(from_euler_pi2.RotateVector(fvec3(1, 0, 0)), fvec3(0, 0, -1), SMALL_EPSILON));

        // TRS
        fvec3 trans(1, 2, 3), scale(4, 5, 6);
        fquat rot = fquat::fromIdentity();
        fmat4_t<V> from_trs(trans, rot, scale);
        assert(vec3_approx_equal(from_trs.ExtractPosition(), trans));
        assert(quat_approx_equal(from_trs.ExtractRotation(), rot));
        assert(vec3_approx_equal(from_trs.ExtractScale(), scale));

        // Cross-SIMD
        fmat4_t<!V> other = fmat4_t<!V>(trans, rot, scale);
        fmat4_t<V> from_other(other);
        assert(mat4_approx_equal(from_other, from_trs));
    }

    template <bool V>
    void test_static_constructors() {
        // fromIdentity
        auto id = fmat4_t<V>::fromIdentity();
        assert(approx_equal(id(0, 0), 1.0f) && approx_equal(id(1, 1), 1.0f) && approx_equal(id(2, 2), 1.0f) && approx_equal(id(3, 3), 1.0f));
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) if (i != j) assert(approx_equal(id(i, j), 0.0f));

        // fromZero
        auto zero = fmat4_t<V>::fromZero();
        for (float e : zero.m_Elements) assert(approx_equal(e, 0.0f));

        // fromTranslation
        fvec3 t(1, 2, 3);
        auto trans_mat = fmat4_t<V>::fromTranslation(t);
        assert(vec3_approx_equal(trans_mat.ExtractPosition(), t));
        assert(mat4_approx_equal(trans_mat.ClearTranslation(), fmat4_t<V>::fromIdentity()));

        // fromRotation (quat)
        auto rot_q = fmat4_t<V>::fromRotation(fquat::fromIdentity());
        assert(mat4_approx_equal(rot_q, id));

        fquat q_rot = fquat::fromAxisAngle(fvec3(1, 0, 0), xmath::pi_v);
        auto rot_q_pi = fmat4_t<V>::fromRotation(q_rot);
        assert(vec3_approx_equal(rot_q_pi.RotateVector(fvec3(0, 1, 0)), fvec3(0, -1, 0), SMALL_EPSILON));

        // fromRotation (axis-angle)
        fvec3 axis(0, 1, 0);
        radian angle = xmath::pi_over2_v;
        auto rot_aa = fmat4_t<V>::fromRotation(axis, angle);
        assert(vec3_approx_equal(rot_aa.RotateVector(fvec3(1, 0, 0)), fvec3(0, 0, -1), SMALL_EPSILON));

        // Edge case: zero angle
        auto rot_aa_zero = fmat4_t<V>::fromRotation(axis, radian{ 0.f });
        assert(mat4_approx_equal(rot_aa_zero, id));

        // fromScale
        fvec3 s(2, 3, 4);
        auto scale_mat = fmat4_t<V>::fromScale(s);
        assert(vec3_approx_equal(scale_mat.ExtractScale(), s));

        // fromPerspective (fov)
        radian fov{ 1.0472f }; // 60 deg
        float aspect = 1.333f, near_p = 0.1f, far_p = 100.0f;
        auto pers_fov = fmat4_t<V>::fromPerspective(fov, aspect, near_p, far_p);
        assert(approx_equal(pers_fov(0, 0), 1.0f / (aspect * std::tan(fov.m_Value / 2.0f)), SMALL_EPSILON));
        assert(approx_equal(pers_fov(1, 1), 1.0f / std::tan(fov.m_Value / 2.0f), SMALL_EPSILON));
//assert(approx_equal(pers_fov(2, 2), -(far_p + near_p) / (far_p - near_p), SMALL_EPSILON));
//assert(approx_equal(pers_fov(3, 2), -2.0f * far_p * near_p / (far_p - near_p), SMALL_EPSILON));
//assert(approx_equal(pers_fov(2, 3), -1.0f));

        // fromPerspective (frustum)
        float left = -1, right = 1, bottom = -1, top = 1;
        auto pers_frust = fmat4_t<V>::fromPerspective(left, right, bottom, top, near_p, far_p);
//assert(approx_equal(pers_frust(0, 0), 2.0f * near_p / (right - left), SMALL_EPSILON));
//assert(approx_equal(pers_frust(2, 3), -1.0f));

        // fromOrtho (frustum)
        auto ortho_frust = fmat4_t<V>::fromOrtho(left, right, bottom, top, near_p, far_p);
        assert(approx_equal(ortho_frust(0, 0), 2.0f / (right - left), SMALL_EPSILON));
// assert(approx_equal(ortho_frust(2, 2), -2.0f / (far_p - near_p), SMALL_EPSILON));
// assert(approx_equal(ortho_frust(3, 2), -(far_p + near_p) / (far_p - near_p), SMALL_EPSILON));

        // fromOrtho (width/height)
        float width = 2.0f, height = 2.0f;
        auto ortho_wh = fmat4_t<V>::fromOrtho(width, height, near_p, far_p);
        assert(approx_equal(ortho_wh(0, 0), 2.0f / width, SMALL_EPSILON));

        // fromLookAt
        fvec3 eye(0, 0, 5), target(0, 0, 0), up(0, 1, 0);
        auto look = fmat4_t<V>::fromLookAt(eye, target, up);
        assert(vec3_approx_equal(look.Forward(), fvec3(0, 0, -1).Normalize(), SMALL_EPSILON));
        assert(vec3_approx_equal(look.Up(), up.Normalize(), SMALL_EPSILON));
        assert(vec3_approx_equal(look.ExtractPosition(), eye));

        // fromBillboard
        fvec3 from(0, 0, 0), to(0, 0, 1);
        auto bill = fmat4_t<V>::fromBillboard(from, to, up);
        assert(vec3_approx_equal(bill.Forward(), (to - from).Normalize(), SMALL_EPSILON));
        assert(vec3_approx_equal(bill.ExtractPosition(), from));

        // fromRotationX
        auto rot_x = fmat4_t<V>::fromRotationX(xmath::pi_over2_v);
        assert(vec3_approx_equal(rot_x.RotateVector(fvec3(0, 1, 0)), fvec3(0, 0, 1), SMALL_EPSILON));

        // fromRotationY
        auto rot_y = fmat4_t<V>::fromRotationY(xmath::pi_over2_v);
        assert(vec3_approx_equal(rot_y.RotateVector(fvec3(0, 0, 1)), fvec3(1, 0, 0), SMALL_EPSILON));

        // fromRotationZ
        auto rot_z = fmat4_t<V>::fromRotationZ(xmath::pi_over2_v);
        assert(vec3_approx_equal(rot_z.RotateVector(fvec3(1, 0, 0)), fvec3(0, 1, 0), SMALL_EPSILON));
    }

    template <bool V>
    void test_setup_methods() {
        fmat4_t<V> m;

        // setupIdentity
        m.setupIdentity();
        assert(m.isIdentity());

        // setupZero
        m.setupZero();
        assert(mat4_approx_equal(m, fmat4_t<V>::fromZero()));

        // setupTranslation
        fvec3 t(1, 2, 3);
        m.setupTranslation(t);
        assert(vec3_approx_equal(m.ExtractPosition(), t));
        assert(vec3_approx_equal(m.ExtractScale(), fvec3(1.0f)));

        // setupRotation (quat)
        fquat q = fquat::fromAxisAngle(fvec3(1, 0, 0), xmath::pi_over4_v);
        m.setupRotation(q);
        assert(quat_approx_equal(m.ExtractRotation(), q, SMALL_EPSILON));

        // setupRotation (euler)
        radian3 euler(xmath::pi_over4_v, xmath::pi_over2_v, xmath::pi_v);
        m.setupRotation(euler);
        fmat4_t<V> expected_euler = fmat4_t<V>::fromRotationY(euler.m_Yaw) * fmat4_t<V>::fromRotationX(euler.m_Pitch) * fmat4_t<V>::fromRotationZ(euler.m_Roll);
        assert(mat4_approx_equal(m, expected_euler, SMALL_EPSILON));

        // setupScale (vec3)
        fvec3 s(2, 3, 4);
        m.setupScale(s);
        assert(vec3_approx_equal(m.ExtractScale(), s));

        // setupScale (float)
        m.setupScale(5.0f);
        assert(vec3_approx_equal(m.ExtractScale(), fvec3(5.0f)));

        // setup (TRS)
        fvec3 trans(1, 2, 3), scale(4, 5, 6);
        fquat rot = fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v);
        m.setup(trans, rot, scale);
        assert(vec3_approx_equal(m.ExtractPosition(), trans));
        assert(quat_approx_equal(m.ExtractRotation(), rot, SMALL_EPSILON));
        assert(vec3_approx_equal(m.ExtractScale(), scale));
    }

    template <bool V>
    void test_accessors()
    {
        auto m = fmat4_t<V>::fromIdentity();
        assert(vec4_approx_equal(m[0], fvec4(1, 0, 0, 0)));
        assert(vec4_approx_equal(m[1], fvec4(0, 1, 0, 0)));
        assert(vec4_approx_equal(m[2], fvec4(0, 0, 1, 0)));
        assert(vec4_approx_equal(m[3], fvec4(0, 0, 0, 1)));

        assert(approx_equal(m(0, 0), 1.0f));
        m(0, 0) = 2.0f;
        assert(approx_equal(m(0, 0), 2.0f));

        std::span<const float, 16> span = m;
        assert(approx_equal(span[0], 2.0f));

        fquat q_from_m = static_cast<fquat>(m);
        assert(quat_approx_equal(q_from_m, fquat::fromIdentity()));
    }

    template <bool V>
    void test_operations()
    {
        auto id = fmat4_t<V>::fromIdentity();
        auto trans = fmat4_t<V>::fromTranslation(fvec3(1, 1, 1));
        auto scale = fmat4_t<V>::fromScale(fvec3(2, 2, 2));

        // +
        auto sum = id + trans;
        assert(mat4_approx_equal(sum, trans + id));
        assert(vec3_approx_equal(sum.ExtractPosition(), fvec3(1, 1, 1)));

        // -
        auto diff = trans - id;
        assert(mat4_approx_equal(diff * id, diff));

        // *
        auto prod = id * trans;
        assert(mat4_approx_equal(prod, trans));

        auto prod_scale_trans = scale * trans;
        assert(vec3_approx_equal(prod_scale_trans.TransformPosition(fvec3(0)), fvec3(2, 2, 2)));

        // += 
        auto m1 = id;
        m1 += trans;
        assert(mat4_approx_equal(m1, id + trans));

        // -=
        m1 = trans;
        m1 -= id;
        assert(mat4_approx_equal(m1, trans - id));

        // *=
        m1 = scale;
        m1 *= trans;
        assert(mat4_approx_equal(m1, scale * trans));

        // * vec4
        fvec4 v(1, 2, 3, 1);
        auto transformed = trans * v;
        assert(vec4_approx_equal(transformed, fvec4(2, 3, 4, 1)));

        // * vec3 (implies w=1)
        fvec3 p(1, 2, 3);
        auto transformed_p = trans * p;
        assert(vec3_approx_equal(transformed_p, fvec3(2, 3, 4)));

        // Equals
        assert(id.Equals(id, 0.0f));
        assert(id.Equals(fmat4_t<V>(1.000001f), EPSILON));
        assert(!id.Equals(trans, 0.0f));
    }

    template <bool V>
    void test_math_functions() {
        auto id = fmat4_t<V>::fromIdentity();
        assert(mat4_approx_equal(id.Transpose(), id));
        assert(mat4_approx_equal(id.Inverse(), id));
        assert(approx_equal(id.Determinant(), 1.0f));

        auto scale_mat = fmat4_t<V>::fromScale(fvec3(2, 3, 4));
        assert(approx_equal(scale_mat.Determinant(), 2 * 3 * 4, SMALL_EPSILON));
        assert(mat4_approx_equal(scale_mat.Inverse(), fmat4_t<V>::fromScale(fvec3(1 / 2.f, 1 / 3.f, 1 / 4.f)), SMALL_EPSILON));

        // SRT inverse
        auto srt = fmat4_t<V>::fromScale(fvec3(2.0f)) * fmat4_t<V>::fromRotation(fvec3(0, 1, 0), xmath::pi_over2_v) * fmat4_t<V>::fromTranslation(fvec3(1, 0, 0));
        auto inv_srt = srt.InverseSRT();
        assert(mat4_approx_equal(inv_srt * srt, id, SMALL_EPSILON));

        // RT inverse
        auto rt = fmat4_t<V>::fromRotation(fvec3(0, 1, 0), xmath::pi_over2_v) * fmat4_t<V>::fromTranslation(fvec3(1, 0, 0));
        auto inv_rt = rt.InverseRT();
        assert(mat4_approx_equal(inv_rt * rt, id, SMALL_EPSILON));

        // Orthogonalize
        auto ortho = fmat4_t<V>::fromIdentity();
        ortho(0, 1) = 0.5f; // Skew
        ortho.Orthogonalize();
        assert(approx_equal(ortho.Right().Dot(ortho.Up()), 0.0f, SMALL_EPSILON));
        assert(ortho.Right().isNormalized(SMALL_EPSILON));

        // Transpose non-symmetric
        auto non_sym = fmat4_t<V>(std::array<float, 16>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16});
        auto transp = non_sym.Transpose();
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) assert(approx_equal(non_sym(i, j), transp(j, i)));
    }

    template <bool V>
    void test_geometry_helpers() {
        fvec3 pos(1, 2, 3), sc(4, 5, 6);
        fquat rot = fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v);
        auto m = fmat4_t<V>(pos, rot, sc);
        assert(vec3_approx_equal(m.ExtractPosition(), pos));
        assert(quat_approx_equal(m.ExtractRotation(), rot, SMALL_EPSILON));
        assert(vec3_approx_equal(m.ExtractScale(), sc));

        assert(vec3_approx_equal(m.Forward(), rot.RotateVector(fvec3(0, 0, 1)), SMALL_EPSILON));
        assert(vec3_approx_equal(m.Back(), rot.RotateVector(fvec3(0, 0, -1)), SMALL_EPSILON));
        assert(vec3_approx_equal(m.Up(), rot.RotateVector(fvec3(0, 1, 0)), SMALL_EPSILON));
        assert(vec3_approx_equal(m.Down(), rot.RotateVector(fvec3(0, -1, 0)), SMALL_EPSILON));
        assert(vec3_approx_equal(m.Left(), rot.RotateVector(fvec3(-1, 0, 0)), SMALL_EPSILON));
        assert(vec3_approx_equal(m.Right(), rot.RotateVector(fvec3(1, 0, 0)), SMALL_EPSILON));

        fvec3 v(1, 0, 0);
        assert(vec3_approx_equal(m.ExtractRotation().RotateVector(v), rot.RotateVector(v), SMALL_EPSILON));
        assert(vec3_approx_equal(m.InvRotateVector(m.RotateVector(v)), v, SMALL_EPSILON));

        assert(vec3_approx_equal(m.TransformPosition(fvec3(0)), pos));
        assert(vec3_approx_equal(m.TransformDirection(v), m.RotateVector(v), SMALL_EPSILON));
    }

    template <bool V>
    void test_chaining() {
        auto m = fmat4_t<V>::fromIdentity();
        m.Scale(fvec3(2.0f)).Rotate(fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v)).Translate(fvec3(1, 0, 0));
        assert(vec3_approx_equal(m.ExtractPosition(), fvec3(1, 0, 0)));
        assert(vec3_approx_equal(m.RotateVector(fvec3(1, 0, 0)), fvec3(0, 0, -2), SMALL_EPSILON));

        m = fmat4_t<V>::fromIdentity();
        m.RotateX(xmath::pi_over2_v);
        assert(vec3_approx_equal(m * fvec3(0, 1, 0), fvec3(0, 0, 1), SMALL_EPSILON));

        m = fmat4_t<V>::fromIdentity();
        m.RotateY(xmath::pi_over2_v);
        assert(vec3_approx_equal(m * fvec3(0, 0, 1), fvec3(1, 0, 0), SMALL_EPSILON));

        m = fmat4_t<V>::fromIdentity();
        m.RotateZ(xmath::pi_over2_v);
        assert(vec3_approx_equal(m * fvec3(1, 0, 0), fvec3(0, 1, 0), SMALL_EPSILON));

        // Combined
        m = fmat4_t<V>::fromIdentity();
        m.RotateZ(xmath::pi_over2_v).RotateY(xmath::pi_over2_v).RotateX(xmath::pi_over2_v);
        assert(vec3_approx_equal(m * fvec3(0, 1, 0), fvec3(0, -1, 0), SMALL_EPSILON));

        m.Rotate(fvec3(1, 0, 0), xmath::pi_over2_v);
        assert(vec3_approx_equal(m * fvec3(0, 1, 0), fvec3(1, 0, 0), SMALL_EPSILON));
    }

    template <bool V>
    void test_pre_chaining() {
        auto m = fmat4_t<V>::fromIdentity();
        m.Scale(fvec3(2.0f)).Rotate(fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v)).Translate(fvec3(1, 0, 0));
        assert(vec3_approx_equal(m.TransformPosition(fvec3(0)), fvec3(1, 0, 0)));
        assert(vec3_approx_equal(m * fvec3(1, 0, 0), fvec3(0, 0, -2) + fvec3(1, 0, 0), SMALL_EPSILON));

        // Similar for PreRotateX etc.
        m = fmat4_t<V>::fromIdentity();
        m.PreRotateX(xmath::pi_over2_v);
        assert(vec3_approx_equal(m * fvec3(0, 1, 0), fvec3(0, 0, 1), SMALL_EPSILON));
    }

    template <bool V>
    void test_clear_methods() {
        auto m = fmat4_t<V>(fvec3(1, 2, 3), fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v), fvec3(4, 5, 6));
        m.ClearTranslation();
        assert(vec3_approx_equal(m.ExtractPosition(), fvec3(0)));

        m.ClearRotation();
        assert(quat_approx_equal(m.ExtractRotation(), fquat::fromIdentity()));

        m.ClearScale();
        assert(vec3_approx_equal(m.ExtractScale(), fvec3(1.0f)));
    }

    template <bool V>
    void test_copy_methods() {
        auto m = fmat4_t<V>::fromIdentity();
        auto trans_copy = m.TranslateCopy(fvec3(1, 0, 0));
        assert(vec3_approx_equal(trans_copy.ExtractPosition(), fvec3(1, 0, 0)));
        assert(m.isIdentity());

        auto rotate_copy = m.RotateCopy(fquat::fromAxisAngle(fvec3(0, 1, 0), xmath::pi_over2_v));
        assert(vec3_approx_equal(rotate_copy.RotateVector(fvec3(1, 0, 0)), fvec3(0, 0, -1), SMALL_EPSILON));
        assert(m.isIdentity());

        auto rotate_aa_copy = m.RotateCopy(fvec3(0, 1, 0), xmath::pi_over2_v);
        assert(vec3_approx_equal(rotate_aa_copy.RotateVector(fvec3(1, 0, 0)), fvec3(0, 0, -1), SMALL_EPSILON));

        auto scale_copy = m.ScaleCopy(fvec3(2.0f));
        assert(vec3_approx_equal(scale_copy.ExtractScale(), fvec3(2.0f)));

        auto pre_trans_copy = m.PreTranslateCopy(fvec3(1, 0, 0));
        assert(vec3_approx_equal(pre_trans_copy.ExtractPosition(), fvec3(1, 0, 0)));

        // Similar for others
    }

    template <bool V>
    void test_safety_validation() {
        auto m = fmat4_t<V>::fromIdentity();
        assert(m.isFinite());
        assert(m.isIdentity());

        m(0, 0) = std::numeric_limits<float>::infinity();
        assert(!m.isFinite());

        auto non_id = fmat4_t<V>::fromTranslation(fvec3(1));
        assert(!non_id.isIdentity());

        // SanityCheck: assuming it throws or asserts on invalid, but since noexcept, perhaps no-op or check NaN/inf
        // Test if it doesn't crash
        m.SanityCheck(); // If it throws, test will fail anyway
    }

    template <bool V>
    void run_all_tests()
    {
        test_constructors<V>();
        test_static_constructors<V>();
        test_setup_methods<V>();
        test_accessors<V>();
        test_operations<V>();
        test_math_functions<V>();
        test_geometry_helpers<V>();
        test_chaining<V>();
        test_pre_chaining<V>();
        test_clear_methods<V>();
        test_copy_methods<V>();
        //test_safety_validation<V>();
    }

    int RunTests()
    {
        run_all_tests<true>();
        run_all_tests<false>();
        return 0;
    }

} // end of namespace
