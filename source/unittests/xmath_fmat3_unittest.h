#include <cassert>
#include <array>
#include <span>
#include <cmath>

namespace xmath
{
    namespace unit_test::_fmat3
    {
        template <bool V>
        void TestConstructorDiagonal() {
            fmat3_t<V> m(2.0f);
            assert(m(0, 0) == 2.0f);
            assert(m(0, 1) == 0.0f);
            assert(m(0, 2) == 0.0f);
            assert(m(1, 0) == 0.0f);
            assert(m(1, 1) == 2.0f);
            assert(m(1, 2) == 0.0f);
            assert(m(2, 0) == 0.0f);
            assert(m(2, 1) == 0.0f);
            assert(m(2, 2) == 2.0f);
        }

        template <bool V>
        void TestConstructorArray() {
            std::array<float, 9> arr = { 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f };
            fmat3_t<V> m(arr);
            assert(m(0, 0) == 1.f);
            assert(m(0, 1) == 4.f);
            assert(m(0, 2) == 7.f);
            assert(m(1, 0) == 2.f);
            assert(m(1, 1) == 5.f);
            assert(m(1, 2) == 8.f);
            assert(m(2, 0) == 3.f);
            assert(m(2, 1) == 6.f);
            assert(m(2, 2) == 9.f);
        }

        template <bool V>
        void TestConstructorSpan() {
            std::array<float, 9> arr = { 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f };
            auto xx = std::span< const float, 9 >(arr);
            fmat3_t<V> Mat(xx);

            assert(Mat(0, 0) == 1.f);
            assert(Mat(0, 1) == 4.f);
            assert(Mat(0, 2) == 7.f);
            assert(Mat(1, 0) == 2.f);
            assert(Mat(1, 1) == 5.f);
            assert(Mat(1, 2) == 8.f);
            assert(Mat(2, 0) == 3.f);
            assert(Mat(2, 1) == 6.f);
            assert(Mat(2, 2) == 9.f);
        }

        template <bool V>
        void TestConstructorQuaternion() {
            fquat q(0.f, 0.f, 0.f, 1.f);  // identity
            fmat3_t<V> m(q);
            assert(m.isIdentity());

            fquat q2(1.f, 0.f, 0.f, 0.f);  // 180 deg around x
            fmat3_t<V> m2(q2);
            assert(m2(0, 0) == 1.f);
            assert(m2(0, 1) == 0.f);
            assert(m2(0, 2) == 0.f);
            assert(m2(1, 0) == 0.f);
            assert(m2(1, 1) == -1.f);
            assert(m2(1, 2) == 0.f);
            assert(m2(2, 0) == 0.f);
            assert(m2(2, 1) == 0.f);
            assert(m2(2, 2) == -1.f);
        }

        template <bool V>
        void TestConstructorEuler()
        {
            radian3 euler(radian3::fromZero());
            fmat3_t<V> m(euler);
            assert(m.isIdentity());

            // Test a known Euler rotation, e.g., 90 deg around x
            radian3 euler_x(radian(1.5707963f), radian(0.f), radian(0.f));
            fmat3_t<V> mx(euler_x);
            assert(std::abs(mx(0, 0) - 1.f) < 1e-6f);
            assert(std::abs(mx(0, 1) - 0.f) < 1e-6f);
            assert(std::abs(mx(0, 2) - 0.f) < 1e-6f);
            assert(std::abs(mx(1, 0) - 0.f) < 1e-6f);
            assert(std::abs(mx(1, 1) - 0.f) < 1e-6f);
            assert(std::abs(mx(1, 2) + 1.f) < 1e-6f);
            assert(std::abs(mx(2, 0) - 0.f) < 1e-6f);
            assert(std::abs(mx(2, 1) - 1.f) < 1e-6f);
            assert(std::abs(mx(2, 2) - 0.f) < 1e-6f);
        }

        template <bool V>
        void TestConstructorRotationScale() {
            fquat q(0.f, 0.f, 0.f, 1.f);
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m(q, s);
            assert(m(0, 0) == 2.f);
            assert(m(1, 1) == 3.f);
            assert(m(2, 2) == 4.f);
            assert(m.ExtractScale() == s);
            assert(m.ExtractRotation() == q);
        }

        template <bool V>
        void TestConstructorCopyOtherSIMD() {
            fmat3_t<!V> other = fmat3_t<!V>::fromIdentity();
            fmat3_t<V> m(other);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestFromIdentity() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            assert(m.isIdentity());
        }

        template <bool V>
        void TestFromZero() {
            fmat3_t<V> m = fmat3_t<V>::fromZero();
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    assert(m(i, j) == 0.0f);
                }
            }
        }

        template <bool V>
        void TestFromRotationQuaternion() {
            fquat q(0.f, 0.f, 0.f, 1.f);
            fmat3_t<V> m = fmat3_t<V>::fromRotation(q);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestFromRotationAxisAngle() {
            fvec3 axis(1.f, 0.f, 0.f);
            radian angle(1.5707963f);
            fmat3_t<V> m = fmat3_t<V>::fromRotation(axis, angle);
            assert(std::abs(m(0, 0) - 1.f) < 1e-6f);
            assert(std::abs(m(0, 1) - 0.f) < 1e-6f);
            assert(std::abs(m(0, 2) - 0.f) < 1e-6f);
            assert(std::abs(m(1, 0) - 0.f) < 1e-6f);
            assert(std::abs(m(1, 1) - 0.f) < 1e-6f);
            assert(std::abs(m(1, 2) + 1.f) < 1e-6f);
            assert(std::abs(m(2, 0) - 0.f) < 1e-6f);
            assert(std::abs(m(2, 1) - 1.f) < 1e-6f);
            assert(std::abs(m(2, 2) - 0.f) < 1e-6f);
        }

        template <bool V>
        void TestFromScale() {
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m = fmat3_t<V>::fromScale(s);
            assert(m(0, 0) == 2.f);
            assert(m(1, 1) == 3.f);
            assert(m(2, 2) == 4.f);
            assert(m(0, 1) == 0.f);
        }

        template <bool V>
        void TestFromRotationX() {
            radian angle(1.5707963f); // pi/2
            fmat3_t<V> m = fmat3_t<V>::fromRotationX(angle);
            assert(std::abs(m(0, 0) - 1.f) < 1e-6f);
            assert(std::abs(m(1, 1) - 0.f) < 1e-6f);
            assert(std::abs(m(1, 2) + 1.f) < 1e-6f);
            assert(std::abs(m(2, 1) - 1.f) < 1e-6f);
            assert(std::abs(m(2, 2) - 0.f) < 1e-6f);
        }

        template <bool V>
        void TestFromRotationY() {
            radian angle(1.5707963f);
            fmat3_t<V> m = fmat3_t<V>::fromRotationY(angle);
            assert(std::abs(m(0, 0) - 0.f) < 1e-6f);
            assert(std::abs(m(0, 2) + 1.f) < 1e-6f);
            assert(std::abs(m(1, 1) - 1.f) < 1e-6f);
            assert(std::abs(m(2, 0) - 1.f) < 1e-6f);
            assert(std::abs(m(2, 2) - 0.f) < 1e-6f);
        }

        template <bool V>
        void TestFromRotationZ() {
            radian angle(1.5707963f);
            fmat3_t<V> m = fmat3_t<V>::fromRotationZ(angle);
            assert(std::abs(m(0, 0) - 0.f) < 1e-6f);
            assert(std::abs(m(0, 1) + 1.f) < 1e-6f);
            assert(std::abs(m(1, 0) - 1.f) < 1e-6f);
            assert(std::abs(m(1, 1) - 0.f) < 1e-6f);
            assert(std::abs(m(2, 2) - 1.f) < 1e-6f);
        }

        template <bool V>
        void TestSetupIdentity() {
            fmat3_t<V> m;
            m.setupIdentity();
            assert(m.isIdentity());
        }

        template <bool V>
        void TestSetupZero() {
            fmat3_t<V> m;
            m.setupZero();
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    assert(m(i, j) == 0.0f);
                }
            }
        }

        template <bool V>
        void TestSetupRotationQuaternion() {
            fquat q(0.f, 0.f, 0.f, 1.f);
            fmat3_t<V> m;
            m.setupRotation(q);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestSetupRotationEuler() {
            radian3 euler(radian3::fromZero());
            fmat3_t<V> m;
            m.setupRotation(euler);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestSetupScaleVector() {
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m;
            m.setupScale(s);
            assert(m(0, 0) == 2.f);
            assert(m(1, 1) == 3.f);
            assert(m(2, 2) == 4.f);
        }

        template <bool V>
        void TestSetupScaleFloat() {
            fmat3_t<V> m;
            m.setupScale(5.f);
            assert(m(0, 0) == 5.f);
            assert(m(1, 1) == 5.f);
            assert(m(2, 2) == 5.f);
        }

        template <bool V>
        void TestSetupRotationScale() {
            fquat q(0.f, 0.f, 0.f, 1.f);
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m;
            m.setup(q, s);
            assert(m(0, 0) == 2.f);
            assert(m(1, 1) == 3.f);
            assert(m(2, 2) == 4.f);
        }

        template <bool V>
        void TestOperatorBracket() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 row0 = m[0];
            assert(row0.m_X == 1.f);
            assert(row0.m_Y == 0.f);
            assert(row0.m_Z == 0.f);
        }

        template <bool V>
        void TestOperatorParen() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            assert(m(0, 0) == 1.f);
            m(0, 0) = 5.f;
            assert(m(0, 0) == 5.f);
        }

        template <bool V>
        void TestOperatorSpan() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            std::span<const float, 9> sp = m;
            assert(std::abs(sp[0] - 1.f) < 1e-6f); // m_00
            assert(std::abs(sp[1] - 0.f) < 1e-6f); // m_10
            assert(std::abs(sp[2] - 0.f) < 1e-6f); // m_20
            if constexpr (V)
            {
                assert(std::abs(sp[3] - 0.f) < 1e-6f); // padding
                assert(std::abs(sp[4] - 0.f) < 1e-6f); // m_01
                assert(std::abs(sp[5] - 1.f) < 1e-6f); // m_11
                assert(std::abs(sp[6] - 0.f) < 1e-6f); // m_21
                assert(std::abs(sp[7] - 0.f) < 1e-6f); // padding
                assert(std::abs(sp[8] - 0.f) < 1e-6f); // m_02
            }
            else
            {
                assert(std::abs(sp[3] - 0.f) < 1e-6f); // m_01
                assert(std::abs(sp[4] - 1.f) < 1e-6f); // m_11
                assert(std::abs(sp[5] - 0.f) < 1e-6f); // m_21
                assert(std::abs(sp[6] - 0.f) < 1e-6f); // m_02
                assert(std::abs(sp[7] - 0.f) < 1e-6f); // m_12
                assert(std::abs(sp[8] - 1.f) < 1e-6f); // m_22
            }
        }

        template <bool V>
        void TestOperatorFQuat() {
            fquat q(0.f, 0.f, 0.f, 1.f);
            fmat3_t<V> m(q);
            fquat q2 = m;
            assert(q2 == q);
        }

        template <bool V>
        void TestOperatorPlus() {
            fmat3_t<V> m1 = fmat3_t<V>::fromIdentity();
            fmat3_t<V> m2 = fmat3_t<V>::fromIdentity();
            fmat3_t<V> sum = m1 + m2;
            assert(sum(0, 0) == 2.f);
            assert(sum(1, 1) == 2.f);
            assert(sum(2, 2) == 2.f);
        }

        template <bool V>
        void TestOperatorMinus() {
            fmat3_t<V> m1 = fmat3_t<V>::fromIdentity();
            fmat3_t<V> m2 = fmat3_t<V>::fromIdentity();
            fmat3_t<V> diff = m1 - m2;
            assert(diff(0, 0) == 0.f);
            assert(diff(1, 1) == 0.f);
            assert(diff(2, 2) == 0.f);
        }

        template <bool V>
        void TestOperatorMultiplyMatrix() {
            fmat3_t<V> m1 = fmat3_t<V>::fromIdentity();
            fmat3_t<V> m2 = fmat3_t<V>::fromScale(fvec3(2.f));
            fmat3_t<V> prod = m1 * m2;
            assert(prod(0, 0) == 2.f);
            assert(prod(1, 1) == 2.f);
            assert(prod(2, 2) == 2.f);
        }

        template <bool V>
        void TestOperatorMultiplyVector() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 v(1.f, 2.f, 3.f);
            fvec3 res = m * v;
            assert(res == v);
        }

        template <bool V>
        void TestOperatorPlusEqual() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m += fmat3_t<V>::fromIdentity();
            assert(m(0, 0) == 2.f);
        }

        template <bool V>
        void TestOperatorMinusEqual() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m -= fmat3_t<V>::fromIdentity();
            assert(m(0, 0) == 0.f);
        }

        template <bool V>
        void TestOperatorMultiplyEqual() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m *= fmat3_t<V>::fromScale(fvec3(2.f));
            assert(m(0, 0) == 2.f);
        }

        template <bool V>
        void TestEquals() {
            fmat3_t<V> m1 = fmat3_t<V>::fromIdentity();
            fmat3_t<V> m2 = fmat3_t<V>::fromIdentity();
            assert(m1.Equals(m2, 1e-6f));
            m2(0, 0) = 1.00001f;
            assert(!m1.Equals(m2, 1e-6f));
            assert(m1.Equals(m2, 1e-4f));
        }

        template <bool V>
        void TestTranspose() {
            std::array<float, 9> arr = { 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f };
            fmat3_t<V> m(arr);
            fmat3_t<V> mt = m.Transpose();
            assert(mt(0, 0) == 1.f);
            assert(mt(0, 1) == 2.f);
            assert(mt(0, 2) == 3.f);
            assert(mt(1, 0) == 4.f);
            assert(mt(1, 1) == 5.f);
            assert(mt(1, 2) == 6.f);
            assert(mt(2, 0) == 7.f);
            assert(mt(2, 1) == 8.f);
            assert(mt(2, 2) == 9.f);
        }

        template <bool V>
        void TestDeterminant() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            assert(m.Determinant() == 1.f);

            fmat3_t<V> m2 = fmat3_t<V>::fromScale(fvec3(2.f, 3.f, 4.f));
            assert(m2.Determinant() == 24.f);
        }

        template <bool V>
        void TestInverse() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fmat3_t<V> inv = m.Inverse();
            assert(inv.isIdentity());

            std::array<float, 9> arr = { 1.f, 0.f, 5.f, 2.f, 1.f, 6.f, 3.f, 4.f, 0.f };
            fmat3_t<V> m3(arr);
            fmat3_t<V> inv3 = m3.Inverse();
            fmat3_t<V> prod = m3 * inv3;
            assert(prod.Equals(fmat3_t<V>::fromIdentity(), 1e-5f));
        }

        template <bool V>
        void TestOrthogonalize() {
            fmat3_t<V> m;
            m(0, 0) = 1.f; m(0, 1) = 0.1f; m(0, 2) = 0.f;
            m(1, 0) = 0.f; m(1, 1) = 1.f; m(1, 2) = 0.f;
            m(2, 0) = 0.f; m(2, 1) = 0.f; m(2, 2) = 1.f;
            m.Orthogonalize();
            assert(std::abs(m.Right().Dot(m.Up())) < 1e-6f);
            assert(std::abs(m.Right().Length() - 1.f) < 1e-6f);
            assert(std::abs(m.Up().Length() - 1.f) < 1e-6f);
            assert(std::abs(m.Forward().Length() - 1.f) < 1e-6f);
        }

        template <bool V>
        void TestExtractRotation() {
            fquat q(0.f, 0.f, 0.f, 1.f);
            fmat3_t<V> m(q);
            fquat extracted = m.ExtractRotation();
            assert(extracted == q);
        }

        template <bool V>
        void TestExtractScale() {
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m = fmat3_t<V>::fromScale(s);
            fvec3 extracted = m.ExtractScale();
            assert(extracted == s);
        }

        template <bool V>
        void TestDirections() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            assert(m.Forward() == fvec3(0.f, 0.f, 1.f));
            assert(m.Back() == fvec3(0.f, 0.f, -1.f));
            assert(m.Up() == fvec3(0.f, 1.f, 0.f));
            assert(m.Down() == fvec3(0.f, -1.f, 0.f));
            assert(m.Left() == fvec3(-1.f, 0.f, 0.f));
            assert(m.Right() == fvec3(1.f, 0.f, 0.f));
        }

        template <bool V>
        void TestRotateVector() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 v(1.f, 2.f, 3.f);
            assert(m.RotateVector(v) == v);
        }

        template <bool V>
        void TestInvRotateVector() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 v(1.f, 2.f, 3.f);
            assert(m.InvRotateVector(v) == v);
        }

        template <bool V>
        void TestTransformDirection() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 d(1.f, 0.f, 0.f);
            assert(m.TransformDirection(d) == d);
        }

        template <bool V>
        void TestRotateChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fquat q(0.f, 0.f, 0.f, 1.f);
            m.Rotate(q);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestRotateAxisAngleChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 axis(1.f, 0.f, 0.f);
            radian angle(0.f);
            m.Rotate(axis, angle);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestRotateXChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.RotateX(radian(0.f));
            assert(m.isIdentity());
        }

        template <bool V>
        void TestRotateYChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.RotateY(radian(0.f));
            assert(m.isIdentity());
        }

        template <bool V>
        void TestRotateZChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.RotateZ(radian(0.f));
            assert(m.isIdentity());
        }

        template <bool V>
        void TestScaleVectorChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 s(2.f, 2.f, 2.f);
            m.Scale(s);
            assert(m(0, 0) == 2.f);
        }

        template <bool V>
        void TestScaleFloatChaining() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.Scale(3.f);
            assert(m(0, 0) == 3.f);
        }

        template <bool V>
        void TestPreRotate() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fquat q(0.f, 0.f, 0.f, 1.f);
            m.PreRotate(q);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestPreRotateAxisAngle() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 axis(1.f, 0.f, 0.f);
            radian angle(0.f);
            m.PreRotate(axis, angle);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestPreRotateX() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.PreRotateX(radian(0.f));
            assert(m.isIdentity());
        }

        template <bool V>
        void TestPreRotateY() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.PreRotateY(radian(0.f));
            assert(m.isIdentity());
        }

        template <bool V>
        void TestPreRotateZ() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.PreRotateZ(radian(0.f));
            assert(m.isIdentity());
        }

        template <bool V>
        void TestPreScaleVector() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 s(2.f, 2.f, 2.f);
            m.PreScale(s);
            assert(m(0, 0) == 2.f);
        }

        template <bool V>
        void TestPreScaleFloat() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.PreScale(3.f);
            assert(m(0, 0) == 3.f);
        }

        template <bool V>
        void TestClearRotation() {
            fquat q(1.f, 0.f, 0.f, 0.f);
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m(q, s);
            m.ClearRotation();
            assert(m.Equals(fmat3_t<V>::fromScale(s), 1e-6f));
        }

        template <bool V>
        void TestClearScale() {
            fquat q(1.f, 0.f, 0.f, 0.f);
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m(q, s);
            m.ClearScale();
            assert(m.Equals(fmat3_t<V>::fromRotation(q), 1e-6f));
        }

        template <bool V>
        void TestRotateCopy() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fquat q(0.f, 0.f, 0.f, 1.f);
            fmat3_t<V> copy = m.RotateCopy(q);
            assert(copy.isIdentity());
            assert(m.isIdentity());  // original unchanged
        }

        template <bool V>
        void TestRotateCopyAxisAngle() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 axis(1.f, 0.f, 0.f);
            radian angle(0.f);
            fmat3_t<V> copy = m.RotateCopy(axis, angle);
            assert(copy.isIdentity());
        }

        template <bool V>
        void TestScaleCopy() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 s(2.f, 2.f, 2.f);
            fmat3_t<V> copy = m.ScaleCopy(s);
            assert(copy(0, 0) == 2.f);
            assert(m.isIdentity());
        }

        template <bool V>
        void TestPreRotateCopy() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fquat q(0.f, 0.f, 0.f, 1.f);
            fmat3_t<V> copy = m.PreRotateCopy(q);
            assert(copy.isIdentity());
        }

        template <bool V>
        void TestPreRotateCopyAxisAngle() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 axis(1.f, 0.f, 0.f);
            radian angle(0.f);
            fmat3_t<V> copy = m.PreRotateCopy(axis, angle);
            assert(copy.isIdentity());
        }

        template <bool V>
        void TestPreScaleCopy() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fvec3 s(2.f, 2.f, 2.f);
            fmat3_t<V> copy = m.PreScaleCopy(s);
            assert(copy(0, 0) == 2.f);
        }

        template <bool V>
        void TestIsFinite() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            assert(m.isFinite());
            m(0, 0) = std::numeric_limits<float>::infinity();
            assert(!m.isFinite());
        }

        template <bool V>
        void TestIsIdentity() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            assert(m.isIdentity());
            m(0, 0) = 1.00001f;
            assert(!m.isIdentity());
        }

        template <bool V>
        void TestSanityCheck() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            m.SanityCheck();  // Should not assert
        }

        template <bool V>
        void TestInverseIdentity() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fmat3_t<V> inv = m.Inverse();
            assert(inv.isIdentity());
        }

        template <bool V>
        void TestTransposeIdentity() {
            fmat3_t<V> m = fmat3_t<V>::fromIdentity();
            fmat3_t<V> trans = m.Transpose();
            assert(trans.isIdentity());
        }

        template <bool V>
        void TestDeterminantScale() {
            fvec3 s(2.f, 3.f, 4.f);
            fmat3_t<V> m = fmat3_t<V>::fromScale(s);
            assert(std::abs(m.Determinant() - 24.f) < 1e-6f);
        }

        int RunTests()
        {
            // Run tests for V=false
            TestConstructorDiagonal<false>();
            TestConstructorArray<false>();
            TestConstructorSpan<false>();
            TestConstructorQuaternion<false>();
            TestConstructorEuler<false>();
            TestConstructorRotationScale<false>();
            TestConstructorCopyOtherSIMD<false>();
            TestFromIdentity<false>();
            TestFromZero<false>();
            TestFromRotationQuaternion<false>();
            TestFromRotationAxisAngle<false>();
            TestFromScale<false>();
            TestFromRotationX<false>();
            TestFromRotationY<false>();
            TestFromRotationZ<false>();
            TestSetupIdentity<false>();
            TestSetupZero<false>();
            TestSetupRotationQuaternion<false>();
            TestSetupRotationEuler<false>();
            TestSetupScaleVector<false>();
            TestSetupScaleFloat<false>();
            TestSetupRotationScale<false>();
            TestOperatorBracket<false>();
            TestOperatorParen<false>();
            TestOperatorSpan<false>();
            TestOperatorFQuat<false>();
            TestOperatorPlus<false>();
            TestOperatorMinus<false>();
            TestOperatorMultiplyMatrix<false>();
            TestOperatorMultiplyVector<false>();
            TestOperatorPlusEqual<false>();
            TestOperatorMinusEqual<false>();
            TestOperatorMultiplyEqual<false>();
            TestEquals<false>();
            TestTranspose<false>();
            TestDeterminant<false>();
            TestInverse<false>();
            TestOrthogonalize<false>();
            TestExtractRotation<false>();
            TestExtractScale<false>();
            TestDirections<false>();
            TestRotateVector<false>();
            TestInvRotateVector<false>();
            TestTransformDirection<false>();
            TestRotateChaining<false>();
            TestRotateAxisAngleChaining<false>();
            TestRotateXChaining<false>();
            TestRotateYChaining<false>();
            TestRotateZChaining<false>();
            TestScaleVectorChaining<false>();
            TestScaleFloatChaining<false>();
            TestPreRotate<false>();
            TestPreRotateAxisAngle<false>();
            TestPreRotateX<false>();
            TestPreRotateY<false>();
            TestPreRotateZ<false>();
            TestPreScaleVector<false>();
            TestPreScaleFloat<false>();
            TestClearRotation<false>();
            TestClearScale<false>();
            TestRotateCopy<false>();
            TestRotateCopyAxisAngle<false>();
            TestScaleCopy<false>();
            TestPreRotateCopy<false>();
            TestPreRotateCopyAxisAngle<false>();
            TestPreScaleCopy<false>();
            TestIsFinite<false>();
            TestIsIdentity<false>();
            TestSanityCheck<false>();
            TestInverseIdentity<false>();
            TestTransposeIdentity<false>();
            TestDeterminantScale<false>();

            // Run tests for V=true
            TestConstructorDiagonal<true>();
            TestConstructorArray<true>();
            TestConstructorSpan<true>();
            TestConstructorQuaternion<true>();
            TestConstructorEuler<true>();
            TestConstructorRotationScale<true>();
            TestConstructorCopyOtherSIMD<true>();
            TestFromIdentity<true>();
            TestFromZero<true>();
            TestFromRotationQuaternion<true>();
            TestFromRotationAxisAngle<true>();
            TestFromScale<true>();
            TestFromRotationX<true>();
            TestFromRotationY<true>();
            TestFromRotationZ<true>();
            TestSetupIdentity<true>();
            TestSetupZero<true>();
            TestSetupRotationQuaternion<true>();
            TestSetupRotationEuler<true>();
            TestSetupScaleVector<true>();
            TestSetupScaleFloat<true>();
            TestSetupRotationScale<true>();
            TestOperatorBracket<true>();
            TestOperatorParen<true>();
            TestOperatorSpan<true>();
            TestOperatorFQuat<true>();
            TestOperatorPlus<true>();
            TestOperatorMinus<true>();
            TestOperatorMultiplyMatrix<true>();
            TestOperatorMultiplyVector<true>();
            TestOperatorPlusEqual<true>();
            TestOperatorMinusEqual<true>();
            TestOperatorMultiplyEqual<true>();
            TestEquals<true>();
            TestTranspose<true>();
            TestDeterminant<true>();
            TestInverse<true>();
            TestOrthogonalize<true>();
            TestExtractRotation<true>();
            TestExtractScale<true>();
            TestDirections<true>();
            TestRotateVector<true>();
            TestInvRotateVector<true>();
            TestTransformDirection<true>();
            TestRotateChaining<true>();
            TestRotateAxisAngleChaining<true>();
            TestRotateXChaining<true>();
            TestRotateYChaining<true>();
            TestRotateZChaining<true>();
            TestScaleVectorChaining<true>();
            TestScaleFloatChaining<true>();
            TestPreRotate<true>();
            TestPreRotateAxisAngle<true>();
            TestPreRotateX<true>();
            TestPreRotateY<true>();
            TestPreRotateZ<true>();
            TestPreScaleVector<true>();
            TestPreScaleFloat<true>();
            TestClearRotation<true>();
            TestClearScale<true>();
            TestRotateCopy<true>();
            TestRotateCopyAxisAngle<true>();
            TestScaleCopy<true>();
            TestPreRotateCopy<true>();
            TestPreRotateCopyAxisAngle<true>();
            TestPreScaleCopy<true>();
            TestIsFinite<true>();
            TestIsIdentity<true>();
            TestSanityCheck<true>();
            TestInverseIdentity<true>();
            TestTransposeIdentity<true>();
            TestDeterminantScale<true>();

            return 0;
        }

    }  // namespace unit_test::_fmat3
}  // namespace xmath