#pragma once

namespace xmath::unit_test::_ivec2
{
    void RunTests()
    {
        //------------------------------------------------------------------------------
        // Constructors
        //------------------------------------------------------------------------------
        {
            ivec2 v(1, 2);
            assert(v.m_X == 1 && v.m_Y == 2);
        }
        {
            ivec2 v(3);
            assert(v.m_X == 3 && v.m_Y == 3);
        }
        {
            std::array<int, 2> arr = { 4, 5 };
            std::span<int> sp(arr);
            ivec2 v(sp);
            assert(v.m_X == 4 && v.m_Y == 5);
        }
        {
            std::array<double, 2> arr = { 6.1, 7.9 };
            ivec2 v(arr);
            assert(v.m_X == 6 && v.m_Y == 7);
        }
        {
            std::array<float, 2> arr = { 8.2f, 9.8f };
            ivec2 v(arr);
            assert(v.m_X == 8 && v.m_Y == 9);
        }

        //------------------------------------------------------------------------------
        // Conversions
        //------------------------------------------------------------------------------
        {
            ivec2 v(10, 11);
            std::array<int, 2> arr = static_cast<std::array<int, 2>>(v);
            assert(arr[0] == 10 && arr[1] == 11);
        }
        {
            ivec2 v(12, 13);
            std::array<float, 2> arr = static_cast<std::array<float, 2>>(v);
            assert(arr[0] == 12.f && arr[1] == 13.f);
        }
        {
            ivec2 v(14, 15);
            std::string s = static_cast<std::string>(v);
            assert(s == "(14, 15)");
        }
        {
            ivec2 v(16, 17);
            assert(v.ToString() == "(16, 17)");
        }
        {
            ivec2 v(18, 19);
            std::ostringstream oss;
            oss << v;
            assert(oss.str() == "(18, 19)");
        }

        //------------------------------------------------------------------------------
        // Static properties
        //------------------------------------------------------------------------------
        {
            assert(ivec2::fromZero() == ivec2(0, 0));
            assert(ivec2::fromOne() == ivec2(1, 1));
            assert(ivec2::fromUnit() == ivec2(1, 1));
            assert(ivec2::fromUp() == ivec2(0, 1));
            assert(ivec2::fromDown() == ivec2(0, -1));
            assert(ivec2::fromLeft() == ivec2(-1, 0));
            assert(ivec2::fromRight() == ivec2(1, 0));
        }

        //------------------------------------------------------------------------------
        // Setup functions
        //------------------------------------------------------------------------------
        {
            ivec2 v;
            v.setupMax();
            assert(v.m_X == std::numeric_limits<int>::max() && v.m_Y == std::numeric_limits<int>::max());
        }
        {
            ivec2 v(20, 21);
            v.setupZero();
            assert(v == ivec2(0, 0));
        }
        {
            ivec2 v;
            v.setup(22, 23);
            assert(v == ivec2(22, 23));
        }

        //------------------------------------------------------------------------------
        // Static methods
        //------------------------------------------------------------------------------
        {
            ivec2 a(1, 2), b(3, 4);
            assert(ivec2::Dot(a, b) == 11LL);
        }
        {
            ivec2 a(1, 4), b(3, 2);
            assert(ivec2::Min(a, b) == ivec2(1, 2));
        }
        {
            ivec2 a(1, 4), b(3, 2);
            assert(ivec2::Max(a, b) == ivec2(3, 4));
        }
        {
            ivec2 a(0, 3), b(4, 0);
            assert(std::abs(ivec2::Distance(a, b) - 5.f) < 1e-5f);
        }
        {
            ivec2 a(1, 2), b(3, 4);
            assert(ivec2::Cross(a, b) == -2LL);
        }

        //------------------------------------------------------------------------------
        // Static methods as members
        //------------------------------------------------------------------------------
        {
            ivec2 v(1, 2), a(3, 4);
            assert(v.Dot(a) == 11LL);
        }
        {
            ivec2 v(1, 4), a(3, 2);
            assert(v.Min(a) == ivec2(1, 2));
        }
        {
            ivec2 v(1, 4), a(3, 2);
            assert(v.Max(a) == ivec2(3, 4));
        }
        {
            ivec2 v(0, 3), a(4, 0);
            assert(std::abs(v.Distance(a) - 5.f) < 1e-5f);
        }

        //------------------------------------------------------------------------------
        // Basic operations
        //------------------------------------------------------------------------------
        {
            ivec2 v(3, 4);
            assert(std::abs(v.Length() - 5.f) < 1e-5f);
        }
        {
            ivec2 v(3, 4);
            assert(v.LengthSquared() == 25LL);
        }
        {
            ivec2 v(3, 4);
            assert(v.LengthManhattan() == 7);
        }
        {
            ivec2 v(5, 10);
            assert(v.isInRange(0, 10));
            assert(!v.isInRange(6, 9));
        }
        {
            ivec2 v(1, 2), a(1, 2);
            assert(v.Equals(a));
            assert(!v.Equals(ivec2(1, 3)));
        }

        //------------------------------------------------------------------------------
        // Component-wise math
        //------------------------------------------------------------------------------
        {
            ivec2 v(-1, 2);
            assert(v.AbsCopy() == ivec2(1, 2));
        }
        {
            ivec2 v(-1, 2);
            v.Abs();
            assert(v == ivec2(1, 2));
        }
        {
            ivec2 v(-3, 0);
            assert(v.SignCopy() == ivec2(-1, 0));
        }
        {
            ivec2 v(4, -5);
            v.Sign();
            assert(v == ivec2(1, -1));
        }
        {
            ivec2 v(10, 7);
            assert(v.ModCopy(3) == ivec2(1, 1));
        }
        {
            ivec2 v(10, 7);
            v.Mod(3);
            assert(v == ivec2(1, 1));
        }
        {
            ivec2 v(5, 15);
            assert(v.ClampCopy(0, 10) == ivec2(5, 10));
        }
        {
            ivec2 v(5, 15);
            v.Clamp(0, 10);
            assert(v == ivec2(5, 10));
        }
        {
            ivec2 v(5, 15), min(0, 5), max(10, 20);
            assert(v.ClampCopy(min, max) == ivec2(5, 15));
        }
        {
            ivec2 v(5, 15), min(0, 5), max(10, 20);
            v.Clamp(min, max);
            assert(v == ivec2(5, 15));
        }

        //------------------------------------------------------------------------------
        // Geometry
        //------------------------------------------------------------------------------
        {
            ivec2 v(1, 2), a(3, 4);
            assert(v.DistanceSquare(a) == 8LL);
        }
        {
            ivec2 v(5, 9);
            v.GridSnap(2, 3);
            assert(v == ivec2(4, 9));
        }
        {
            ivec2 v(1, 2);
            assert(v.Perp() == ivec2(-2, 1));
        }
        {
            ivec2 p(3, 3), v0(0, 0), v1(2, 0);
            assert(p.WhichSideOfLine(v0, v1) == 6LL);
        }

        //------------------------------------------------------------------------------
        // Swizzle methods
        //------------------------------------------------------------------------------
        {
            ivec2 v(1, 2);
            assert(v.x() == 1);
            assert(v.y() == 2);
        }
        {
            ivec2 v(1, 2);
            assert(v.xx() == ivec2(1, 1));
            assert(v.xy() == ivec2(1, 2));
            assert(v.yx() == ivec2(2, 1));
            assert(v.yy() == ivec2(2, 2));
        }

        //------------------------------------------------------------------------------
        // Operator overloads
        //------------------------------------------------------------------------------
        {
            ivec2 a(1, 2), b(3, 4);
            assert(a + b == ivec2(4, 6));
        }
        {
            ivec2 a(1, 2), b(3, 4);
            assert(a - b == ivec2(-2, -2));
        }
        {
            ivec2 a(1, 2);
            assert(a * 3 == ivec2(3, 6));
        }
        {
            ivec2 a(6, 8);
            assert(a / 2 == ivec2(3, 4));
        }
        {
            ivec2 v(1, 2);
            v += ivec2(3, 4);
            assert(v == ivec2(4, 6));
        }
        {
            ivec2 v(1, 2);
            v -= ivec2(3, 4);
            assert(v == ivec2(-2, -2));
        }
        {
            ivec2 v(1, 2);
            v *= 3;
            assert(v == ivec2(3, 6));
        }
        {
            ivec2 v(6, 8);
            v /= 2;
            assert(v == ivec2(3, 4));
        }
        {
            ivec2 a(1, 2), b(1, 2);
            assert(a == b);
            assert(!(a != b));
        }
        {
            ivec2 v(1, 2);
            assert(v[0] == 1);
            assert(v[1] == 2);
            v[0] = 3;
            assert(v[0] == 3);
        }
        {
            ivec2 v(1, 2);
            assert(3 * v == ivec2(3, 6));
        }
        {
            ivec2 v(1, 2);
            assert(-v == ivec2(-1, -2));
        }
    }
}
