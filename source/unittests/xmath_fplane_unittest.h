#include <cassert>
#include <cmath>
#include <array>
#include <string>

namespace xmath::unit_test::_fplane
{
    // Helper functions (reused from fbbox unit tests)
    constexpr float EPSILON = 1e-5f;
    constexpr float SMALL_EPSILON = 1e-6f;
    bool approx_equal(float a, float b, float eps = EPSILON)
    {
        return std::fabs(a - b) < eps;
    }
    bool vec3_approx_equal(const fvec3& a, const fvec3& b, float eps = EPSILON)
    {
        return approx_equal(a.m_X, b.m_X, eps) && approx_equal(a.m_Y, b.m_Y, eps) && approx_equal(a.m_Z, b.m_Z, eps);
    }
    bool plane_approx_equal(const fplane_t<false>& a, const fplane_t<false>& b, float eps = EPSILON)
    {
        return vec3_approx_equal(a.Normal(), b.Normal(), eps) && approx_equal(a.Offset(), b.Offset(), eps);
    }

    // Tests for constructors
    template <bool V>
    void test_constructors()
    {
        // Default constructor (uninitialized, no checks as memory is undefined)
        fplane_t<V> default_plane{};

        // Components constructor
        fplane_t<V> comp_plane(0.f, 1.f, 0.f, 0.f); // y=0 plane
        assert(vec3_approx_equal(comp_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(comp_plane.Offset(), 0.f));

        // Normal and distance constructor
        fvec3 normal(0.f, 1.f, 0.f);
        fplane_t<V> norm_dist_plane(normal, -1.f);
        assert(vec3_approx_equal(norm_dist_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(norm_dist_plane.Offset(), 1.f));

        // Normal and point constructor
        fvec3 point(0.f, 2.f, 0.f);
        fplane_t<V> norm_point_plane(normal, point);
        assert(vec3_approx_equal(norm_point_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(norm_point_plane.Offset(), -2.f));

        // Three points constructor (right-handed system: points in CCW order)
        fvec3 p1(0.f, 0.f, 0.f), p2(1.f, 0.f, 0.f), p3(0.f, 0.f, 1.f);
        fplane_t<V> three_points_plane(p1, p2, p3); // Normal should be (0, 1, 0)
        assert(vec3_approx_equal(three_points_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(three_points_plane.Offset(), 0.f));

        // SIMD constructor (SIMD only)
        if constexpr (V)
        {
            floatx4 reg = _mm_set_ps(0.f, 0.f, 1.f, 0.f); // (x, y, z, d)
            fplane_t<V> simd_plane(reg);
            assert(vec3_approx_equal(simd_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
            assert(approx_equal(simd_plane.Offset(), 0.f));
        }

        // Cross-SIMD constructor
        fplane_t<!V> other_plane(0.f, 1.f, 0.f, -1.f); // y = -1
        fplane_t<V> cross_plane(other_plane);
        assert(vec3_approx_equal(cross_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(cross_plane.Offset(), 1.f));

        // Array of floats constructor
        std::array<float, 4> arr = { 0.f, 1.f, 0.f, -1.f };
        fplane_t<V> arr_plane(arr);
        assert(vec3_approx_equal(arr_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(arr_plane.Offset(), -1.f));
    }

    // Tests for static properties
    template <bool V>
    void test_static_properties()
    {
        // fromZero
        auto zero_plane = fplane_t<V>::fromZero();
        assert(vec3_approx_equal(zero_plane.Normal(), fvec3(0.f)));
        assert(approx_equal(zero_plane.Offset(), 0.f));
        // fromIdentity
        auto id_plane = fplane_t<V>::fromIdentity();
        assert(vec3_approx_equal(id_plane.Normal(), fvec3(0.f, 0.f, 1.f)));
        assert(approx_equal(id_plane.Offset(), 0.f));

        // fromNormalDistance
        auto norm_dist_plane = fplane_t<V>::fromNormalDistance(fvec3(0.f, 1.f, 0.f), -1.f);
        assert(vec3_approx_equal(norm_dist_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(norm_dist_plane.Offset(), 1.f)); // Expect m_D = -(-1.f) = 1.f

        // fromNormalPoint
        fvec3 point(0.f, 2.f, 0.f);
        auto norm_point_plane = fplane_t<V>::fromNormalPoint(fvec3(0.f, 1.f, 0.f), point);
        assert(vec3_approx_equal(norm_point_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(norm_point_plane.Offset(), -2.f));

        // fromThreePoints
        fvec3 p1(0.f, 0.f, 0.f), p2(1.f, 0.f, 0.f), p3(0.f, 0.f, 1.f);
        auto three_points_plane = fplane_t<V>::fromThreePoints(p1, p2, p3);
        assert(vec3_approx_equal(three_points_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(three_points_plane.Offset(), 0.f));
    }

    // Tests for static methods
    template <bool V>
    void test_static_methods()
    {
        // IntersectThreePlanes
        fplane_t<V> p1(fvec3(0.f, 1.f, 0.f), 0.f); // y=0
        fplane_t<V> p2(fvec3(1.f, 0.f, 0.f), 0.f); // x=0
        fplane_t<V> p3(fvec3(0.f, 0.f, 1.f), 0.f); // z=0
        fvec3 result;
        assert(fplane_t<V>::IntersectThreePlanes(p1, p2, p3, result));
        assert(vec3_approx_equal(result, fvec3(0.f)));

        // Dot
        fplane_t<V> plane(fvec3(0.f, 1.f, 0.f), -1.f); // y=1
        fvec3 point(0.f, 2.f, 0.f);
        assert(approx_equal(fplane_t<V>::Dot(plane, point), 3.f));

        // Distance
        assert(approx_equal(fplane_t<V>::Distance(plane, point), 3.f));

        // Side
        assert(fplane_t<V>::Side(plane, point) == 1); // Above
        assert(fplane_t<V>::Side(plane, fvec3(0.f, 0.f, 0.f)) == 1); // Above
        assert(fplane_t<V>::Side(plane, fvec3(0.f, -1.f, 0.f)) == 0); // On plane
        assert(fplane_t<V>::Side(plane, fvec3(0.f, -3.f, 0.f)) == -1); // Below

        // IntersectLine
        float t;
        fvec3 start(0.f, 0.f, 0.f), dir(0.f, 1.f, 0.f);
        assert(fplane_t<V>::IntersectLine(plane, t, start, dir));
        assert(approx_equal(t, -1.f, SMALL_EPSILON));

        // IntersectLineSegment
        fvec3 p0(0.f, 0.f, 0.f), p10(0.f, 2.f, 0.f);
        assert(!fplane_t<V>::IntersectLineSegment(plane, t, p0, p10));

        {
            fvec3 p0(0.f, -2.f, 0.f), p10(0.f, 0.f, 0.f);
            assert(fplane_t<V>::IntersectLineSegment(plane, t, p0, p10)); // Intersects at y = -1
            assert(approx_equal(t, 0.5f, SMALL_EPSILON)); // T = 0.5 -> (0, -1, 0)
        }
    }

    // Tests for instance methods
    template <bool V>
    void test_instance_methods()
    {
        fplane_t<V> plane(fvec3(0.f, 1.f, 0.f), -1.f); // y = -1
        // Normal
        assert(vec3_approx_equal(plane.Normal(), fvec3(0.f, 1.f, 0.f)));

        // Offset
        assert(approx_equal(plane.Offset(), 1.f));

        // GetOrigin
        assert(vec3_approx_equal(plane.GetOrigin(), fvec3(0.f, -1.f, 0.f)));

        // NormalLength
        assert(approx_equal(plane.NormalLength(), 1.f));

        // NormalLengthSq
        assert(approx_equal(plane.NormalLengthSq(), 1.f));

        // isFinite
        assert(plane.isFinite());
        fplane_t<V> inf_plane(std::numeric_limits<float>::infinity(), 0.f, 0.f, 0.f);
        assert(!inf_plane.isFinite());

        // isNormalized
        assert(plane.isNormalized());
        fplane_t<V> unnorm_plane(0.f, 2.f, 0.f, -1.f);
        assert(!unnorm_plane.isNormalized());

        // NormalizeCopy
        {
            fplane_t<V> unnorm_plane(0.f, 2.f, 0.f, -1.f); // y = -0.5
            auto norm_plane = unnorm_plane.NormalizeCopy();
            assert(vec3_approx_equal(norm_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
            assert(approx_equal(norm_plane.Offset(), 0.5f));
        }

        // Normalize
        unnorm_plane.Normalize();
        assert(vec3_approx_equal(unnorm_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(unnorm_plane.Offset(), 0.5f));

        // isNearlyZero
        assert(!plane.isNearlyZero());
        assert(fplane_t<V>::fromZero().isNearlyZero());

        // Equals
        assert(plane.Equals(fplane_t<V>(0.f, 1.f, 0.f, -1.f)));
        assert(!plane.Equals(fplane_t<V>(0.f, 1.f, 0.f, 0.f)));

        // Dot
        assert(approx_equal(plane.Dot(fvec3(0.f, 2.f, 0.f)), 3.f));
        // Distance
        assert(approx_equal(plane.Distance(fvec3(0.f, 2.f, 0.f)), 3.f));
        // Side
        assert(plane.Side(fvec3(0.f, 2.f, 0.f)) == 1);
        // IntersectLine
        float t;
        assert(plane.IntersectLine(t, fvec3(0.f, 0.f, 0.f), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(t, -1.f, SMALL_EPSILON));
        // IntersectLineSegment
        assert(!plane.IntersectLineSegment(t, fvec3(0.f, 0.f, 0.f), fvec3(0.f, 2.f, 0.f)));
        assert(plane.IntersectLineSegment(t, fvec3(0.f, -2.f, 0.f), fvec3(0.f, 0.f, 0.f)));
        assert(approx_equal(t, 0.5f, SMALL_EPSILON));
        // Project
        assert(vec3_approx_equal(plane.Project(fvec3(1.f, 2.f, 3.f)), fvec3(1.f, -1.f, 3.f)));
        // IsPointOnPlane
        assert(plane.IsPointOnPlane(fvec3(0.f, -1.f, 0.f)));
        assert(!plane.IsPointOnPlane(fvec3(0.f, 2.f, 0.f)));
        // SameSide
        assert(plane.SameSide(fvec3(0.f, 2.f, 0.f), fvec3(0.f, 3.f, 0.f)));
        assert(!plane.SameSide(fvec3(0.f, 2.f, 0.f), fvec3(0.f, -2.f, 0.f)));
        // setupFromNormalDistance
        plane.setupFromNormalDistance(fvec3(1.f, 0.f, 0.f), -2.f);
        assert(vec3_approx_equal(plane.Normal(), fvec3(1.f, 0.f, 0.f)));
        assert(approx_equal(plane.Offset(), 2.f));
        // setupFromNormalPoint
        plane.setupFromNormalPoint(fvec3(0.f, 1.f, 0.f), fvec3(0.f, 2.f, 0.f));
        assert(vec3_approx_equal(plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(plane.Offset(), -2.f));
        // setupFromThreePoints
        plane.setupFromThreePoints(fvec3(0.f, 0.f, 0.f), fvec3(1.f, 0.f, 0.f), fvec3(0.f, 0.f, 1.f));
        assert(vec3_approx_equal(plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(plane.Offset(), 0.f));
        // Translate
        auto trans_plane = plane.TranslateCopy(fvec3(0.f, 1.f, 0.f));
        assert(vec3_approx_equal(trans_plane.Normal(), fvec3(0.f, 1.f, 0.f)));
        assert(approx_equal(trans_plane.Offset(), -1.f));
        // Flip
        auto flip_plane = plane.FlipCopy();
        assert(vec3_approx_equal(flip_plane.Normal(), fvec3(0.f, -1.f, 0.f)));
        assert(approx_equal(flip_plane.Offset(), 0.f));
    }
    // Tests for operator overloads
    template <bool V>
    void test_operators()
    {
        fplane_t<V> plane(fvec3(0.f, 1.f, 0.f), -1.f);

        // operator-
        auto neg_plane = -plane;
        assert(vec3_approx_equal(neg_plane.Normal(), fvec3(0.f, -1.f, 0.f)));
        assert(approx_equal(neg_plane.Offset(), -1.f));

        // operator[]
        assert(approx_equal(plane[1], 1.f)); // Normal.y
        assert(approx_equal(plane[3], 1.f)); // Offset
        plane[3] = -2.f;
        assert(approx_equal(plane.Offset(), -2.f));

        // operator==
        assert(!(plane == fplane_t<V>(fvec3(0.f, 1.f, 0.f), -2.f)));
        assert((plane == fplane_t<V>(0.f, 1.f, 0.f, 2.f)));

        // operator!=
        assert((plane != fplane_t<V>(fvec3(0.f, 1.f, 0.f), -2.f)));
        assert(!(plane != fplane_t<V>(0.f, 1.f, 0.f, 2.f)));
    }

    // Tests for conversion operators
    template <bool V>
    void test_conversions()
    {
        fplane_t<V> plane(fvec3(0.f, 1.f, 0.f), 1.f);

        // operator std::array<float, 4>
        std::array<float, 4> arr = plane;
        assert(approx_equal(arr[0], 0.f));
        assert(approx_equal(arr[1], 1.f));
        assert(approx_equal(arr[2], 0.f));
        assert(approx_equal(arr[3], -1.f));

        // operator std::string
        std::string str = plane;
        assert(str == plane.ToString());

        // ToString
        assert(plane.ToString().find("(0, 1, 0,") != std::string::npos);
        assert(plane.ToString().find("-1)") != std::string::npos);
    }

    // Runs all unit tests for fplane_t
    template <bool V>
    void run_all_tests()
    {
        test_constructors<V>();
        test_static_properties<V>();
        test_static_methods<V>();
        test_instance_methods<V>();
        test_operators<V>();
        test_conversions<V>();
    }

    // Entry point for running all tests
    int RunTests()
    {
        run_all_tests<true>();
        run_all_tests<false>();
        return 0;
    }
} // namespace xmath::unit_test::_fplane