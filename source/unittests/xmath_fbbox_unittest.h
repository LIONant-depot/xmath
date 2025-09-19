#include <cassert>
#include <cmath>
#include <array>

namespace xmath::unit_test::_fbbox
{
    // Helper functions (reused from fmat4 unit tests)
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

    bool bbox_approx_equal(const std::array<float, 6>& a, const std::array<float, 6>& b, float eps = EPSILON)
    {
        for (size_t i = 0; i < 6; ++i)
        {
            if (!approx_equal(a[i], b[i], eps)) return false;
        }
        return true;
    }

    //------------------------------------------------------------------------------
    // Tests for constructors
    //------------------------------------------------------------------------------
    template <bool V>
    void test_constructors()
    {
        // Default constructor (uninitialized)
        fbbox_t<V> default_bbox;
        assert( default_bbox.isValid() == true );

        // No specific checks for default, as memory is uninitialized

#if 0
        // Single point constructor
        fvec3 point(1.f, 2.f, 3.f);
        fbbox_t<V> point_bbox(point);
        assert(vec3_approx_equal(point_bbox.m_Min, point));
        assert(vec3_approx_equal(point_bbox.m_Max, point));
#endif

        // Min/max constructor
        fvec3 min(0.f, 0.f, 0.f), max(2.f, 2.f, 2.f);
        fbbox_t<V> minmax_bbox(min, max);
        assert(vec3_approx_equal(minmax_bbox.m_Min, min));
        assert(vec3_approx_equal(minmax_bbox.m_Max, max));

        // Array of vertices constructor
        std::array verts = { fvec3(0.f, 0.f, 0.f), fvec3(1.f, 2.f, 3.f), fvec3(-1.f, -2.f, -3.f) };
        fbbox_t<V> verts_bbox(verts);
        assert(vec3_approx_equal(verts_bbox.m_Min, fvec3(-1.f, -2.f, -3.f)));
        assert(vec3_approx_equal(verts_bbox.m_Max, fvec3(1.f, 2.f, 3.f)));

        // Center and radius constructor
        fbbox_t<V> center_bbox(fvec3(0.f, 0.f, 0.f), 1.f);
        assert(vec3_approx_equal(center_bbox.m_Min, fvec3(-1.f, -1.f, -1.f)));
        assert(vec3_approx_equal(center_bbox.m_Max, fvec3(1.f, 1.f, 1.f)));

        // SIMD min/max constructor (SIMD only)
        if constexpr (V)
        {
            std::array minmax = { _mm_set_ps(1.f, 0.f, 0.f, 0.f), _mm_set_ps(1.f, 2.f, 2.f, 2.f) };
            fbbox_t<V> simd_bbox(minmax);
            assert(vec3_approx_equal(simd_bbox.m_Min, fvec3(0.f, 0.f, 0.f)));
            assert(vec3_approx_equal(simd_bbox.m_Max, fvec3(2.f, 2.f, 2.f)));
        }

        // Cross-SIMD constructor
        fbbox_t<!V> other_bbox(fvec3(1.f, 1.f, 1.f), fvec3(2.f, 2.f, 2.f));
        fbbox_t<V> cross_bbox(other_bbox);
        assert(vec3_approx_equal(cross_bbox.m_Min, fvec3(1.f, 1.f, 1.f)));
        assert(vec3_approx_equal(cross_bbox.m_Max, fvec3(2.f, 2.f, 2.f)));

        // Array of floats constructor
        std::array arr = { -1.f, -2.f, -3.f, 1.f, 2.f, 3.f };
        fbbox_t<V> arr_bbox(arr);
        assert(vec3_approx_equal(arr_bbox.m_Min, fvec3(-1.f, -2.f, -3.f)));
        assert(vec3_approx_equal(arr_bbox.m_Max, fvec3(1.f, 2.f, 3.f)));
    }

    //------------------------------------------------------------------------------
    // Tests for static properties
    //------------------------------------------------------------------------------
    template <bool V>
    void test_static_properties()
    {
        // fromZero
        auto zero_bbox = fbbox_t<V>::fromZero();
        assert(vec3_approx_equal(zero_bbox.m_Min, fvec3::fromZero()));
        assert(vec3_approx_equal(zero_bbox.m_Max, fvec3::fromZero()));

        // fromIdentity
        auto id_bbox = fbbox_t<V>::fromIdentity();
        assert(vec3_approx_equal(id_bbox.m_Min, fvec3(std::numeric_limits<float>::max())));
        assert(vec3_approx_equal(id_bbox.m_Max, fvec3(std::numeric_limits<float>::min())));

        // fromPoint
        fvec3 point(1.f, 2.f, 3.f);
        auto point_bbox = fbbox_t<V>::fromPoint(point);
        assert(vec3_approx_equal(point_bbox.m_Min, point));
        assert(vec3_approx_equal(point_bbox.m_Max, point));

        // fromMinMax
        fvec3 min(0.f, 0.f, 0.f), max(2.f, 2.f, 2.f);
        auto minmax_bbox = fbbox_t<V>::fromMinMax(min, max);
        assert(vec3_approx_equal(minmax_bbox.m_Min, min));
        assert(vec3_approx_equal(minmax_bbox.m_Max, max));

        // fromCenterRadius
        auto center_bbox = fbbox_t<V>::fromCenterRadius(fvec3(0.f), 1.f);
        assert(vec3_approx_equal(center_bbox.m_Min, fvec3(-1.f)));
        assert(vec3_approx_equal(center_bbox.m_Max, fvec3(1.f)));

        // fromVerts
        std::array verts = { fvec3(0.f), fvec3(1.f, 2.f, 3.f), fvec3(-1.f, -2.f, -3.f) };
        auto verts_bbox = fbbox_t<V>::fromVerts(std::span(verts));
        assert(vec3_approx_equal(verts_bbox.m_Min, fvec3(-1.f, -2.f, -3.f)));
        assert(vec3_approx_equal(verts_bbox.m_Max, fvec3(1.f, 2.f, 3.f)));
    }

    //------------------------------------------------------------------------------
    // Tests for static methods
    //------------------------------------------------------------------------------
    template <bool V>
    void test_static_methods()
    {
        // Intersect (boxes)
        fbbox_t<V> box1(fvec3(0.f), fvec3(1.f));
        fbbox_t<V> box2(fvec3(0.5f), fvec3(1.5f));
        fbbox_t<V> box3(fvec3(2.f), fvec3(3.f));
        assert(fbbox_t<V>::Intersect(box1, box2));
        assert(!fbbox_t<V>::Intersect(box1, box3));

        // Intersect (point)
        fvec3 point(0.5f);
        fvec3 outside(2.f);
        assert(fbbox_t<V>::Intersect(box1, point));
        assert(!fbbox_t<V>::Intersect(box1, outside));

        // Intersect (plane)
        fplane_t<V> plane(fvec3(0.f, 1.f, 0.f), 0.f); // y=0 plane
        fbbox_t<V> above(fvec3(0.f, 1.f, 0.f), fvec3(1.f, 2.f, 1.f));
        fbbox_t<V> straddling(fvec3(0.f, -1.f, 0.f), fvec3(1.f, 1.f, 1.f));
        assert(!fbbox_t<V>::Intersect(above, plane));
        assert(fbbox_t<V>::Intersect(straddling, plane));

        // Intersect (line segment)
        float t;
        fvec3 p0(0.f, 0.f, 0.f), p1(1.f, 1.f, 1.f);
        fvec3 p2(2.f, 2.f, 2.f), p3(3.f, 3.f, 3.f);
        assert(fbbox_t<V>::Intersect(box1, t, p0, p1));
        assert(approx_equal(t, 0.f, SMALL_EPSILON));
        assert(!fbbox_t<V>::Intersect(box1, t, p2, p3));

        // Intersect (vertices)
        std::array verts = { fvec3(0.5f), fvec3(2.f) };
        assert(fbbox_t<V>::Intersect(box1, std::span(verts) ));
        fvec3 verts_out[] = { fvec3(2.f), fvec3(3.f) };
        assert(!fbbox_t<V>::Intersect(box1, std::span(verts_out, verts_out + 2)));

        // IntersectTriangle
        fvec3 tri1[] = { fvec3(0.5f), fvec3(0.6f, 0.6f, 0.5f), fvec3(0.5f, 0.5f, 0.6f) };
        fvec3 tri2[] = { fvec3(2.f), fvec3(2.1f, 2.1f, 2.f), fvec3(2.f, 2.f, 2.1f) };
        assert(fbbox_t<V>::IntersectTriangle(box1, tri1[0], tri1[1], tri1[2]));
        assert(!fbbox_t<V>::IntersectTriangle(box1, tri2[0], tri2[1], tri2[2]));

        // Intersection
        auto inter_bbox = fbbox_t<V>::Intersection(box1, box2);
        assert(vec3_approx_equal(inter_bbox.m_Min, fvec3(0.5f)));
        assert(vec3_approx_equal(inter_bbox.m_Max, fvec3(1.f)));
    }

    //------------------------------------------------------------------------------
    // Tests for instance methods
    //------------------------------------------------------------------------------
    template <bool V>
    void test_instance_methods()
    {
        fbbox_t<V> box(fvec3(0.f), fvec3(1.f));

        // isFinite
        assert(box.isFinite());
        fbbox_t<V> inf_box(fvec3(std::numeric_limits<float>::infinity()), fvec3(1.f));
        assert(!inf_box.isFinite());

        // isValid
        assert(box.isValid());
        fbbox_t<V> invalid_box(fvec3(2.f), fvec3(1.f));
        assert(!invalid_box.isValid());

        // getSize
        assert(vec3_approx_equal(box.getSize(), fvec3(1.f)));

        // getCenter
        assert(vec3_approx_equal(box.getCenter(), fvec3(0.5f)));

        // getRadius
        assert(approx_equal(box.getRadius(), std::sqrt(3.f) / 2.f, SMALL_EPSILON));

        // getRadiusSquared
        assert(approx_equal(box.getRadiusSquared(), 3.f / 4.f, SMALL_EPSILON));

        // getSurfaceArea
        assert(approx_equal(box.getSurfaceArea(), 6.f, SMALL_EPSILON));

        // isInRange
        assert(box.isInRange(0.f, 1.f));
        assert(!box.isInRange(0.1f, 0.9f));

        // Contains
        fbbox_t<V> smaller(fvec3(0.1f), fvec3(0.9f));
        fbbox_t<V> larger(fvec3(-1.f), fvec3(2.f));
        assert(box.Contains(smaller));
        assert(!box.Contains(larger));

        // ContainsPoint
        fvec3 point(0.5f);
        fvec3 outside(2.f);
        assert(box.ContainsPoint(point));
        assert(!box.ContainsPoint(outside));

        // getClosestPoint
        fvec3 far_point(2.f, 2.f, 2.f);
        assert(vec3_approx_equal(box.getClosestPoint(far_point), fvec3(1.f)));

        // getVerts
        std::array<fvec3,8> verts;
        box.getVerts(verts);
        assert(vec3_approx_equal(verts[0], fvec3(0.f)));
        assert(vec3_approx_equal(verts[7], fvec3(1.f)));
        assert(vec3_approx_equal(verts[3], fvec3(0.f, 0.f, 1.f)));

        // setupFromPoint
        box.setupFromPoint(fvec3(1.f));
        assert(vec3_approx_equal(box.m_Min, fvec3(1.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(1.f)));

        // setupFromMinMax
        box.setupFromMinMax(fvec3(0.f), fvec3(2.f));
        assert(vec3_approx_equal(box.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(2.f)));

        // setupFromCenterRadius
        box.setupFromCenterRadius(fvec3(0.f), 1.f);
        assert(vec3_approx_equal(box.m_Min, fvec3(-1.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(1.f)));

        // setupFromVerts
        std::array setup_verts = { fvec3(0.f), fvec3(1.f, 2.f, 3.f) };
        box.setupFromVerts(setup_verts);
        assert(vec3_approx_equal(box.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(1.f, 2.f, 3.f)));

        // setZero
        box.setZero();
        assert(vec3_approx_equal(box.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(0.f)));

        // setIdentity
        box.setIdentity();
        assert(vec3_approx_equal(box.m_Min, fvec3(-1.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(1.f)));

        // Inflate
        box.Inflate(fvec3(0.5f));
        assert(vec3_approx_equal(box.m_Min, fvec3(-1.5f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(1.5f)));

        // Deflate
        box.Deflate(fvec3(0.5f));
        assert(vec3_approx_equal(box.m_Min, fvec3(-1.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(1.f)));

        // Translate
        box.Translate(fvec3(1.f));
        assert(vec3_approx_equal(box.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(2.f)));

        // Transform and TransformCopy
        fmat4 trans_mat = fmat4::fromTranslation(fvec3(1.f, 2.f, 3.f));
        auto Transformed = box.TransformCopy(trans_mat);
        assert(vec3_approx_equal(Transformed.m_Min, fvec3(1.f, 2.f, 3.f)));
        assert(vec3_approx_equal(Transformed.m_Max, fvec3(3.f, 4.f, 5.f)));
        box.Transform(trans_mat);
        assert(vec3_approx_equal(box.m_Min, fvec3(1.f, 2.f, 3.f)));
        assert(vec3_approx_equal(box.m_Max, fvec3(3.f, 4.f, 5.f)));
    }

    //------------------------------------------------------------------------------
    // Tests for instance methods (Intersection-related)
    //------------------------------------------------------------------------------
    template <bool V>
    void test_instance_Intersection_methods()
    {
        fbbox_t<V> box(fvec3(0.f), fvec3(1.f));

        // Intersect (box)
        fbbox_t<V> box2(fvec3(0.5f), fvec3(1.5f));
        fbbox_t<V> box3(fvec3(2.f), fvec3(3.f));
        assert(box.Intersect(box2));
        assert(!box.Intersect(box3));

        // Intersect (point)
        fvec3 point(0.5f);
        fvec3 outside(2.f);
        assert(box.Intersect(point));
        assert(!box.Intersect(outside));

        // Intersect (plane)
        fplane_t<V> plane(fvec3(0.f, 1.f, 0.f), 0.f);
        assert(box.Intersect(plane));
        fplane_t<V> plane_above(fvec3(0.f, 1.f, 0.f), 2.f);
        assert(!box.Intersect(plane_above));

        // Intersect (line segment)
        float t;
        fvec3 p0(0.f), p1(1.f);
        fvec3 p2(2.f), p3(3.f);
        assert(box.Intersect(t, p0, p1));
        assert(approx_equal(t, 0.f, SMALL_EPSILON));
        assert(!box.Intersect(t, p2, p3));

        // Intersect (vertices)
        std::array verts = { fvec3(0.5f), fvec3(2.f) };
        assert(box.Intersect(verts));
        std::array verts_out = { fvec3(2.f), fvec3(3.f) };
        assert(!box.Intersect(verts_out));

        // IntersectTriangle
        fvec3 tri1[] = { fvec3(0.5f), fvec3(0.6f, 0.6f, 0.5f), fvec3(0.5f, 0.5f, 0.6f) };
        fvec3 tri2[] = { fvec3(2.f), fvec3(2.1f, 2.1f, 2.f), fvec3(2.f, 2.f, 2.1f) };
        assert(box.IntersectTriangle(tri1[0], tri1[1], tri1[2]));
        assert(!box.IntersectTriangle(tri2[0], tri2[1], tri2[2]));

        // Intersection
        auto inter_bbox = box.Intersection(box2);
        assert(vec3_approx_equal(inter_bbox.m_Min, fvec3(0.5f)));
        assert(vec3_approx_equal(inter_bbox.m_Max, fvec3(1.f)));
    }

    //------------------------------------------------------------------------------
    // Tests for operator overloads
    //------------------------------------------------------------------------------
    template <bool V>
    void test_operators()
    {
        fbbox_t<V> box(fvec3(0.f), fvec3(1.f));

        // operator+= (box)
        fbbox_t<V> box2(fvec3(0.5f), fvec3(1.5f));
        auto copy = box;
        copy += box2;
        assert(vec3_approx_equal(copy.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(copy.m_Max, fvec3(1.5f)));

        // operator+= (point)
        fvec3 point(2.f);
        copy = box;
        copy += point;
        assert(vec3_approx_equal(copy.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(copy.m_Max, fvec3(2.f)));

        // operator+ (box)
        auto sum = box + box2;
        assert(vec3_approx_equal(sum.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(sum.m_Max, fvec3(1.5f)));

        // operator+ (point)
        sum = box + point;
        assert(vec3_approx_equal(sum.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(sum.m_Max, fvec3(2.f)));

        // operator+ (point + box)
        sum = point + box;
        assert(vec3_approx_equal(sum.m_Min, fvec3(0.f)));
        assert(vec3_approx_equal(sum.m_Max, fvec3(2.f)));

        // operator==
        assert(box == fbbox_t<V>(fvec3(0.f), fvec3(1.f)));
        assert(!(box == box2));

        // operator!=
        assert(box != box2);
        assert(!(box != fbbox_t<V>(fvec3(0.f), fvec3(1.f))));

        // operator[]
        std::array<float, 6> arr = box;
        assert(approx_equal(box[0], 0.f));
        assert(approx_equal(box[3], 1.f));
        box[0] = 2.f;
        assert(approx_equal(box.m_Min.m_X, 2.f));
    }

    //------------------------------------------------------------------------------
    // Tests for conversion operators
    //------------------------------------------------------------------------------
    template <bool V>
    void test_conversions()
    {
        fbbox_t<V> box(fvec3(0.f), fvec3(1.f));

        // operator std::array<float, 6>
        std::array<float, 6> arr = box;
        assert(bbox_approx_equal(arr, { 0.f, 0.f, 0.f, 1.f, 1.f, 1.f }));

        // operator std::string
        std::string str = box;
        assert(str == box.toString());

        // toString
        assert(box.toString().find("(0, 0, 0)") != std::string::npos);
        assert(box.toString().find("(1, 1, 1)") != std::string::npos);
    }

    //------------------------------------------------------------------------------
    // Runs all unit tests for fbbox_t
    //------------------------------------------------------------------------------
    template <bool V>
    void run_all_tests()
    {
        test_constructors<V>();
        test_static_properties<V>();
        test_static_methods<V>();
        test_instance_methods<V>();
        test_instance_Intersection_methods<V>();
        test_operators<V>();
        test_conversions<V>();
    }

    //------------------------------------------------------------------------------
    // Entry point for running all tests
    //------------------------------------------------------------------------------
    int RunTests()
    {
        run_all_tests<true>();
        run_all_tests<false>();
        return 0;
    }
} // namespace xmath::unit_test::_fbbox
