// Include the irect implementations here (paste the class definition and all function implementations from previous)
namespace xmath::unit_test::_irect
{
    int RunTests()
    {
        // Test fromZero
        auto r0 = irect::fromZero();
        assert(r0.m_Min == ivec2(0, 0));
        assert(r0.m_Max == ivec2(0, 0));

        // Test fromIdentity
        auto r1 = irect::fromIdentity();
        assert(r1.m_Min == ivec2(0, 0));
        assert(r1.m_Max == ivec2(1, 1));

        // Test fromPoint
        auto r2 = irect::fromPoint(ivec2(5, 10));
        assert(r2.m_Min == ivec2(5, 10));
        assert(r2.m_Max == ivec2(5, 10));

        // Test fromMinMax
        auto r3 = irect::fromMinMax(ivec2(1, 2), ivec2(3, 4));
        assert(r3.m_Min == ivec2(1, 2));
        assert(r3.m_Max == ivec2(3, 4));

        // Test fromCenterRadius
        auto r4 = irect::fromCenterRadius(ivec2(0, 0), 5);
        assert(r4.m_Min == ivec2(-5, -5));
        assert(r4.m_Max == ivec2(5, 5));

        // Test fromVerts
        std::array<ivec2, 3> verts = { ivec2(0,0), ivec2(10,20), ivec2(-5,15) };
        auto r5 = irect::fromVerts(verts);
        assert(r5.m_Min == ivec2(-5, 0));
        assert(r5.m_Max == ivec2(10, 20));

        // Test Intersect static rect-rect
        assert(irect::Intersect(r3, r4));
        assert(!irect::Intersect(r3, irect::fromMinMax(ivec2(4, 5), ivec2(6, 7))));

        // Test Intersect static rect-point
        assert(irect::Intersect(r3, ivec2(2, 3)));
        assert(!irect::Intersect(r3, ivec2(0, 0)));

        // Test Intersect static rect-line
        float t;
        assert(irect::Intersect(r3, t, ivec2(0, 0), ivec2(4, 4)));
        assert(t >= 0.0f && t <= 1.0f);
        assert(!irect::Intersect(r3, t, ivec2(4, 5), ivec2(6, 7)));

        // Test Intersect static rect-verts
        assert(irect::Intersect(r3, verts));
        std::array<ivec2, 1> no_verts = { ivec2(0,0) };
        assert(!irect::Intersect(r3, no_verts));

        // Test Intersection static
        auto inter = irect::Intersection(r3, r4);
        assert(inter.m_Min == ivec2(1, 2));
        assert(inter.m_Max == ivec2(3, 4));

        // Test constructors
        irect rc1(ivec2(1, 1));
        assert(rc1.m_Min == ivec2(1, 1));
        irect rc2(ivec2(0, 0), ivec2(10, 10));
        assert(rc2.m_Max == ivec2(10, 10));
        irect rc3(verts);
        assert(rc3.m_Min == ivec2(-5, 0));
        irect rc4(ivec2(2, 2), 3);
        assert(rc4.m_Min == ivec2(-1, -1));
        std::array<int, 4> arr = { 1,2,3,4 };
        irect rc5(arr);
        assert(rc5.m_Min.m_X == 1);
        irect rc6(1, 2, 3, 4);
        assert(rc6.m_Max.m_Y == 4);

        // Test conversions
        auto arr_out = static_cast<std::array<int, 4>>(rc6);
        assert(arr_out[0] == 1);
        std::string str = rc6;
        assert(!str.empty());
        std::ostringstream oss;
        oss << rc6;
        assert(oss.str() == "[(1, 2), (3, 4)]");  // Adjust based on ivec2 << impl

        // Test member Intersect
        assert(rc2.Intersect(rc4));
        assert(rc2.Intersect(ivec2(0, 0)));
        assert(rc2.Intersect(t, ivec2(-1, -1), ivec2(11, 11)));
        assert(rc2.Intersect(verts));
        auto inter_m = rc2.Intersection(rc4);
        assert(inter_m.m_Min == ivec2(0, 0));

        // Test isFinite (always true)
        assert(rc2.isFinite());

        // Test isValid
        assert(rc2.isValid());
        irect invalid(ivec2(5, 5), ivec2(0, 0));
        assert(!invalid.isValid());

        // Test getSize
        assert(rc2.getSize() == ivec2(10, 10));

        // Test getCenter
        assert(rc2.getCenter() == ivec2(5, 5));

        // Test getRadius
        assert(std::abs(rc2.getRadius() - 7.0710678f) < 1e-5f);

        // Test getRadiusSquared
        assert(rc2.getRadiusSquared() == 50);  // 200/4=50

        // Test getArea
        assert(rc2.getArea() == 100);

        // Test isInRange
        assert(rc2.isInRange(0, 10));
        assert(!rc2.isInRange(1, 9));

        // Test Contains
        assert(!rc2.Contains(rc4));
        assert(!rc4.Contains(rc2));

        // Test ContainsPoint
        assert(rc2.ContainsPoint(ivec2(5, 5)));

        // Test getClosestPoint
        assert(rc2.getClosestPoint(ivec2(-1, -1)) == ivec2(0, 0));
        assert(rc2.getClosestPoint(ivec2(15, 15)) == ivec2(10, 10));

        // Test getVerts
        std::array<ivec2, 4> dst;
        rc2.getVerts(dst);
        assert(dst[0] == ivec2(0, 0));
        assert(dst[1] == ivec2(10, 0));
        assert(dst[2] == ivec2(10, 10));
        assert(dst[3] == ivec2(0, 10));

        // Test getWidth
        assert(rc2.getWidth() == 10);

        // Test getHeight
        assert(rc2.getHeight() == 10);

        // Test Difference
        assert(rc2.Difference(rc4) > 0.0f);

        // Test isEmpty
        assert(!rc2.isEmpty());
        assert(r0.isEmpty());

        // Test setupFromPoint
        irect rs;
        rs.setupFromPoint(ivec2(3, 4));
        assert(rs.m_Min == ivec2(3, 4));

        // Test setupFromMinMax
        rs.setupFromMinMax(ivec2(0, 0), ivec2(5, 5));
        assert(rs.m_Max == ivec2(5, 5));

        // Test setupFromCenterRadius
        rs.setupFromCenterRadius(ivec2(0, 0), 2);
        assert(rs.m_Min == ivec2(-2, -2));

        // Test setupFromVerts
        rs.setupFromVerts(verts);
        assert(rs.m_Min == ivec2(-5, 0));

        // Test setup
        rs.setup(1, 2, 3, 4);
        assert(rs.m_Min == ivec2(1, 2));

        // Test setZero
        rs.setZero();
        assert(rs.m_Min == ivec2(0, 0));

        // Test setIdentity
        rs.setIdentity();
        assert(rs.m_Max == ivec2(1, 1));

        // Test AddVerts
        rs.AddVerts(verts);
        assert(rs.m_Min.m_X <= -5);

        // Test Inflate
        rs = rc2;
        rs.Inflate(ivec2(1, 1));
        assert(rs.m_Min == ivec2(-1, -1));
        assert(rs.m_Max == ivec2(11, 11));

        // Test Deflate
        rs.Deflate(ivec2(1, 1));
        assert(rs.m_Min == ivec2(0, 0));
        assert(rs.m_Max == ivec2(10, 10));

        // Test Translate
        rs.Translate(ivec2(5, 5));
        assert(rs.m_Min == ivec2(5, 5));

        // Test setWidth
        rs.setWidth(20);
        assert(rs.getWidth() == 20);

        // Test setHeight
        rs.setHeight(30);
        assert(rs.getHeight() == 30);

        // Test setSize
        rs.setSize(40, 50);
        assert(rs.getSize() == ivec2(40, 50));

        // Test setMax
        rs.setMax();
        assert(rs.m_Min.m_X == std::numeric_limits<int>::min());
        assert(rs.m_Max.m_X == std::numeric_limits<int>::max());

        // Test operator+= rect
        rs = rc2;
        rs += rc4;
        assert(rs.m_Min == ivec2(-1, -1));

        // Test operator+= point
        rs += ivec2(-2, -2);
        assert(rs.m_Min == ivec2(-2, -2));

        // Test operator+ rect
        auto ru = rc2 + rc4;
        assert(ru.m_Min == ivec2(-1, -1));

        // Test operator+ point
        auto rp = rc2 + ivec2(1, 1);
        assert(rp.m_Max == ivec2(10, 10));

        // Test operator==
        assert(rc2 == rc2);
        assert(!(rc2 == rc4));

        // Test operator!=
        assert(rc2 != rc4);

        // Test operator[] const
        assert(rc2[0] == 0);
        assert(rc2[2] == 10);

        // Test operator[]
        irect& ref = rc2;
        ref[0] = 2;
        assert(ref.m_Min.m_X == 2);

        // Test friend operator+
        auto rf = ivec2(11, 11) + rc2;
        assert(rf.m_Max == ivec2(11, 11));

        return 0;
    }
}
