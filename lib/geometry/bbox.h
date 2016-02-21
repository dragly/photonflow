#ifndef BBOX_H
#define BBOX_H

#include "../core/units.h"
#include "point3d.h"

using namespace boost::units::photonflow::literals;

class Ray;

class BBox {
public:
    // BBox Public Methods
    BBox() {
        pMin = Point3D( INFINITY * boost::units::photonflow::micrometer,  INFINITY * boost::units::photonflow::micrometer,  INFINITY * boost::units::photonflow::micrometer);
        pMax = Point3D(-INFINITY * boost::units::photonflow::micrometer, -INFINITY * boost::units::photonflow::micrometer, -INFINITY * boost::units::photonflow::micrometer);
    }
    BBox(const Point3D &p) : pMin(p), pMax(p) { }
    BBox(const Point3D &p1, const Point3D &p2) {
        pMin = Point3D(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
        pMax = Point3D(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
    }
    friend BBox makeUnion(const BBox &b, const Point3D &p);
    friend BBox makeUnion(const BBox &b, const BBox &b2);
    bool overlaps(const BBox &b) const {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }
    bool inside(const Point3D &pt) const {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }
    bool fuzzyInside(const Point3D &pt) const {
        const boost::units::photonflow::Length eps = 0.01_um;
        return (pt.x >= pMin.x - eps && pt.x <= pMax.x + eps &&
                pt.y >= pMin.y - eps && pt.y <= pMax.y + eps &&
                pt.z >= pMin.z - eps && pt.z <= pMax.z + eps);
    }
    void expand(boost::units::photonflow::Length delta) {
        pMin -= Length3D(delta, delta, delta);
        pMax += Length3D(delta, delta, delta);
    }
    auto surfaceArea() const {
        Length3D d = pMax - pMin;
        return 2.0 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    auto volume() const {
        Length3D d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int maximumExtent() const {
        Length3D diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }
    const Point3D &operator[](int i) const;
    Point3D &operator[](int i);
    Point3D lerp(double tx, double ty, double tz) const {
        return Point3D(::lerp(tx, pMin.x, pMax.x), ::lerp(ty, pMin.y, pMax.y),
                       ::lerp(tz, pMin.z, pMax.z));
    }
    Length3D offset(const Point3D &p) const {
        return Length3D((p.x - pMin.x) / (pMax.x - pMin.x) * 1.0_um,
                        (p.y - pMin.y) / (pMax.y - pMin.y) * 1.0_um,
                        (p.z - pMin.z) / (pMax.z - pMin.z) * 1.0_um);
    }
    void boundingSphere(Point3D *c, boost::units::photonflow::Length *rad) const;
    bool intersectP(const Ray &ray, double *hitt0 = NULL, double *hitt1 = NULL) const;

    bool operator==(const BBox &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const BBox &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }

    // BBox Public Data
    Point3D pMin, pMax;
};

#endif // BBOX_H
