#ifndef BBOX_H
#define BBOX_H

#include "../core/units.h"
#include "point3d.h"

using namespace photonflow::literals;

class Ray;

class BBox {
public:
    // BBox Public Methods
    BBox();
    BBox(const Point3D &p);
    BBox(const Point3D &p1, const Point3D &p2);
    friend BBox makeUnion(const BBox &b, const Point3D &p);
    friend BBox makeUnion(const BBox &b, const BBox &b2);
    bool overlaps(const BBox &b) const;
    bool inside(const Point3D &pt) const;
    bool fuzzyInside(const Point3D &pt) const;
    void expand(photonflow::Length delta);
    photonflow::Area surfaceArea() const;
    photonflow::Volume volume() const;
    int maximumExtent() const;
    const Point3D &operator[](int i) const;
    Point3D &operator[](int i);
    Point3D lerp(double tx, double ty, double tz) const;
    Length3D offset(const Point3D &p) const;
    void boundingSphere(Point3D *c, photonflow::Length *rad) const;
    bool intersectP(const Ray &ray, double *hitt0 = NULL, double *hitt1 = NULL) const;

    bool operator==(const BBox &b) const;
    bool operator!=(const BBox &b) const;

    // BBox Public Data
    Point3D pMin, pMax;
};

inline BBox::BBox() {
    pMin = Point3D( INFINITY * photonflow::micrometer,  INFINITY * photonflow::micrometer,  INFINITY * photonflow::micrometer);
    pMax = Point3D(-INFINITY * photonflow::micrometer, -INFINITY * photonflow::micrometer, -INFINITY * photonflow::micrometer);
}

inline BBox::BBox(const Point3D &p) : pMin(p), pMax(p) { }

inline BBox::BBox(const Point3D &p1, const Point3D &p2) {
    pMin = Point3D(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
    pMax = Point3D(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
}

inline bool BBox::overlaps(const BBox &b) const {
    bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
    bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
    bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
    return (x && y && z);
}

inline bool BBox::inside(const Point3D &pt) const {
    return (pt.x >= pMin.x && pt.x <= pMax.x &&
            pt.y >= pMin.y && pt.y <= pMax.y &&
            pt.z >= pMin.z && pt.z <= pMax.z);
}

inline bool BBox::fuzzyInside(const Point3D &pt) const {
    const photonflow::Length eps = 0.01_um;
    return (pt.x >= pMin.x - eps && pt.x <= pMax.x + eps &&
            pt.y >= pMin.y - eps && pt.y <= pMax.y + eps &&
            pt.z >= pMin.z - eps && pt.z <= pMax.z + eps);
}

inline void BBox::expand(photonflow::Length delta) {
    pMin -= Length3D(delta, delta, delta);
    pMax += Length3D(delta, delta, delta);
}

inline photonflow::Area BBox::surfaceArea() const {
    Length3D d = pMax - pMin;
    return 2.0 * (d.x * d.y + d.x * d.z + d.y * d.z);
}

inline photonflow::Volume BBox::volume() const {
    Length3D d = pMax - pMin;
    return d.x * d.y * d.z;
}

inline int BBox::maximumExtent() const {
    Length3D diag = pMax - pMin;
    if (diag.x > diag.y && diag.x > diag.z)
        return 0;
    else if (diag.y > diag.z)
        return 1;
    else
        return 2;
}

inline Point3D BBox::lerp(double tx, double ty, double tz) const {
    return Point3D(::lerp(tx, pMin.x, pMax.x), ::lerp(ty, pMin.y, pMax.y),
                   ::lerp(tz, pMin.z, pMax.z));
}

inline Length3D BBox::offset(const Point3D &p) const {
    return Length3D((p.x - pMin.x) / (pMax.x - pMin.x) * 1.0_um,
                    (p.y - pMin.y) / (pMax.y - pMin.y) * 1.0_um,
                    (p.z - pMin.z) / (pMax.z - pMin.z) * 1.0_um);
}

inline bool BBox::operator==(const BBox &b) const {
    return b.pMin == pMin && b.pMax == pMax;
}

inline bool BBox::operator!=(const BBox &b) const {
    return b.pMin != pMin || b.pMax != pMax;
}

#endif // BBOX_H
