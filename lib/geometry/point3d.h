#ifndef POINT3D_H
#define POINT3D_H

#include "../core/common.h"
#include "../core/units.h"

#include <iostream>

using namespace std;
using namespace boost::units::photonflow::literals;

class Point3D {
public:
    // Point Public Methods
    Point3D();
    Point3D(boost::units::photonflow::Length xx, boost::units::photonflow::Length yy, boost::units::photonflow::Length zz);
#ifndef NDEBUG
    Point3D(const Point3D &p);
    Point3D &operator=(const Point3D &p);
#endif // !NDEBUG
    Point3D operator+(const Length3D &v) const;

    Point3D &operator+=(const Length3D &v);
    Length3D operator-(const Point3D &p) const;
    Point3D operator-(const Length3D &v) const;
    Point3D &operator-=(const Length3D &v);
    Point3D &operator+=(const Point3D &p);
    Point3D operator+(const Point3D &p) const;
    Point3D operator* (double f) const;
    Point3D &operator*=(double f);
    Point3D operator/ (double f) const;
    Point3D &operator/=(double f);
    boost::units::photonflow::Length operator[](int i) const;

    boost::units::photonflow::Length &operator[](int i);
    bool hasNaNs() const;

    bool operator==(const Point3D &p) const;
    bool operator!=(const Point3D &p) const;
    friend std::ostream& operator<< (std::ostream &out, const Point3D &point);

    // Point Public Data
    boost::units::photonflow::Length x, y, z;
};

#include "vector3d.h"

inline Point3D::Point3D()
    : x(0.0_m)
    , y(0.0_m)
    , z(0.0_m)
{

}

inline Point3D::Point3D(boost::units::photonflow::Length xx, boost::units::photonflow::Length yy, boost::units::photonflow::Length zz)
    : x(xx),
      y(yy),
      z(zz)
{
    photonflowAssert(!hasNaNs());
}

#ifndef NDEBUG
inline Point3D::Point3D(const Point3D &p) {
    photonflowAssert(!p.hasNaNs());
    x = p.x; y = p.y; z = p.z;
}
inline Point3D &Point3D::operator=(const Point3D &p) {
    photonflowAssert(!p.hasNaNs());
    x = p.x; y = p.y; z = p.z;
    return *this;
}
#endif

inline Point3D Point3D::operator+(const Length3D &v) const {
    photonflowAssert(!v.hasNaNs());
    return Point3D(x + v.x, y + v.y, z + v.z);
}

inline Point3D &Point3D::operator+=(const Length3D &v) {
    photonflowAssert(!v.hasNaNs());
    x += v.x; y += v.y; z += v.z;
    return *this;
}

inline Length3D Point3D::operator-(const Point3D &p) const {
    photonflowAssert(!p.hasNaNs());
    return Length3D(x - p.x, y - p.y, z - p.z);
}

inline Point3D Point3D::operator-(const Length3D &v) const {
    photonflowAssert(!v.hasNaNs());
    return Point3D(x - v.x, y - v.y, z - v.z);
}

inline Point3D &Point3D::operator-=(const Length3D &v) {
    photonflowAssert(!v.hasNaNs());
    x -= v.x; y -= v.y; z -= v.z;
    return *this;
}

inline Point3D &Point3D::operator+=(const Point3D &p) {
    photonflowAssert(!p.hasNaNs());
    x += p.x; y += p.y; z += p.z;
    return *this;
}

inline Point3D Point3D::operator+(const Point3D &p) const {
    photonflowAssert(!p.hasNaNs());
    return Point3D(x + p.x, y + p.y, z + p.z);
}

inline Point3D Point3D::operator*(double f) const {
    return Point3D(f*x, f*y, f*z);
}

inline Point3D &Point3D::operator*=(double f) {
    x *= f; y *= f; z *= f;
    return *this;
}

inline Point3D Point3D::operator/(double f) const {
    double inv = 1.f/f;
    return Point3D(inv*x, inv*y, inv*z);
}

inline Point3D &Point3D::operator/=(double f) {
    double inv = 1.f/f;
    x *= inv; y *= inv; z *= inv;
    return *this;
}

inline boost::units::photonflow::Length Point3D::operator[](int i) const {
    photonflowAssert(i >= 0 && i <= 2);
    return (&x)[i];
}

inline boost::units::photonflow::Length &Point3D::operator[](int i) {
    photonflowAssert(i >= 0 && i <= 2);
    return (&x)[i];
}

inline bool Point3D::hasNaNs() const {
    return isnan(x) || isnan(y) || isnan(z);
}

inline bool Point3D::operator==(const Point3D &p) const {
    return x == p.x && y == p.y && z == p.z;
}

inline bool Point3D::operator!=(const Point3D &p) const {
    return x != p.x || y != p.y || z != p.z;
}

#endif // POINT3D_H
