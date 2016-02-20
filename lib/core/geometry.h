
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_GEOMETRY_H
#define PBRT_CORE_GEOMETRY_H

#include <cmath>
#include <boost/units/cmath.hpp>
#include <ostream>
#include <type_traits>

#include "../core/common.h"
#include "../core/units.h"

//using namespace boost::units::photonflow;
using namespace boost::units::photonflow::literals;

class Point3D;

class Normal;

// Geometry Declarations
template<typename T = double>
class GeneralVector3D {
public:
    // Vector Public Methods
    //    GeneralVector3D() { x = y = z = 0.0; }
    GeneralVector3D() {} // TODO: Proper initialization to zero for all types
    GeneralVector3D(T xx, T yy, T zz)
        : x(xx), y(yy), z(zz) {
        photonFlowAssert(!hasNaNs());
    }
    //    bool hasNaNs() const { return isnan(x) || isnan(y) || isnan(z); }
    bool hasNaNs() const { return false; }
    explicit GeneralVector3D(const Point3D &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    GeneralVector3D(const GeneralVector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x = v.x; y = v.y; z = v.z;
    }

    GeneralVector3D &operator=(const GeneralVector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
#endif // !NDEBUG
    GeneralVector3D operator+(const GeneralVector3D &v) const {
        photonFlowAssert(!v.hasNaNs());
        return GeneralVector3D(x + v.x, y + v.y, z + v.z);
    }

    GeneralVector3D& operator+=(const GeneralVector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    GeneralVector3D operator-(const GeneralVector3D &v) const {
        photonFlowAssert(!v.hasNaNs());
        return GeneralVector3D(x - v.x, y - v.y, z - v.z);
    }

    GeneralVector3D& operator-=(const GeneralVector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    template<typename U>
    auto operator*(U f) const {
        return GeneralVector3D<decltype(f*x)>(f*x, f*y, f*z);
    }

    GeneralVector3D &operator*=(double f) {
        photonFlowAssert(!isnan(f));
        x *= f; y *= f; z *= f;
        return *this;
    }

    template<typename U>
    auto operator/(U f) const {
        photonFlowAssert(f != 0);
        auto inv = 1.0 / f;
        return (*this)*inv;
    }

    GeneralVector3D &operator/=(double f) {
        photonFlowAssert(f != 0);
        auto inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    GeneralVector3D operator-() const { return GeneralVector3D(-x, -y, -z); }
    T operator[](int i) const {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    T &operator[](int i)
    {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    auto lengthSquared() const
    {
        return x*x + y*y + z*z;
    }

    auto length() const
    {
        return sqrt(lengthSquared());
    }

    explicit GeneralVector3D(const Normal &n);

    bool operator==(const GeneralVector3D &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    bool operator!=(const GeneralVector3D &v) const {
        return x != v.x || y != v.y || z != v.z;
    }

    template<typename U>
    friend std::ostream& operator<< (std::ostream &out, const GeneralVector3D<U> &vector);

    GeneralVector3D perpendicular() const
    {
        if(x == 0.0 * T() && y == 0.0 * T()) {
            if(z == 0.0 * T()) {
                return GeneralVector3D(0.0 * T(), 0.0 * T(), 0.0 * T());
            }
            return GeneralVector3D(0.0 * T(), 1.0 * T(), 0.0 * T());
        }
        return GeneralVector3D(-y, x, 0.0 * T());
    }

    GeneralVector3D normalized() {
        return *this / length().value();
    }

    // Vector Public Data
    T x, y, z;
};
using Vector3D = GeneralVector3D<boost::units::photonflow::length>;

class Point3D {
public:
    // Point Public Methods
    Point3D() {  }
    Point3D(boost::units::photonflow::length xx, boost::units::photonflow::length yy, boost::units::photonflow::length zz)
        : x(xx), y(yy), z(zz) {
        photonFlowAssert(!hasNaNs());
    }
#ifndef NDEBUG
    Point3D(const Point3D &p) {
        photonFlowAssert(!p.hasNaNs());
        x = p.x; y = p.y; z = p.z;
    }

    Point3D &operator=(const Point3D &p) {
        photonFlowAssert(!p.hasNaNs());
        x = p.x; y = p.y; z = p.z;
        return *this;
    }
#endif // !NDEBUG
    Point3D operator+(const GeneralVector3D<boost::units::photonflow::length> &v) const {
        photonFlowAssert(!v.hasNaNs());
        return Point3D(x + v.x, y + v.y, z + v.z);
    }

    Point3D &operator+=(const GeneralVector3D<boost::units::photonflow::length> &v) {
        photonFlowAssert(!v.hasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    GeneralVector3D<boost::units::photonflow::length> operator-(const Point3D &p) const {
        photonFlowAssert(!p.hasNaNs());
        return GeneralVector3D<boost::units::photonflow::length>(x - p.x, y - p.y, z - p.z);
    }

    Point3D operator-(const GeneralVector3D<boost::units::photonflow::length> &v) const {
        photonFlowAssert(!v.hasNaNs());
        return Point3D(x - v.x, y - v.y, z - v.z);
    }

    Point3D &operator-=(const GeneralVector3D<boost::units::photonflow::length> &v) {
        photonFlowAssert(!v.hasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Point3D &operator+=(const Point3D &p) {
        photonFlowAssert(!p.hasNaNs());
        x += p.x; y += p.y; z += p.z;
        return *this;
    }
    Point3D operator+(const Point3D &p) const {
        photonFlowAssert(!p.hasNaNs());
        return Point3D(x + p.x, y + p.y, z + p.z);
    }
    Point3D operator* (double f) const {
        return Point3D(f*x, f*y, f*z);
    }
    Point3D &operator*=(double f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Point3D operator/ (double f) const {
        double inv = 1.f/f;
        return Point3D(inv*x, inv*y, inv*z);
    }
    Point3D &operator/=(double f) {
        double inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    boost::units::photonflow::length operator[](int i) const {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    boost::units::photonflow::length &operator[](int i) {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    bool hasNaNs() const {
        return isnan(x) || isnan(y) || isnan(z);
    }

    bool operator==(const Point3D &p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(const Point3D &p) const {
        return x != p.x || y != p.y || z != p.z;
    }
    friend std::ostream& operator<< (std::ostream &out, const Point3D &point);

    // Point Public Data
    boost::units::photonflow::length x, y, z;
};

class Normal {
public:
    // Normal Public Methods
    Normal() { x = y = z = 0.0_um; }
    Normal(boost::units::photonflow::length xx, boost::units::photonflow::length yy, boost::units::photonflow::length zz)
        : x(xx), y(yy), z(zz) {
        photonFlowAssert(!hasNaNs());
    }
    Normal operator-() const {
        return Normal(-x, -y, -z);
    }
    Normal operator+ (const Normal &n) const {
        photonFlowAssert(!n.hasNaNs());
        return Normal(x + n.x, y + n.y, z + n.z);
    }

    Normal& operator+=(const Normal &n) {
        photonFlowAssert(!n.hasNaNs());
        x += n.x; y += n.y; z += n.z;
        return *this;
    }
    Normal operator- (const Normal &n) const {
        photonFlowAssert(!n.hasNaNs());
        return Normal(x - n.x, y - n.y, z - n.z);
    }

    Normal& operator-=(const Normal &n) {
        photonFlowAssert(!n.hasNaNs());
        x -= n.x; y -= n.y; z -= n.z;
        return *this;
    }
    bool hasNaNs() const {
        return isnan(x) || isnan(y) || isnan(z);
    }
    Normal operator*(double f) const {
        return Normal(f*x, f*y, f*z);
    }

    Normal &operator*=(double f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Normal operator/(double f) const {
        photonFlowAssert(f != 0);
        double inv = 1.f/f;
        return Normal(x * inv, y * inv, z * inv);
    }

    Normal &operator/=(double f) {
        photonFlowAssert(f != 0);
        double inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    auto lengthSquared() const {
        return x*x + y*y + z*z;
    }
    auto length() const {
        return sqrt(lengthSquared());
    }

#ifndef NDEBUG
    Normal(const Normal &n) {
        photonFlowAssert(!n.hasNaNs());
        x = n.x; y = n.y; z = n.z;
    }

    Normal &operator=(const Normal &n) {
        photonFlowAssert(!n.hasNaNs());
        x = n.x; y = n.y; z = n.z;
        return *this;
    }
#endif // !NDEBUG
    explicit Normal(const GeneralVector3D<boost::units::photonflow::length> &v)
        : x(v.x), y(v.y), z(v.z) {
        photonFlowAssert(!v.hasNaNs());
    }
    boost::units::photonflow::length operator[](int i) const {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    boost::units::photonflow::length &operator[](int i) {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    bool operator==(const Normal &n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(const Normal &n) const {
        return x != n.x || y != n.y || z != n.z;
    }

    // Normal Public Data
    boost::units::photonflow::length x, y, z;
};

class Ray {
public:
    // Ray Public Methods
    Ray()
        : m_mint(0.0)
        , m_maxt(INFINITY)
        , m_time(0.0_us)
        , m_depth(0)
    { }

    Ray(const Point3D &origin, const GeneralVector3D<boost::units::photonflow::length> &direction,
        double start = 0.0,
        double end = INFINITY,
        boost::units::photonflow::time t = 0.0_us,
        int d = 0)
        : m_origin(origin)
        , m_direction(direction)
        , m_mint(start)
        , m_maxt(end)
        , m_time(t)
        , m_depth(d)
    { }

    Ray(const Point3D &origin, const GeneralVector3D<boost::units::photonflow::length> &direction, const Ray &parent,
        double start = 0.0,
        double end = INFINITY)
        : m_origin(origin)
        , m_direction(direction)
        , m_mint(start)
        , m_maxt(end)
        , m_time(parent.m_time)
        , m_depth(parent.m_depth+1)
    { }

    Point3D operator()(double t) const {
        return m_origin + m_direction * t;
    }
    bool hasNaNs() const {
        return (m_origin.hasNaNs() || m_direction.hasNaNs() ||
                isnan(m_mint) || isnan(m_maxt));
    }

    Point3D origin() const {
        return m_origin;
    }
    GeneralVector3D<boost::units::photonflow::length> direction() const {
        return m_direction;
    }

    // Ray Public Data
    Point3D m_origin;
    GeneralVector3D<boost::units::photonflow::length> m_direction;
    mutable double m_mint;
    mutable double m_maxt;
    boost::units::photonflow::time m_time;
    int m_depth;
};


class RayDifferential : public Ray {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point3D &org, const Vector3D &dir,
                    double start = 0.0,
                    double end = INFINITY,
                    boost::units::photonflow::time t = 0.0_us,
                    int d = 0)
        : Ray(org, dir, start, end, t, d) {
        hasDifferentials = false;
    }
    //    RayDifferential(const Point3D &org, const GeneralVector3D<boost::units::photonflow::length> &dir, const Ray &parent,
    //        double start, double end = INFINITY)
    //            : Ray(org, dir, start, end, parent.m_time, parent.m_depth+1) {
    //        hasDifferentials = false;
    //    }
    explicit RayDifferential(const Ray &ray) : Ray(ray) {
        hasDifferentials = false;
    }
    bool hasNaNs() const {
        return Ray::hasNaNs() ||
                (hasDifferentials && (rxOrigin.hasNaNs() || ryOrigin.hasNaNs() ||
                                      rxDirection.hasNaNs() || ryDirection.hasNaNs()));
    }
    void scaleDifferentials(double s) {
        rxOrigin = m_origin + (rxOrigin - m_origin) * s;
        ryOrigin = m_origin + (ryOrigin - m_origin) * s;
        rxDirection = m_direction + (rxDirection - m_direction) * s;
        ryDirection = m_direction + (ryDirection - m_direction) * s;
    }

    // RayDifferential Public Data
    bool hasDifferentials;
    Point3D rxOrigin, ryOrigin;
    GeneralVector3D<boost::units::photonflow::length> rxDirection, ryDirection;
};


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
        const boost::units::photonflow::length eps = 0.01_um;
        return (pt.x >= pMin.x - eps && pt.x <= pMax.x + eps &&
                pt.y >= pMin.y - eps && pt.y <= pMax.y + eps &&
                pt.z >= pMin.z - eps && pt.z <= pMax.z + eps);
    }
    void expand(boost::units::photonflow::length delta) {
        pMin -= GeneralVector3D<boost::units::photonflow::length>(delta, delta, delta);
        pMax += GeneralVector3D<boost::units::photonflow::length>(delta, delta, delta);
    }
    auto surfaceArea() const {
        GeneralVector3D<boost::units::photonflow::length> d = pMax - pMin;
        return 2.0 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    auto volume() const {
        GeneralVector3D<boost::units::photonflow::length> d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int maximumExtent() const {
        GeneralVector3D<boost::units::photonflow::length> diag = pMax - pMin;
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
    Vector3D offset(const Point3D &p) const {
        return Vector3D((p.x - pMin.x) / (pMax.x - pMin.x) * 1.0_um,
                        (p.y - pMin.y) / (pMax.y - pMin.y) * 1.0_um,
                        (p.z - pMin.z) / (pMax.z - pMin.z) * 1.0_um);
    }
    void boundingSphere(Point3D *c, boost::units::photonflow::length *rad) const;
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

//BBox makeUnion(const BBox &b, const Point3D &p);
//BBox makeUnion(const BBox &b, const BBox &b2);

// Geometry Inline Functions
template<typename T>
inline GeneralVector3D<T>::GeneralVector3D(const Point3D &p)
    : x(p.x), y(p.y), z(p.z) {
    photonFlowAssert(!hasNaNs());
}

template<typename T>
inline GeneralVector3D<T> operator*(double f, const GeneralVector3D<T> &v) {
    return v*f;
}

template<typename T>
inline auto dot(const GeneralVector3D<T> &v1, const GeneralVector3D<T> &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


inline auto absDot(const auto &v1, const auto &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    return abs(dot(v1, v2));
}

template<typename T>
auto dummy(T a, T b) {
    return a*b;
}

template<typename T>
inline auto cross(const GeneralVector3D<T> &v1, const GeneralVector3D<T> &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    auto v1x = v1.x, v1y = v1.y, v1z = v1.z;
    auto v2x = v2.x, v2y = v2.y, v2z = v2.z;
    auto v1yv2z = v1y * v2z;
    auto v1zv2y = v1z * v2y;
    auto v1zv2x = v1z * v2x;
    auto v1xv2z = v1x * v2z;
    auto v1xv2y = v1x * v2y;
    auto v1yv2x = v1y * v2x;
    return GeneralVector3D<decltype(v1x*v2x)>(
                (v1yv2z - v1zv2y),
                (v1zv2x - v1xv2z),
                (v1xv2y - v1yv2x));
}

inline auto normalize(const auto &v) {
    return v / v.length().value();
}

//inline void coordinateSystem(const GeneralVector3D &v1, GeneralVector3D *v2, GeneralVector3D *v3) {
//    if (fabsf(v1.x) > fabsf(v1.y)) {
//        double invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
//        *v2 = GeneralVector3D(-v1.z * invLen, 0.0, v1.x * invLen);
//    }
//    else {
//        double invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
//        *v2 = GeneralVector3D(0.0, v1.z * invLen, -v1.y * invLen);
//    }
//    *v3 = cross(v1, *v2);
//}

inline auto distance(const Point3D &p1, const Point3D &p2) {
    return (p1 - p2).length();
}


inline auto distanceSquared(const Point3D &p1, const Point3D &p2) {
    return (p1 - p2).lengthSquared();
}


inline Point3D operator*(double f, const Point3D &p) {
    photonFlowAssert(!p.hasNaNs());
    return p*f;
}


inline Normal operator*(double f, const Normal &n) {
    return Normal(f*n.x, f*n.y, f*n.z);
}


inline Normal normalize(const Normal &n) {
    return n / n.length().value();
}

template<typename T>
inline GeneralVector3D<T>::GeneralVector3D(const Normal &n)
    : x(n.x), y(n.y), z(n.z) {
    photonFlowAssert(!n.hasNaNs());
}


inline auto faceforward(const auto &n, const auto &v) {
    return (dot(n, v) < 0.0) ? -n : n;
}


inline const Point3D &BBox::operator[](int i) const {
    photonFlowAssert(i == 0 || i == 1);
    return (&pMin)[i];
}



inline Point3D &BBox::operator[](int i) {
    photonFlowAssert(i == 0 || i == 1);
    return (&pMin)[i];
}


//inline GeneralVector3D sphericalDirection(double sintheta,
//                                 double costheta, double phi) {
//    return GeneralVector3D(sintheta * cosf(phi),
//                  sintheta * sinf(phi),
//                  costheta);
//}


//inline GeneralVector3D sphericalDirection(double sintheta, double costheta,
//                                 double phi, const GeneralVector3D &x,
//                                 const GeneralVector3D &y, const GeneralVector3D &z) {
//    return sintheta * cosf(phi) * x +
//           sintheta * sinf(phi) * y + costheta * z;
//}


//inline double sphericalTheta(const GeneralVector3D &v) {
//    return acosf(clamp(v.z, -1.f, 1.f));
//}


//inline double sphericalPhi(const GeneralVector3D &v) {
//    double p = atan2f(v.y, v.x);
//    return (p < 0.0) ? p + 2.f*M_PI : p;
//}


#endif // PBRT_CORE_GEOMETRY_H
