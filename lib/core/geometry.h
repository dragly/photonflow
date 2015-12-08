
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
#include <ostream>
#include "../core/common.h"

class Point3D;
class Normal;

// Geometry Declarations
class Vector3D {
public:
    // Vector Public Methods
    Vector3D() { x = y = z = 0.0; }
    Vector3D(double xx, double yy, double zz)
        : x(xx), y(yy), z(zz) {
        photonFlowAssert(!hasNaNs());
    }
    bool hasNaNs() const { return isnan(x) || isnan(y) || isnan(z); }
    explicit Vector3D(const Point3D &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    Vector3D(const Vector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x = v.x; y = v.y; z = v.z;
    }

    Vector3D &operator=(const Vector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
#endif // !NDEBUG
    Vector3D operator+(const Vector3D &v) const {
        photonFlowAssert(!v.hasNaNs());
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    Vector3D& operator+=(const Vector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector3D operator-(const Vector3D &v) const {
        photonFlowAssert(!v.hasNaNs());
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    Vector3D& operator-=(const Vector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Vector3D operator*(double f) const { return Vector3D(f*x, f*y, f*z); }

    Vector3D &operator*=(double f) {
        photonFlowAssert(!isnan(f));
        x *= f; y *= f; z *= f;
        return *this;
    }
    Vector3D operator/(double f) const {
        photonFlowAssert(f != 0);
        double inv = 1.f / f;
        return Vector3D(x * inv, y * inv, z * inv);
    }

    Vector3D &operator/=(double f) {
        photonFlowAssert(f != 0);
        double inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    Vector3D operator-() const { return Vector3D(-x, -y, -z); }
    double operator[](int i) const {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    double &operator[](int i)
    {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    double lengthSquared() const
    {
        return x*x + y*y + z*z;
    }

    double length() const
    {
        return sqrtf(lengthSquared());
    }

    explicit Vector3D(const Normal &n);

    bool operator==(const Vector3D &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    bool operator!=(const Vector3D &v) const {
        return x != v.x || y != v.y || z != v.z;
    }

    friend std::ostream& operator<< (std::ostream &out, const Vector3D &vector);

    Vector3D perpendicular() const
    {
        if(x == 0.0 && y == 0.0) {
            if(z == 0.0) {
                return Vector3D(0.0, 0.0, 0.0);
            }
            return Vector3D(0.0, 1.0, 0.0);
        }
        return Vector3D(-y, x, 0.0);
    }

    Vector3D normalized() {
        return *this / length();
    }

    // Vector Public Data
    double x, y, z;
};

class Point3D {
public:
    // Point Public Methods
    Point3D() { x = y = z = 0.0; }
    Point3D(double xx, double yy, double zz)
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
    Point3D operator+(const Vector3D &v) const {
        photonFlowAssert(!v.hasNaNs());
        return Point3D(x + v.x, y + v.y, z + v.z);
    }

    Point3D &operator+=(const Vector3D &v) {
        photonFlowAssert(!v.hasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector3D operator-(const Point3D &p) const {
        photonFlowAssert(!p.hasNaNs());
        return Vector3D(x - p.x, y - p.y, z - p.z);
    }

    Point3D operator-(const Vector3D &v) const {
        photonFlowAssert(!v.hasNaNs());
        return Point3D(x - v.x, y - v.y, z - v.z);
    }

    Point3D &operator-=(const Vector3D &v) {
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
    double operator[](int i) const {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    double &operator[](int i) {
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
    double x, y, z;
};


class Normal {
public:
    // Normal Public Methods
    Normal() { x = y = z = 0.0; }
    Normal(double xx, double yy, double zz)
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
    double LengthSquared() const { return x*x + y*y + z*z; }
    double Length() const        { return sqrtf(LengthSquared()); }

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
    explicit Normal(const Vector3D &v)
      : x(v.x), y(v.y), z(v.z) {
        photonFlowAssert(!v.hasNaNs());
    }
    double operator[](int i) const {
        photonFlowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    double &operator[](int i) {
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
    double x, y, z;
};


class Ray {
public:
    // Ray Public Methods
    Ray() : m_mint(0.0), m_maxt(INFINITY), m_time(0.0), m_depth(0) { }
    Ray(const Point3D &origin, const Vector3D &direction,
        double start = 0.0, double end = INFINITY, double t = 0.0, int d = 0)
        : m_origin(origin), m_direction(direction), m_mint(start), m_maxt(end), m_time(t), m_depth(d) { }
    Ray(const Point3D &origin, const Vector3D &direction, const Ray &parent,
        double start = 0.0, double end = INFINITY)
        : m_origin(origin), m_direction(direction), m_mint(start), m_maxt(end),
          m_time(parent.m_time), m_depth(parent.m_depth+1) { }
    Point3D operator()(double t) const { return m_origin + m_direction * t; }
    bool hasNaNs() const {
        return (m_origin.hasNaNs() || m_direction.hasNaNs() ||
                isnan(m_mint) || isnan(m_maxt));
    }

    Point3D origin() const {
        return m_origin;
    }
    Vector3D direction() const {
        return m_direction;
    }

    // Ray Public Data
    Point3D m_origin;
    Vector3D m_direction;
    mutable double m_mint, m_maxt;
    double m_time;
    int m_depth;
};


class RayDifferential : public Ray {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point3D &org, const Vector3D &dir, double start,
        double end = INFINITY, double t = 0.0, int d = 0)
            : Ray(org, dir, start, end, t, d) {
        hasDifferentials = false;
    }
    RayDifferential(const Point3D &org, const Vector3D &dir, const Ray &parent,
        double start, double end = INFINITY)
            : Ray(org, dir, start, end, parent.m_time, parent.m_depth+1) {
        hasDifferentials = false;
    }
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
    Vector3D rxDirection, ryDirection;
};


class BBox {
public:
    // BBox Public Methods
    BBox() {
        pMin = Point3D( INFINITY,  INFINITY,  INFINITY);
        pMax = Point3D(-INFINITY, -INFINITY, -INFINITY);
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
        const double eps = 0.01;
        return (pt.x >= pMin.x - eps && pt.x <= pMax.x + eps &&
                pt.y >= pMin.y - eps && pt.y <= pMax.y + eps &&
                pt.z >= pMin.z - eps && pt.z <= pMax.z + eps);
    }
    void expand(double delta) {
        pMin -= Vector3D(delta, delta, delta);
        pMax += Vector3D(delta, delta, delta);
    }
    double surfaceArea() const {
        Vector3D d = pMax - pMin;
        return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    double volume() const {
        Vector3D d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int maximumExtent() const {
        Vector3D diag = pMax - pMin;
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
        return Vector3D((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }
    void boundingSphere(Point3D *c, double *rad) const;
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



// Geometry Inline Functions
inline Vector3D::Vector3D(const Point3D &p)
    : x(p.x), y(p.y), z(p.z) {
    photonFlowAssert(!hasNaNs());
}


inline Vector3D operator*(double f, const Vector3D &v) { return v*f; }
inline double dot(const Vector3D &v1, const Vector3D &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


inline double absDot(const Vector3D &v1, const Vector3D &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    return fabsf(dot(v1, v2));
}


inline Vector3D cross(const Vector3D &v1, const Vector3D &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3D((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector3D cross(const Vector3D &v1, const Normal &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3D((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector3D cross(const Normal &v1, const Vector3D &v2) {
    photonFlowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3D((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector3D normalize(const Vector3D &v) { return v / v.length(); }
inline void CoordinateSystem(const Vector3D &v1, Vector3D *v2, Vector3D *v3) {
    if (fabsf(v1.x) > fabsf(v1.y)) {
        double invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
        *v2 = Vector3D(-v1.z * invLen, 0.0, v1.x * invLen);
    }
    else {
        double invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
        *v2 = Vector3D(0.0, v1.z * invLen, -v1.y * invLen);
    }
    *v3 = cross(v1, *v2);
}


inline double distance(const Point3D &p1, const Point3D &p2) {
    return (p1 - p2).length();
}


inline double distanceSquared(const Point3D &p1, const Point3D &p2) {
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
    return n / n.Length();
}


inline Vector3D::Vector3D(const Normal &n)
  : x(n.x), y(n.y), z(n.z) {
    photonFlowAssert(!n.hasNaNs());
}


inline double dot(const Normal &n1, const Vector3D &v2) {
    photonFlowAssert(!n1.hasNaNs() && !v2.hasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}


inline double dot(const Vector3D &v1, const Normal &n2) {
    photonFlowAssert(!v1.hasNaNs() && !n2.hasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}


inline double dot(const Normal &n1, const Normal &n2) {
    photonFlowAssert(!n1.hasNaNs() && !n2.hasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}


inline double absDot(const Normal &n1, const Vector3D &v2) {
    photonFlowAssert(!n1.hasNaNs() && !v2.hasNaNs());
    return fabsf(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}


inline double absDot(const Vector3D &v1, const Normal &n2) {
    photonFlowAssert(!v1.hasNaNs() && !n2.hasNaNs());
    return fabsf(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}


inline double absDot(const Normal &n1, const Normal &n2) {
    photonFlowAssert(!n1.hasNaNs() && !n2.hasNaNs());
    return fabsf(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}


inline Normal faceforward(const Normal &n, const Vector3D &v) {
    return (dot(n, v) < 0.0) ? -n : n;
}


inline Normal faceforward(const Normal &n, const Normal &n2) {
    return (dot(n, n2) < 0.0) ? -n : n;
}



inline Vector3D faceforward(const Vector3D &v, const Vector3D &v2) {
    return (dot(v, v2) < 0.0) ? -v : v;
}



inline Vector3D faceforward(const Vector3D &v, const Normal &n2) {
    return (dot(v, n2) < 0.0) ? -v : v;
}


inline const Point3D &BBox::operator[](int i) const {
    photonFlowAssert(i == 0 || i == 1);
    return (&pMin)[i];
}



inline Point3D &BBox::operator[](int i) {
    photonFlowAssert(i == 0 || i == 1);
    return (&pMin)[i];
}


inline Vector3D sphericalDirection(double sintheta,
                                 double costheta, double phi) {
    return Vector3D(sintheta * cosf(phi),
                  sintheta * sinf(phi),
                  costheta);
}


inline Vector3D sphericalDirection(double sintheta, double costheta,
                                 double phi, const Vector3D &x,
                                 const Vector3D &y, const Vector3D &z) {
    return sintheta * cosf(phi) * x +
           sintheta * sinf(phi) * y + costheta * z;
}


inline double sphericalTheta(const Vector3D &v) {
    return acosf(clamp(v.z, -1.f, 1.f));
}


inline double sphericalPhi(const Vector3D &v) {
    double p = atan2f(v.y, v.x);
    return (p < 0.0) ? p + 2.f*M_PI : p;
}


#endif // PBRT_CORE_GEOMETRY_H
