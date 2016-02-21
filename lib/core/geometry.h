
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
#include "../geometry/bbox.h"
#include "../geometry/normal.h"
#include "../geometry/ray.h"
#include "../geometry/raydifferential.h"
#include "../geometry/point3d.h"
#include "../geometry/vector3d.h"

// Geometry Inline Functions

template<typename T>
inline auto dot(const BaseVector3D<T> &v1, const BaseVector3D<T> &v2) {
    photonflowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline auto absDot(const auto &v1, const auto &v2) {
    photonflowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    return abs(dot(v1, v2));
}

template<typename T>
auto dummy(T a, T b) {
    return a*b;
}

template<typename T>
inline auto cross(const BaseVector3D<T> &v1, const BaseVector3D<T> &v2) {
    photonflowAssert(!v1.hasNaNs() && !v2.hasNaNs());
    auto v1x = v1.x, v1y = v1.y, v1z = v1.z;
    auto v2x = v2.x, v2y = v2.y, v2z = v2.z;
    auto v1yv2z = v1y * v2z;
    auto v1zv2y = v1z * v2y;
    auto v1zv2x = v1z * v2x;
    auto v1xv2z = v1x * v2z;
    auto v1xv2y = v1x * v2y;
    auto v1yv2x = v1y * v2x;
    return BaseVector3D<decltype(v1x*v2x)>(
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
    photonflowAssert(!p.hasNaNs());
    return p*f;
}


inline Normal operator*(double f, const Normal &n) {
    return Normal(f*n.x, f*n.y, f*n.z);
}


inline Normal normalize(const Normal &n) {
    return n / n.length().value();
}

template<typename T>
inline BaseVector3D<T>::BaseVector3D(const Normal &n)
    : x(n.x), y(n.y), z(n.z) {
    photonflowAssert(!n.hasNaNs());
}


inline auto faceforward(const auto &n, const auto &v) {
    return (dot(n, v) < 0.0) ? -n : n;
}


inline const Point3D &BBox::operator[](int i) const {
    photonflowAssert(i == 0 || i == 1);
    return (&pMin)[i];
}



inline Point3D &BBox::operator[](int i) {
    photonflowAssert(i == 0 || i == 1);
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
