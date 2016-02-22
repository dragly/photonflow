#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "../core/common.h"
#include <iostream>

namespace photonflow {

class Point3D;
class Normal;

template<typename T>
class BaseVector3D {
public:
    BaseVector3D()
        : x(0.0*T())
        , y(0.0*T())
        , z(0.0*T())
    {

    }

    BaseVector3D(T xx, T yy, T zz)
        : x(xx)
        , y(yy)
        , z(zz) {
        photonflowAssert(!hasNaNs());
    }

    //    bool hasNaNs() const { return isnan(x) || isnan(y) || isnan(z); }
    bool hasNaNs() const { return false; }
    explicit BaseVector3D(const Point3D &p);

#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    BaseVector3D(const BaseVector3D &v) {
        photonflowAssert(!v.hasNaNs());
        x = v.x; y = v.y; z = v.z;
    }

    template<typename U>
    BaseVector3D &operator=(const BaseVector3D<U> &v) {
        photonflowAssert(!v.hasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
#endif // !NDEBUG
    BaseVector3D operator+(const BaseVector3D &v) const {
        photonflowAssert(!v.hasNaNs());
        return BaseVector3D(x + v.x, y + v.y, z + v.z);
    }

    BaseVector3D& operator+=(const BaseVector3D &v) {
        photonflowAssert(!v.hasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    BaseVector3D operator-(const BaseVector3D &v) const {
        photonflowAssert(!v.hasNaNs());
        return BaseVector3D(x - v.x, y - v.y, z - v.z);
    }

    BaseVector3D& operator-=(const BaseVector3D &v) {
        photonflowAssert(!v.hasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    template<typename U>
    auto operator*(U f) const {
        return BaseVector3D<decltype(f*x)>(f*x, f*y, f*z);
    }

    template<typename U>
    BaseVector3D &operator*=(U f) {
        photonflowAssert(!isnan(f));
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }

    template<typename U>
    auto operator/(U f) const {
        photonflowAssert(f != 0.0 * U());
        auto inv = 1.0 / f;
        return (*this)*inv;
    }

    template<typename U>
    BaseVector3D &operator/=(U f) {
        photonflowAssert(f != 0.0 * U());
        auto inv = 1.0 / f;
        *this *= inv;
        return *this;
    }
    BaseVector3D operator-() const { return BaseVector3D(-x, -y, -z); }
    T operator[](int i) const {
        photonflowAssert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    T &operator[](int i)
    {
        photonflowAssert(i >= 0 && i <= 2);
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

    explicit BaseVector3D(const Normal &n);

    bool operator==(const BaseVector3D &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    bool operator!=(const BaseVector3D &v) const {
        return x != v.x || y != v.y || z != v.z;
    }

    template<typename U>
    friend std::ostream& operator<< (std::ostream &out, const BaseVector3D<U> &vector);

    BaseVector3D perpendicular() const
    {
        if(x == 0.0 * T() && y == 0.0 * T()) {
            if(z == 0.0 * T()) {
                return BaseVector3D(0.0 * T(), 0.0 * T(), 0.0 * T());
            }
            return BaseVector3D(0.0 * T(), 1.0 * T(), 0.0 * T());
        }
        return BaseVector3D(-y, x, 0.0 * T());
    }

    auto normalized() {
        return *this / length();
    }

    // Vector Public Data
    T x = 0.0 * T();
    T y = 0.0 * T();
    T z = 0.0 * T();
};

template<typename T>
inline BaseVector3D<T> operator*(double f, const BaseVector3D<T> &v) {
    return v*f;
}

}

#include "point3d.h"

namespace photonflow {

template<typename T>
inline BaseVector3D<T>::BaseVector3D(const Point3D &p)
    : x(p.x), y(p.y), z(p.z) {
    photonflowAssert(!hasNaNs());
}

} // namespace

#endif // VECTOR3D_H
