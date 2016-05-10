#ifndef NORMAL_H
#define NORMAL_H

#include "../core/units.h"
#include "../core/common.h"

using namespace photonflow::literals;

namespace photonflow {

class Normal
{
public:
    // Normal Public Methods
    Normal();
    Normal(photonflow::Length xx, photonflow::Length yy, photonflow::Length zz);
    Normal operator-() const;
    Normal operator+ (const Normal &n) const;
    Normal& operator+=(const Normal &n);
    Normal operator- (const Normal &n) const;

    Normal& operator-=(const Normal &n);
    bool hasNaNs() const;
    Normal operator*(double f) const;

    Normal &operator*=(double f);
    Normal operator/(double f) const;

    Normal &operator/=(double f);
    auto lengthSquared() const;
    auto length() const;

#ifndef NDEBUG
    Normal(const Normal &n);

    Normal &operator=(const Normal &n);
#endif // !NDEBUG
    explicit Normal(const Length3D &v);
    photonflow::Length operator[](int i) const;

    photonflow::Length &operator[](int i);

    bool operator==(const Normal &n) const;
    bool operator!=(const Normal &n) const;

    // Normal Public Data
    photonflow::Length x, y, z;
};

} // namespace

#include "vector3d.h"

namespace photonflow {

inline Normal::Normal() { x = y = z = 0.0_um; }

inline Normal::Normal(photonflow::Length xx, photonflow::Length yy, photonflow::Length zz)
    : x(xx), y(yy), z(zz) {
    photonflowAssert(!hasNaNs());
}

inline Normal Normal::operator-() const {
    return Normal(-x, -y, -z);
}

inline Normal Normal::operator+(const Normal &n) const {
    photonflowAssert(!n.hasNaNs());
    return Normal(x + n.x, y + n.y, z + n.z);
}

inline Normal &Normal::operator+=(const Normal &n) {
    photonflowAssert(!n.hasNaNs());
    x += n.x; y += n.y; z += n.z;
    return *this;
}

inline Normal Normal::operator-(const Normal &n) const {
    photonflowAssert(!n.hasNaNs());
    return Normal(x - n.x, y - n.y, z - n.z);
}

inline Normal &Normal::operator-=(const Normal &n) {
    photonflowAssert(!n.hasNaNs());
    x -= n.x; y -= n.y; z -= n.z;
    return *this;
}

inline bool Normal::hasNaNs() const {
    return isnan(x) || isnan(y) || isnan(z);
}

inline Normal Normal::operator*(double f) const {
    return Normal(f*x, f*y, f*z);
}

inline Normal &Normal::operator*=(double f) {
    x *= f; y *= f; z *= f;
    return *this;
}

inline Normal Normal::operator/(double f) const {
    photonflowAssert(f != 0);
    double inv = 1.f/f;
    return Normal(x * inv, y * inv, z * inv);
}

inline Normal &Normal::operator/=(double f) {
    photonflowAssert(f != 0);
    double inv = 1.f/f;
    x *= inv; y *= inv; z *= inv;
    return *this;
}

inline auto Normal::lengthSquared() const {
    return x*x + y*y + z*z;
}

inline auto Normal::length() const {
    return sqrt(lengthSquared());
}

#ifndef NDEBUG
inline Normal::Normal(const Normal &n) {
    photonflowAssert(!n.hasNaNs());
    x = n.x; y = n.y; z = n.z;
}
inline Normal &Normal::operator=(const Normal &n) {
    photonflowAssert(!n.hasNaNs());
    x = n.x; y = n.y; z = n.z;
    return *this;
}
#endif

inline Normal::Normal(const Length3D &v)
    : x(v.x), y(v.y), z(v.z) {
    photonflowAssert(!v.hasNaNs());
}

inline photonflow::Length Normal::operator[](int i) const {
    photonflowAssert(i >= 0 && i <= 2);
    return (&x)[i];
}

inline photonflow::Length &Normal::operator[](int i) {
    photonflowAssert(i >= 0 && i <= 2);
    return (&x)[i];
}

inline bool Normal::operator==(const Normal &n) const {
    return x == n.x && y == n.y && z == n.z;
}

inline bool Normal::operator!=(const Normal &n) const {
    return x != n.x || y != n.y || z != n.z;
}

} // namespace

#endif // NORMAL_H
