#ifndef RAY_H
#define RAY_H

#include "../core/common.h"
#include "../core/units.h"

#include "point3d.h"
#include "vector3d.h"

using namespace photonflow::literals;

namespace photonflow {

class Ray {
public:
    // Ray Public Methods
    Ray()
        : m_mint(0.0)
        , m_maxt(INFINITY)
        , m_time(0.0_us)
        , m_depth(0)
    { }

    Ray(const Point3D &origin, const Length3D &direction,
        double start = 0.0,
        double end = INFINITY,
        photonflow::Time t = 0.0_us,
        int d = 0)
        : m_origin(origin)
        , m_direction(direction)
        , m_mint(start)
        , m_maxt(end)
        , m_time(t)
        , m_depth(d)
    { }

    Ray(const Point3D &origin, const Length3D &direction, const Ray &parent,
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
        return (m_origin.hasNaNs() || m_direction.hasNaNs() || std::isnan(m_mint) || std::isnan(m_maxt));
    }

    Point3D origin() const {
        return m_origin;
    }
    Length3D direction() const {
        return m_direction;
    }

    // Ray Public Data
    Point3D m_origin;
    Length3D m_direction;
    mutable double m_mint;
    mutable double m_maxt;
    photonflow::Time m_time;
    int m_depth;
};



} // namespace

#endif // RAY_H
