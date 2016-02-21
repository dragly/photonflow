#ifndef RAYDIFFERENTIAL_H
#define RAYDIFFERENTIAL_H

#include "ray.h"

class RayDifferential : public Ray {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point3D &org, const Length3D &dir,
                    double start = 0.0,
                    double end = INFINITY,
                    boost::units::photonflow::Time t = 0.0_us,
                    int d = 0)
        : Ray(org, dir, start, end, t, d) {
        hasDifferentials = false;
    }
    //    RayDifferential(const Point3D &org, const Vector3D &dir, const Ray &parent,
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
    Length3D rxDirection, ryDirection;
};

#endif // RAYDIFFERENTIAL_H
