
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


// core/geometry.cpp*
#include "stdafx.h"
#include "geometry.h"

// BBox Method Definitions
BoundingBox makeUnion(const BoundingBox &b, const Point3D &p) {
    BoundingBox ret = b;
    ret.pMin.x = min(b.pMin.x, p.x);
    ret.pMin.y = min(b.pMin.y, p.y);
    ret.pMin.z = min(b.pMin.z, p.z);
    ret.pMax.x = max(b.pMax.x, p.x);
    ret.pMax.y = max(b.pMax.y, p.y);
    ret.pMax.z = max(b.pMax.z, p.z);
    return ret;
}

BoundingBox makeUnion(const BoundingBox &b, const BoundingBox &b2) {
    BoundingBox ret;
    ret.pMin.x = min(b.pMin.x, b2.pMin.x);
    ret.pMin.y = min(b.pMin.y, b2.pMin.y);
    ret.pMin.z = min(b.pMin.z, b2.pMin.z);
    ret.pMax.x = max(b.pMax.x, b2.pMax.x);
    ret.pMax.y = max(b.pMax.y, b2.pMax.y);
    ret.pMax.z = max(b.pMax.z, b2.pMax.z);
    return ret;
}

void BoundingBox::boundingSphere(Point3D *c, photonflow::Length *rad) const
{
    *c = .5f * pMin + .5f * pMax;
    *rad = inside(*c) ? distance(*c, pMax) : 0.0_um;
}


bool BoundingBox::intersectP(const Ray &ray,
                      double *hitt0,
                      double *hitt1) const
{
    double t0 = ray.m_mint;
    double t1 = ray.m_maxt;
    for (int i = 0; i < 3; ++i) {
        // Update interval for _i_th bounding box slab
        auto invRayDir =  1.0 / ray.m_direction[i];
        double tNear = (pMin[i] - ray.m_origin[i]) * invRayDir;
        double tFar  = (pMax[i] - ray.m_origin[i]) * invRayDir;

        // Update parametric interval from slab intersection $t$s
        if (tNear > tFar) swap(tNear, tFar);
        t0 = (tNear > t0) ? tNear : t0;
        t1 = (tFar  < t1) ? tFar  : t1;
        if (t0 > t1) return false;
    }
    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;
}

std::ostream& operator<< (std::ostream &out, const Length3D &vector)
{
    (void)vector;
    // TODO Reintroduce output
//    out << vector.x << ", " << vector.y << ", " << vector.z;
    return out;
}

std::ostream& operator<< (std::ostream &out, const Point3D &point)
{
    (void)point;
    out << point.x << ", " << point.y << ", " << point.z;
    return out;
}
