
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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
#include "../core/common.h"
#include "../core/geometry.h"
#include "../core/transform.h"
#include "../core/diffgeom.h"
//#include "../core/memory.h"

// Shape Declarations
class Shape {
public:
    // Shape Interface
    Shape(const Transform *o2w, const Transform *w2o, bool ro);
    virtual ~Shape();
    virtual BoundingBox objectBound() const = 0;
    virtual BoundingBox worldBound() const;
    virtual bool canIntersect() const;
    virtual void refine(vector<Shape> &refined) const;
    virtual bool intersect(const Ray &ray, double *tHit,
                           double *rayEpsilon, DifferentialGeometry *dg) const;
    virtual bool intersectP(const Ray &ray) const;
    virtual void shadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const {
        UNUSED(obj2world);
        *dgShading = dg;
    }
    virtual double area() const;
    virtual Point3D sample(double u1, double u2, Normal *Ns) const {
        UNUSED(u1);
        UNUSED(u2);
        UNUSED(Ns);
        Severe("Unimplemented Shape::Sample() method called");
        return Point3D();
    }
//    virtual double probabilityDistributionFunction(const Point3D &Pshape) const {
//        UNUSED(Pshape);
//        return 1.f / area();
//    }
    virtual Point3D sample(const Point3D &P, double u1, double u2,
                         Normal *Ns) const {
        UNUSED(P);
        return sample(u1, u2, Ns);
    }
//    virtual double probabilityDistributionFunction(const Point3D &p, const Vector3D &wi) const;

    // Shape Public Data
    const Transform *ObjectToWorld, *WorldToObject;
    const bool ReverseOrientation, TransformSwapsHandedness;
    const uint32_t shapeId;
    static uint32_t nextshapeId;
};



#endif // PBRT_CORE_SHAPE_H
