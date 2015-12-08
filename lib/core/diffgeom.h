
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

#ifndef PBRT_CORE_DIFFGEOM_H
#define PBRT_CORE_DIFFGEOM_H

// core/diffgeom.h*
#include "../core/common.h"
#include "geometry.h"

class Shape;

// DifferentialGeometry Declarations
struct DifferentialGeometry {
    DifferentialGeometry() { 
        u = v = dudx = dvdx = dudy = dvdy = 0.; 
        shape = NULL; 
    }
    // DifferentialGeometry Public Methods
    DifferentialGeometry(const Point3D &P, const Vector3D &DPDU,
            const Vector3D &DPDV, const Normal &DNDU,
            const Normal &DNDV, double uu, double vv,
            const Shape *sh);
    void ComputeDifferentials(const RayDifferential &r) const;

    // DifferentialGeometry Public Data
    Point3D p;
    Normal nn;
    double u, v;
    const Shape *shape;
    Vector3D dpdu, dpdv;
    Normal dndu, dndv;
    mutable Vector3D dpdx, dpdy;
    mutable double dudx, dvdx, dudy, dvdy;
};



#endif // PBRT_CORE_DIFFGEOM_H
