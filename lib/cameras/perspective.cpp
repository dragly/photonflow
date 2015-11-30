
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


// cameras/perspective.cpp*
#include "stdafx.h"
#include "perspective.h"
//#include "paramset.h"
#include "../core/sampler.h"
#include "../core/montecarlo.h"

// PerspectiveCamera Method Definitions
PerspectiveCamera:: PerspectiveCamera(const Transform &cam2world,
        const Rectangle &screenWindow, float sopen, float sclose,
        float lensr, float focald, float fov, std::shared_ptr<Film> f)
    : ProjectiveCamera(cam2world, Perspective(fov, 1e-2f, 1000.f),
                       screenWindow, sopen, sclose, lensr, focald, f) {
    // Compute differential changes in origin for perspective camera rays
    dxCamera = RasterToCamera(Point(1,0,0)) - RasterToCamera(Point(0,0,0));
    dyCamera = RasterToCamera(Point(0,1,0)) - RasterToCamera(Point(0,0,0));
}


float PerspectiveCamera::GenerateRay(const CameraSample &sample,
                                     Ray *ray) const {
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    *ray = Ray(Point(0,0,0), Normalize(Vector(Pcamera)), 0.f, INFINITY);
    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->m_direction.z;
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->m_origin = Point(lensU, lensV, 0.f);
        ray->m_direction = Normalize(Pfocus - ray->m_origin);
    }
    ray->m_time = sample.time;
    CameraToWorld(*ray, ray);
    return 1.f;
}


float PerspectiveCamera::GenerateRayDifferential(const CameraSample &sample,
                                                 RayDifferential *ray) const {
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    Vector dir = Normalize(Vector(Pcamera.x, Pcamera.y, Pcamera.z));
    *ray = RayDifferential(Point(0,0,0), dir, 0.f, INFINITY);
    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->m_direction.z;
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->m_origin = Point(lensU, lensV, 0.f);
        ray->m_direction = Normalize(Pfocus - ray->m_origin);
    }

    // Compute offset rays for _PerspectiveCamera_ ray differentials
    if (lensRadius > 0.) {
        // Compute _PerspectiveCamera_ ray differentials with defocus blur

        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        Vector dx = Normalize(Vector(Pcamera + dxCamera));
        float ft = focalDistance / dx.z;
        Point pFocus = Point(0,0,0) + (ft * dx);
        ray->rxOrigin = Point(lensU, lensV, 0.f);
        ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

        Vector dy = Normalize(Vector(Pcamera + dyCamera));
        ft = focalDistance / dy.z;
        pFocus = Point(0,0,0) + (ft * dy);
        ray->ryOrigin = Point(lensU, lensV, 0.f);
        ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
    }
    else {
        ray->rxOrigin = ray->ryOrigin = ray->m_origin;
        ray->rxDirection = Normalize(Vector(Pcamera) + dxCamera);
        ray->ryDirection = Normalize(Vector(Pcamera) + dyCamera);
    }

    ray->m_time = sample.time;
    CameraToWorld(*ray, ray);
    ray->hasDifferentials = true;
    return 1.f;
}
