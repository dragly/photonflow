
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
        const Rectangle &screenWindow, boost::units::photonflow::time sopen, boost::units::photonflow::time sclose,
        boost::units::photonflow::length lensr, boost::units::photonflow::length focald, double fov, std::shared_ptr<Film> f)
    : ProjectiveCamera(cam2world, perspective(fov, 1e-2f, 1000.0),
                       screenWindow, sopen, sclose, lensr, focald, f) {
    // Compute differential changes in origin for perspective camera rays
    dxCamera = RasterToCamera(Point3D(1.0_um,0.0_um,0.0_um)) - RasterToCamera(Point3D(0.0_um,0.0_um,0.0_um));
    dyCamera = RasterToCamera(Point3D(0.0_um,1.0_um,0.0_um)) - RasterToCamera(Point3D(0.0_um,0.0_um,0.0_um));
}


double PerspectiveCamera::generateRay(const CameraSample &sample,
                                     Ray *ray) const {
    auto cameraLength = 1.0_um;

    // Generate raster and camera samples
    Point3D Pras(sample.imageX*cameraLength, sample.imageY*cameraLength, 0);
    Point3D Pcamera;
    RasterToCamera(Pras, &Pcamera);
    *ray = Ray(Point3D(), normalize(Vector3D(Pcamera)));
    // Modify ray for depth of field
    if (lensRadius > 0.0_um) {
        // Sample point on lens
        double lensU, lensV;
        concentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);

        boost::units::photonflow::length x = lensU * lensRadius;
        boost::units::photonflow::length y = lensU * lensRadius;

        // Compute point on plane of focus
        double ft = focalDistance / ray->m_direction.z;
        Point3D Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->m_origin = Point3D(x, y, 0.0_um);
        ray->m_direction = normalize(Pfocus - ray->m_origin);
    }
    ray->m_time = sample.time;
    CameraToWorld(*ray, ray);
    return 1.f;
}


//double PerspectiveCamera::generateRayDifferential(const CameraSample &sample,
//                                                 RayDifferential *ray) const {
//    // Generate raster and camera samples
//    Point3D Pras(sample.imageX*1.0_um, sample.imageY*1.0_um, 0.0_um); // TODO: Consider coordinate system
//    Point3D Pcamera;
//    RasterToCamera(Pras, &Pcamera);
//    Vector3D dir = normalize(Vector3D(Pcamera.x, Pcamera.y, Pcamera.z));
//    *ray = RayDifferential(Point3D(0,0,0), dir, 0.0, INFINITY);
//    // Modify ray for depth of field
//    if (lensRadius > 0.) {
//        // Sample point on lens
//        double lensU, lensV;
//        concentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
//        lensU *= lensRadius;
//        lensV *= lensRadius;

//        // Compute point on plane of focus
//        double ft = focalDistance / ray->m_direction.z;
//        Point3D Pfocus = (*ray)(ft);

//        // Update ray for effect of lens
//        ray->m_origin = Point3D(lensU, lensV, 0.0);
//        ray->m_direction = normalize(Pfocus - ray->m_origin);
//    }

//    // Compute offset rays for _PerspectiveCamera_ ray differentials
//    if (lensRadius > 0.) {
//        // Compute _PerspectiveCamera_ ray differentials with defocus blur

//        // Sample point on lens
//        double lensU, lensV;
//        concentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
//        lensU *= lensRadius;
//        lensV *= lensRadius;

//        Vector3D dx = normalize(Vector3D(Pcamera + dxCamera));
//        double ft = focalDistance / dx.z;
//        Point3D pFocus = Point3D(0,0,0) + (ft * dx);
//        ray->rxOrigin = Point3D(lensU, lensV, 0.0);
//        ray->rxDirection = normalize(pFocus - ray->rxOrigin);

//        Vector3D dy = normalize(Vector3D(Pcamera + dyCamera));
//        ft = focalDistance / dy.z;
//        pFocus = Point3D(0,0,0) + (ft * dy);
//        ray->ryOrigin = Point3D(lensU, lensV, 0.0);
//        ray->ryDirection = normalize(pFocus - ray->ryOrigin);
//    }
//    else {
//        ray->rxOrigin = ray->ryOrigin = ray->m_origin;
//        ray->rxDirection = normalize(Vector3D(Pcamera) + dxCamera);
//        ray->ryDirection = normalize(Vector3D(Pcamera) + dyCamera);
//    }

//    ray->m_time = sample.time;
//    CameraToWorld(*ray, ray);
//    ray->hasDifferentials = true;
//    return 1.f;
//}
