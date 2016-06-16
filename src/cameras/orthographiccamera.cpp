
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

#include "orthographiccamera.h"
#include "../core/sampler.h"
#include "../core/montecarlo.h"

namespace photonflow {

// OrthographicCamera Definitions
OrthographicCamera::OrthographicCamera(const Transform &cam2world,
        const Rectangle &screenWindow, Time sopen, Time sclose,
        Length lensr, Length focald, std::shared_ptr<Film> f)
    : ProjectiveCamera(cam2world, orthographic(0., 1.), screenWindow,
                       sopen, sclose, lensr, focald, f) {
    // Compute differential changes in origin for ortho camera rays
//    dxCamera = RasterToCamera(Point3D(1.0_um, 0.0_um, 0.0_um));
//    dyCamera = RasterToCamera(Point3D(0.0_um, 1.0_um, 0.0_um));
}


double OrthographicCamera::generateRay(const CameraSample &sample, Ray *ray) const {
    auto cameraLength = 1.0_um; // TODO make camera size an option
    // TODO move camera to center, now it appears to be offset
    // Generate raster and camera samples
    Point3D Pras(sample.imageX*cameraLength, sample.imageY*cameraLength, 0);
    Point3D Pcamera;
    RasterToCamera(Pras, &Pcamera);
    *ray = Ray(Pcamera, Length3D(0.0_um, 0.0_um, 1.0_um));
    // Modify ray for depth of field
    if (lensRadius > 0.0_um) {
        // Sample point on lens
        double lensU, lensV;
        concentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);

        Length x = lensU * lensRadius;
        Length y = lensV * lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->m_direction.z;
        Point3D Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->m_origin = Point3D(x, y, 0.0_um);
        ray->m_direction = normalize(Pfocus - ray->m_origin);
    }
    ray->m_time = sample.time;
    CameraToWorld(*ray, ray);
    return 1.f;
}


//float OrthographicCamera::GenerateRayDifferential(const CameraSample &sample,
//        RayDifferential *ray) const {
//    // Compute main orthographic viewing ray

//    // Generate raster and camera samples
//    Point3D Pras(sample.imageX, sample.imageY, 0);
//    Point3D Pcamera;
//    RasterToCamera(Pras, &Pcamera);
//    *ray = RayDifferential(Pcamera, Vector(0,0,1), 0., INFINITY);

//    // Modify ray for depth of field
//    if (lensRadius > 0.) {
//        // Sample point on lens
//        float lensU, lensV;
//        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
//        lensU *= lensRadius;
//        lensV *= lensRadius;

//        // Compute point on plane of focus
//        float ft = focalDistance / ray->d.z;
//        Point3D Pfocus = (*ray)(ft);

//        // Update ray for effect of lens
//        ray->o = Point3D(lensU, lensV, 0.f);
//        ray->d = Normalize(Pfocus - ray->o);
//    }
//    ray->time = sample.time;
//    // Compute ray differentials for _OrthoCamera_
//    if (lensRadius > 0) {
//        // Compute _OrthoCamera_ ray differentials with defocus blur

//        // Sample point on lens
//        float lensU, lensV;
//        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
//        lensU *= lensRadius;
//        lensV *= lensRadius;

//        float ft = focalDistance / ray->d.z;

//        Point3D pFocus = Pcamera + dxCamera + (ft * Vector(0, 0, 1));
//        ray->rxOrigin = Point3D(lensU, lensV, 0.f);
//        ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

//        pFocus = Pcamera + dyCamera + (ft * Vector(0, 0, 1));
//        ray->ryOrigin = Point3D(lensU, lensV, 0.f);
//        ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
//    }
//    else {
//        ray->rxOrigin = ray->o + dxCamera;
//        ray->ryOrigin = ray->o + dyCamera;
//        ray->rxDirection = ray->ryDirection = ray->d;
//    }
//    ray->hasDifferentials = true;
//    CameraToWorld(*ray, ray);
//    return 1.f;
//}


//OrthographicCamera *CreateOrthographicCamera(const ParamSet &params,
//        const AnimatedTransform &cam2world, Film *film) {
//    // Extract common camera parameters from _ParamSet_
//    float shutteropen = params.FindOneFloat("shutteropen", 0.f);
//    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
//    if (shutterclose < shutteropen) {
//        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
//                shutterclose, shutteropen);
//        swap(shutterclose, shutteropen);
//    }
//    float lensradius = params.FindOneFloat("lensradius", 0.f);
//    float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
//    float frame = params.FindOneFloat("frameaspectratio",
//        float(film->xResolution)/float(film->yResolution));
//    float screen[4];
//    if (frame > 1.f) {
//        screen[0] = -frame;
//        screen[1] =  frame;
//        screen[2] = -1.f;
//        screen[3] =  1.f;
//    }
//    else {
//        screen[0] = -1.f;
//        screen[1] =  1.f;
//        screen[2] = -1.f / frame;
//        screen[3] =  1.f / frame;
//    }
//    int swi;
//    const float *sw = params.FindFloat("screenwindow", &swi);
//    if (sw && swi == 4)
//        memcpy(screen, sw, 4*sizeof(float));
//    return new OrthographicCamera(cam2world, screen, shutteropen, shutterclose,
//        lensradius, focaldistance, film);
//}

} // namespace
