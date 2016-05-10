
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


// core/camera.cpp*
#include "stdafx.h"
#include "camera.h"
#include "film.h"
#include "montecarlo.h"
#include "sampler.h"

using namespace std;

namespace photonflow {

// Camera Method Definitions
Camera::~Camera() {
//    delete film;
}


Camera::Camera(const Transform &cam2world,
               photonflow::Time sopen,
               photonflow::Time sclose,
               std::shared_ptr<Film> f)
    : CameraToWorld(cam2world), shutterOpen(sopen), shutterClose(sclose) {
    film = f;
//    if (CameraToWorld.hasScale())
//        Warning("Scaling detected in world-to-camera transformation!\n"
//                "The system has numerous assumptions, implicit and explicit,\n"
//                "that this transform will have no scale factors in it.\n"
//                "Proceed at your own risk; your image may have errors or\n"
//                "the system may crash as a result of this.");
}


//double Camera::generateRayDifferential(const CameraSample &sample,
//                                      RayDifferential *rd) const {
//    double wt = generateRay(sample, rd);
//    // Find ray after shifting one pixel in the $x$ direction
//    CameraSample sshift = sample;
//    ++(sshift.imageX);
//    Ray rx;
//    double wtx = generateRay(sshift, &rx);
//    rd->rxOrigin = rx.m_origin;
//    rd->rxDirection = rx.m_direction;

//    // Find ray after shifting one pixel in the $y$ direction
//    --(sshift.imageX);
//    ++(sshift.imageY);
//    Ray ry;
//    double wty = generateRay(sshift, &ry);
//    rd->ryOrigin = ry.m_origin;
//    rd->ryDirection = ry.m_direction;
//    if (wtx == 0.0 || wty == 0.0) return 0.0;
//    rd->hasDifferentials = true;
//    return wt;
//}


ProjectiveCamera::ProjectiveCamera(const Transform &cam2world,
        const Transform &proj, const Rectangle &screenWindow, photonflow::Time sopen,
        photonflow::Time sclose, photonflow::Length lensr, photonflow::Length focald, shared_ptr<Film> f)
    : Camera(cam2world, sopen, sclose, f) {
    // Initialize depth of field parameters
    lensRadius = lensr;
    focalDistance = focald;

    // Compute projective camera transformations
    CameraToScreen = proj;

    // Compute projective camera screen transformations
    ScreenToRaster = scale(double(film->xResolution),
                           double(film->yResolution), 1.f) *
        scale(1.0 / (screenWindow.width()),
              1.0 / (screenWindow.height()), 1.0) *
        translate(Length3D(-screenWindow.x() * 1.0_um, -screenWindow.y() * 1.0_um, 0.0 * 1.0_um)); // TODO: Should be dimensionless?
    RasterToScreen = inverse(ScreenToRaster);
    RasterToCamera = inverse(CameraToScreen) * RasterToScreen;
}


} // namespace
