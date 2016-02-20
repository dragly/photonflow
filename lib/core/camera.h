
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

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "common.h"
#include "geometry.h"
#include "transform.h"
#include "../geometry/rectangle.h"

#include <memory>

class Film;
class CameraSample;

// Camera Declarations
class Camera {
public:
    // Camera Interface
    Camera(const Transform &cam2world, boost::units::photonflow::Time sopen, boost::units::photonflow::Time sclose,
           std::shared_ptr<Film> film);
    virtual ~Camera();
    virtual double generateRay(const CameraSample &sample,
                              Ray *ray) const = 0;
//    virtual double generateRayDifferential(const CameraSample &sample, RayDifferential *rd) const;

    // Camera Public Data
    Transform CameraToWorld;
    const boost::units::photonflow::Time shutterOpen, shutterClose;
    std::shared_ptr<Film> film;
};


class ProjectiveCamera : public Camera {
public:
    // ProjectiveCamera Public Methods
    ProjectiveCamera(const Transform &cam2world,
                     const Transform &proj,
                     const Rectangle &screenWindow,
                     boost::units::photonflow::Time sopen,
                     boost::units::photonflow::Time sclose,
                     boost::units::photonflow::Length lensr,
                     boost::units::photonflow::Length focald,
                     std::shared_ptr<Film> film);
protected:
    // ProjectiveCamera Protected Data
    Transform CameraToScreen, RasterToCamera;
    Transform ScreenToRaster, RasterToScreen;
    boost::units::photonflow::Length lensRadius, focalDistance;
};



#endif // PBRT_CORE_CAMERA_H
