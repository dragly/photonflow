
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

#ifndef PBRT_CAMERAS_PERSPECTIVE_H
#define PBRT_CAMERAS_PERSPECTIVE_H

// cameras/perspective.h*
#include "../core/common.h"
#include "../core/camera.h"
#include "../core/film.h"

#include <memory>

class ParamSet;

namespace photonflow {

// PerspectiveCamera Declarations
class PerspectiveCamera : public ProjectiveCamera {
public:
    // PerspectiveCamera Public Methods
    PerspectiveCamera(const Transform &cam2world,
        const Rectangle &screenWindow, photonflow::Time sopen, photonflow::Time sclose,
        photonflow::Length lensr, photonflow::Length focald, double fov, std::shared_ptr<Film> film);
    double generateRay(const CameraSample &sample, Ray *) const;
//    double generateRayDifferential(const CameraSample &sample,
//                                  RayDifferential *ray) const;
private:
    // PerspectiveCamera Private Data
    Length3D dxCamera, dyCamera;
};

}

#endif // PBRT_CAMERAS_PERSPECTIVE_H
