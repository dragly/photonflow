
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

#ifndef PBRT_VOLUMES_VOLUMEGRID_H
#define PBRT_VOLUMES_VOLUMEGRID_H

// volumes/volumegrid.h*
#include "../core/volume.h"
#include "armadillo_includer.h"
#include <memory.h>

namespace photonflow {

// VolumeGridDensity Declarations
class VolumeGridDensity : public DensityRegion {
public:
    // VolumeGridDensity Public Methods
    VolumeGridDensity();
    VolumeGridDensity(const Spectrum &sa, const Spectrum &ss, double gg,
            const Spectrum &emita, const BoundingBox &e, const Transform &v2w,
            arma::Cube<short> densitya);
    BoundingBox worldBound() const { return m_worldBound; }
    bool intersectP(const Ray &r, double *t0, double *t1) const;
    double Density(const Point3D &Pobj) const;
    double D(int x, int y, int z) const;
    bool inside(const Point3D &p) const;
    bool fuzzyInside(const Point3D &p) const;
private:
    // VolumeGridDensity Private Data
    BoundingBox extent;
    arma::Cube<short> density;
    BoundingBox m_worldBound;
};

inline bool VolumeGridDensity::intersectP(const Ray &r, double *t0, double *t1) const {
    Ray ray = WorldToVolume(r);
    return extent.intersectP(ray, t0, t1);
}

inline double VolumeGridDensity::D(int x, int y, int z) const {
    x = clamp(x, 0, int(density.n_rows-1));
    y = clamp(y, 0, int(density.n_cols-1));
    z = clamp(z, 0, int(density.n_slices-1));
    double value = 0.0;
    try {
        value = density(x, y, z);
    } catch(std::logic_error &e) {
        std::cout << "Error on fetching value for " << x << " " << y << " " << z << std::endl;
        std::cout << "Size is " << density.n_rows << " " << density.n_cols << " " << density.n_slices << std::endl;
    }

    return value;
}


//VolumeGridDensity *CreateGridVolumeRegion(const Transform &volume2world,
//        const ParamSet &params);

} // namespace

#endif // PBRT_VOLUMES_VOLUMEGRID_H
