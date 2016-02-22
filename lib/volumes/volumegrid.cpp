
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

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


// volumes/volumegrid.cpp*
#include "stdafx.h"
#include "volumegrid.h"
//#include "paramset.h"

using namespace arma;

VolumeGridDensity::VolumeGridDensity()
{

}

VolumeGridDensity::VolumeGridDensity(const Spectrum &sa, const Spectrum &ss, double gg, const Spectrum &emita, const BoundingBox &e, const Transform &v2w, arma::Cube<short> densitya)
    : DensityRegion(sa, ss, gg, emita, v2w)
    , extent(e)
    , m_worldBound(VolumeToWorld(extent))
{
    density = densitya;
}

// VolumeGridDensity Method Definitions
double VolumeGridDensity::Density(const Point3D &Pobj) const {
    Point3D local = WorldToVolume(Pobj);

    if (!extent.inside(local)) {
        return 0;
    }

    // Compute voxel coordinates and offsets for _Pobj_
    Length3D vox = extent.offset(local);
    vox.x = vox.x * double(density.n_rows); // - .5f;
    vox.y = vox.y * double(density.n_cols); // - .5f;
    vox.z = vox.z * double(density.n_slices); // - .5f;
    int vx = Floor2Int(vox.x.value());
    int vy = Floor2Int(vox.y.value());
    int vz = Floor2Int(vox.z.value()); // TODO: Is this proper with dimensions?
//    double dx = vox.x - vx;
//    double dy = vox.y - vy;
//    double dz = vox.z - vz;

    return D(vx, vy, vz);

    // Trilinearly interpolate density values to compute local density
//    double d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
//    double d10 = Lerp(dx, D(vx, vy+1, vz),   D(vx+1, vy+1, vz));
//    double d01 = Lerp(dx, D(vx, vy, vz+1),   D(vx+1, vy, vz+1));
//    double d11 = Lerp(dx, D(vx, vy+1, vz+1), D(vx+1, vy+1, vz+1));
//    double d0 = Lerp(dy, d00, d10);
//    double d1 = Lerp(dy, d01, d11);
    //    return Lerp(dz, d0, d1);
}

bool VolumeGridDensity::inside(const Point3D &p) const
{
    Point3D local = WorldToVolume(p);
    return extent.inside(local);
}

bool VolumeGridDensity::fuzzyInside(const Point3D &p) const
{
    Point3D local = WorldToVolume(p);
    return extent.fuzzyInside(local);
}

