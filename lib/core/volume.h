
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

#ifndef PBRT_CORE_VOLUME_H
#define PBRT_CORE_VOLUME_H

// core/volume.h*
#include "../core/common.h"
#include "../core/spectrum.h"
#include "../core/geometry.h"
#include "../core/transform.h"
#include "../core/integrator.h"

// Volume Scattering Declarations
double PhaseIsotropic(const Length3D &w, const Length3D &wp);
double PhaseRayleigh(const Length3D &w, const Length3D &wp);
double PhaseMieHazy(const Length3D &w, const Length3D &wp);
double PhaseMieMurky(const Length3D &w, const Length3D &wp);
double PhaseHG(const Length3D &w, const Length3D &wp, double g);
double PhaseSchlick(const Length3D &w, const Length3D &wp, double g);
class VolumeRegion {
public:
    // VolumeRegion Interface
    virtual BBox worldBound() const = 0;
    virtual bool intersectP(const Ray &ray, double *t0, double *t1) const = 0;
    virtual Spectrum sigma_a(const Point3D &, const Length3D &,
                             double time) const = 0;
    virtual Spectrum sigma_s(const Point3D &, const Length3D &,
                             double time) const = 0;
    virtual Spectrum Lve(const Point3D &, const Length3D &,
                         double time) const = 0;
    virtual double p(const Point3D &, const Length3D &,
                    const Length3D &, double time) const = 0;
    virtual Spectrum sigma_t(const Point3D &p, const Length3D &wo, double time) const;
//    virtual Spectrum tau(const Ray &ray, double step = 1.f,
//                         double offset = 0.5) const = 0;
};


class DensityRegion : public VolumeRegion {
public:
    // DensityRegion Public Methods
    DensityRegion();
    DensityRegion(const Spectrum &sa, const Spectrum &ss, double gg,
                  const Spectrum &emita, const Transform &VolumeToWorldIn)
        : sig_a(sa), sig_s(ss), le(emita), g(gg),
          WorldToVolume(Inverse(VolumeToWorldIn)),
          VolumeToWorld(VolumeToWorldIn) { }
    virtual double Density(const Point3D &Pobj) const = 0;
    Spectrum sigma_a(const Point3D &p, const Length3D &, double) const {
        UNUSED(p);
        return sig_a;
    }
    Spectrum sigma_s(const Point3D &p, const Length3D &, double) const {
        return Density(p) * sig_s;
    }
    Spectrum sigma_t(const Point3D &p, const Length3D &, double) const {
        return Density(p) * (sig_a + sig_s);
    }
    Spectrum Lve(const Point3D &p, const Length3D &, double) const {
        return Density(p) * le;
    }
    double p(const Point3D &p, const Length3D &w, const Length3D &wp, double) const {
        UNUSED(p);
        return PhaseHG(w, wp, g);
    }
//    Spectrum tau(const Ray &r, double stepSize, double offset) const;
protected:
    // DensityRegion Protected Data
    Spectrum sig_a;
    Spectrum sig_s;
    Spectrum le;
    double g = 0.0;
    Transform WorldToVolume;
    Transform VolumeToWorld;
};


//class AggregateVolume : public VolumeRegion {
//public:
//    // AggregateVolume Public Methods
//    AggregateVolume(const vector<VolumeRegion *> &r);
//    ~AggregateVolume();
//    BBox worldBound() const;
//    bool intersectP(const Ray &ray, double *t0, double *t1) const;
//    Spectrum sigma_a(const Point3D &, const Vector3D &, double) const;
//    Spectrum sigma_s(const Point3D &, const Vector3D &, double) const;
//    Spectrum Lve(const Point3D &, const Vector3D &, double) const;
//    double p(const Point3D &, const Vector3D &, const Vector3D &, double) const;
//    Spectrum sigma_t(const Point3D &, const Vector3D &, double) const;
//    Spectrum tau(const Ray &ray, double, double) const;
//private:
//    // AggregateVolume Private Data
//    vector<VolumeRegion *> regions;
//    BBox bound;
//};


bool volumeScatteringProperties(const std::string &name, Spectrum *sigma_a,
                                   Spectrum *sigma_prime_s);
//class VolumeIntegrator : public Integrator {
//public:
//    // VolumeIntegrator Interface
//    virtual Spectrum Li(const Scene *scene, const Renderer *renderer,
//        const RayDifferential &ray, const Sample *sample, RNG &rng,
//        Spectrum *transmittance) const = 0;
//    virtual Spectrum Transmittance(const Scene *scene,
//        const Renderer *renderer, const RayDifferential &ray,
//        const Sample *sample, RNG &rng) const = 0;
//};


void subsurfaceFromDiffuse(const Spectrum &Kd, double meanPathLength, double eta,
                           Spectrum *sigma_a, Spectrum *sigma_prime_s);

#endif // PBRT_CORE_VOLUME_H
