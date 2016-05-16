
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

namespace photonflow {

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
    virtual BoundingBox worldBound() const = 0;
    virtual bool intersectP(const Ray &ray, double *t0, double *t1) const = 0;
    virtual Spectrum absorption(const Point3D &, const Length3D &,
                             double time) const = 0;
    virtual Spectrum scattering(const Point3D &, const Length3D &,
                             double time) const = 0;
    virtual Spectrum emission(const Point3D &, const Length3D &,
                         double time) const = 0;
    virtual double phase(const Point3D &, const Length3D &,
                    const Length3D &, double time) const = 0;
    virtual Spectrum extinction(const Point3D &p, const Length3D &wo, double time) const;
//    virtual Spectrum tau(const Ray &ray, double step = 1.f,
//                         double offset = 0.5) const = 0;
};


class DensityRegion : public VolumeRegion {
public:
    // DensityRegion Public Methods
    DensityRegion();
    DensityRegion(const Spectrum &sa, const Spectrum &ss, double gg,
                  const Spectrum &emita, const Transform &VolumeToWorldIn)
        : m_absorptionCoefficient(sa)
        , m_scatteringCoefficient(ss)
        , m_emissionCoefficient(emita)
        , m_henyeyGreensteinFactor(gg)
        , WorldToVolume(inverse(VolumeToWorldIn))
        , VolumeToWorld(VolumeToWorldIn) { }
    virtual double Density(const Point3D &Pobj) const = 0;
    Spectrum absorption(const Point3D &p, const Length3D &, double) const {
        UNUSED(p);
        return m_absorptionCoefficient;
    }
    virtual Spectrum scattering(const Point3D &p, const Length3D &, double) const override {
        return Density(p) * m_scatteringCoefficient;
    }
    virtual Spectrum extinction(const Point3D &p, const Length3D &, double) const override {
        return Density(p) * (m_absorptionCoefficient + m_scatteringCoefficient);
    }
    virtual Spectrum emission(const Point3D &p, const Length3D &, double) const override {
        return Density(p) * m_emissionCoefficient;
    }
    virtual double phase(const Point3D &phase, const Length3D &w, const Length3D &wp, double) const override {
        UNUSED(phase);
        return PhaseHG(w, wp, m_henyeyGreensteinFactor);
    }
    double henyeyGreensteinFactor() {
        return m_henyeyGreensteinFactor;
    }

    void setAbsorptionCoefficient(const Spectrum &value) {
        m_absorptionCoefficient = value;
    }
    void setScatteringCoefficient(const Spectrum &value) {
        m_scatteringCoefficient = value;
    }
    void setEmissionCoefficient(const Spectrum &value) {
        m_emissionCoefficient = value;
    }
    void setHenyeyGreensteinFactor(double value) {
        m_henyeyGreensteinFactor = value;
    }

//    Spectrum tau(const Ray &r, double stepSize, double offset) const;
protected:
    // DensityRegion Protected Data
    Spectrum m_absorptionCoefficient;
    Spectrum m_scatteringCoefficient;
    Spectrum m_emissionCoefficient;
    double m_henyeyGreensteinFactor = 0.0;
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

} // namespace

#endif // PBRT_CORE_VOLUME_H
