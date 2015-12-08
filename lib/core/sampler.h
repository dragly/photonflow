
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

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "../core/common.h"
#include "geometry.h"
#include "../core/randomnumbergenerator.h"
#include "memory.h"

class Sample;
class RGBSpectrum;
using Spectrum=RGBSpectrum;
class Intersection;

// Sampling Declarations
class Sampler {
public:
    // Sampler Interface
    virtual ~Sampler();
    Sampler(int xstart, int xend, int ystart, int yend,
            int spp, double sopen, double sclose);
    virtual int GetMoreSamples(Sample *sample, RNG &rng) = 0;
    virtual int MaximumSampleCount() = 0;
    virtual bool ReportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count);
    virtual Sampler *GetSubSampler(int num, int count) = 0;
    virtual int RoundSize(int size) const = 0;

    // Sampler Public Data
    const int xPixelStart, xPixelEnd, yPixelStart, yPixelEnd;
    const int samplesPerPixel;
    const double shutterOpen, shutterClose;
protected:
    // Sampler Protected Methods
    void ComputeSubWindow(int num, int count, int *xstart, int *xend, int *ystart, int *yend) const;
};


struct CameraSample {
    double imageX, imageY;
    double lensU, lensV;
    double time;
};

struct SampleData {
    int size;
    int offset;
};

struct Sample : public CameraSample {
    // Sample Public Methods
    Sample();
    uint32_t Add1D(uint32_t num);
    Sample *Duplicate(int count) const;
    size_t n1Dsize();
    size_t n1Dsize(int i);
    double& get(int i, int j);
private:
    // Sample Private Methods
//    std::vector<uint32_t> n1D;
    std::vector<SampleData> m_metaData;
    std::vector<double> m_data;
};



#endif // PBRT_CORE_SAMPLER_H
