
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

#ifndef PBRT_CORE_MONTECARLO_H
#define PBRT_CORE_MONTECARLO_H

// core/montecarlo.h*
#include "../core/common.h"
#include "geometry.h"
#include "../core/randomnumbergenerator.h"
#include <memory.h>
#include <vector>

using std::vector;

// smallest floating point value less than one; all canonical random samples
// should be <= this.
#ifdef PBRT_IS_WINDOWS
// sadly, MSVC2008 (at least) doesn't support hexidecimal fp constants...
static const double OneMinusEpsilon=0.9999999403953552f;
#else
static const double OneMinusEpsilon=0x1.fffffep-1;
#endif

// Monte Carlo Utility Declarations
struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const double *f, int n) {
        count = n;
        func = new double[n];
        memcpy(func, f, n*sizeof(double));
        cdf = new double[n+1];
        // Compute integral of step function at $x_i$
        cdf[0] = 0.;
        for (int i = 1; i < count+1; ++i)
            cdf[i] = cdf[i-1] + func[i-1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[count];
        if (funcInt == 0.0) {
            for (int i = 1; i < n+1; ++i)
                cdf[i] = double(i) / double(n);
        }
        else {
            for (int i = 1; i < n+1; ++i)
                cdf[i] /= funcInt;
        }
    }
    ~Distribution1D() {
        delete[] func;
        delete[] cdf;
    }
    double sampleContinuous(double u, double *pdf, int *off = NULL) const {
        // Find surrounding CDF segments and _offset_
        double *ptr = std::upper_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));
        if (off) *off = offset;
        photonflowAssert(offset < count);
        photonflowAssert(u >= cdf[offset] && u < cdf[offset+1]);

        // Compute offset along CDF segment
        double du = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
        photonflowAssert(!std::isnan(du));

        // Compute PDF for sampled offset
        if (pdf) *pdf = func[offset] / funcInt;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / count;
    }
    int sampleDiscrete(double u, double *pdf) const {
        // Find surrounding CDF segments and _offset_
        double *ptr = std::upper_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));
        photonflowAssert(offset < count);
        photonflowAssert(u >= cdf[offset] && u < cdf[offset+1]);
        if (pdf) *pdf = func[offset] / (funcInt * count);
        return offset;
    }
private:
    friend struct Distribution2D;
    // Distribution1D Private Data
    double *func, *cdf;
    double funcInt;
    int count;
};


void rejectionSampleDisk(double *x, double *y, RNG &rng);
Length3D uniformSampleHemisphere(double u1, double u2);
double  uniformHemispherePdf();
Length3D uniformSampleSphere(double u1, double u2);
double  uniformSpherePdf();
Length3D uniformSampleCone(double u1, double u2, double thetamax);
Length3D uniformSampleCone(double u1, double u2, double thetamax,
    const Length3D &x, const Length3D &y, const Length3D &z);
double  uniformConePdf(double thetamax);
void uniformSampleDisk(double u1, double u2, double *x, double *y);
void concentricSampleDisk(double u1, double u2, double *dx, double *dy);
inline Length3D cosineSampleHemisphere(double u1, double u2) {
    double x;
    double y;
    concentricSampleDisk(u1, u2, &x, &y);
    Length3D ret(x*1.0_um, y*1.0_um, 0.0_um);
    ret.z = sqrt(max(0.0_um*0.0_um, 1.0_um*1.0_um - ret.x*ret.x - ret.y*ret.y));
    return ret;
}


inline double cosineHemispherePdf(double costheta, double phi) {
    UNUSED((phi));
    return costheta * INV_PI;
}


void uniformSampleTriangle(double ud1, double ud2, double *u, double *v);
struct Distribution2D {
    // Distribution2D Public Methods
    Distribution2D(const double *data, int nu, int nv);
    ~Distribution2D();
    void SampleContinuous(double u0, double u1, double uv[2],
                          double *pdf) const {
        double pdfs[2];
        int v;
        uv[1] = pMarginal->sampleContinuous(u1, &pdfs[1], &v);
        uv[0] = pConditionalV[v]->sampleContinuous(u0, &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
    }
    double Pdf(double u, double v) const {
        int iu = clamp(Float2Int(u * pConditionalV[0]->count), 0,
                       pConditionalV[0]->count-1);
        int iv = clamp(Float2Int(v * pMarginal->count), 0,
                       pMarginal->count-1);
        if (pConditionalV[iv]->funcInt * pMarginal->funcInt == 0.0) return 0.0;
        return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) /
               (pConditionalV[iv]->funcInt * pMarginal->funcInt);
    }
private:
    // Distribution2D Private Data
    vector<Distribution1D *> pConditionalV;
    Distribution1D *pMarginal;
};


void stratifiedSample1D(double *samples, int nsamples, RNG &rng,
                        bool jitter = true);
void stratifiedSample2D(double *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
template <typename T>
void shuffle(T *samp, uint32_t count, uint32_t dims, RNG &rng) {
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t other = i + (rng.randomUInt() % (count - i));
        for (uint32_t j = 0; j < dims; ++j)
            swap(samp[dims*i + j], samp[dims*other + j]);
    }
}


void latinHypercube(double *samples, uint32_t nSamples, uint32_t nDim, RNG &rng);
inline double RadicalInverse(int n, int base) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;
    while (n > 0) {
        // Compute next digit of radical inverse
        int d_i = (n % base);
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


inline void generatePermutation(uint32_t *buf, uint32_t b, RNG &rng) {
    for (uint32_t i = 0; i < b; ++i)
        buf[i] = i;
    shuffle(buf, b, 1, rng);
}


inline double permutedRadicalInverse(uint32_t n, uint32_t base,
                                     const uint32_t *p) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;

    while (n > 0) {
        uint32_t d_i = p[n % base];
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


class PermutedHalton {
public:
    // PermutedHalton Public Methods
    PermutedHalton(uint32_t d, RNG &rng);
    ~PermutedHalton() {
        delete[] b;
        delete[] permute;
    }
    void Sample(uint32_t n, double *out) const {
        uint32_t *p = permute;
        for (uint32_t i = 0; i < dims; ++i) {
            out[i] = min(double(permutedRadicalInverse(n, b[i], p)),
                         OneMinusEpsilon);
            p += b[i];
        }
    }
private:
    // PermutedHalton Private Data
    uint32_t dims;
    uint32_t *b, *permute;
    PermutedHalton(const PermutedHalton &);
    PermutedHalton &operator=(const PermutedHalton &);
};


inline double vanDerCorput(uint32_t n, uint32_t scramble = 0);
inline double sobol2(uint32_t n, uint32_t scramble = 0);
inline double larcherPillichshammer2(uint32_t n, uint32_t scramble = 0);
inline void sample02(uint32_t n, const uint32_t scramble[2], double sample[2]);
//int LDPixelSampleFloatsNeeded(const Sample *sample, int nPixelSamples);
//void LDPixelSample(int xPos, int yPos, double shutterOpen,
//    double shutterClose, int nPixelSamples, Sample *samples, double *buf, RandomNumberGenerator &rng);
Length3D sampleHG(const Length3D &w, double g, double u1, double u2);
double HGPdf(const Length3D &w, const Length3D &wp, double g);

// Monte Carlo Inline Functions
inline double balanceHeuristic(int nf, double fPdf, int ng, double gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


inline double powerHeuristic(int nf, double fPdf, int ng, double gPdf) {
    double f = nf * fPdf, g = ng * gPdf;
    return (f*f) / (f*f + g*g);
}



// Sampling Inline Functions
inline void sample02(uint32_t n, const uint32_t scramble[2],
                     double sample[2]) {
    sample[0] = vanDerCorput(n, scramble[0]);
    sample[1] = sobol2(n, scramble[1]);
}


inline double vanDerCorput(uint32_t n, uint32_t scramble) {
    // Reverse bits of _n_
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n ^= scramble;
    return min(((n>>8) & 0xffffff) / double(1 << 24), OneMinusEpsilon);
}


inline double sobol2(uint32_t n, uint32_t scramble) {
    for (uint32_t v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
        if (n & 0x1) scramble ^= v;
    return min(((scramble>>8) & 0xffffff) / double(1 << 24), OneMinusEpsilon);
}


inline double
larcherPillichshammer2(uint32_t n, uint32_t scramble) {
    for (uint32_t v = 1 << 31; n != 0; n >>= 1, v |= v >> 1)
        if (n & 0x1) scramble ^= v;
    return min(((scramble>>8) & 0xffffff) / double(1 << 24), OneMinusEpsilon);
}


inline void LDShuffleScrambled1D(int nSamples, int nPixel,
                                 double *samples, RNG &rng) {
    uint32_t scramble = rng.randomUInt();
    for (int i = 0; i < nSamples * nPixel; ++i)
        samples[i] = vanDerCorput(i, scramble);
    for (int i = 0; i < nPixel; ++i)
        shuffle(samples + i * nSamples, nSamples, 1, rng);
    shuffle(samples, nPixel, nSamples, rng);
}


inline void LDShuffleScrambled2D(int nSamples, int nPixel,
                                 double *samples, RNG &rng) {
    uint32_t scramble[2] = { rng.randomUInt(), rng.randomUInt() };
    for (int i = 0; i < nSamples * nPixel; ++i)
        sample02(i, scramble, &samples[2*i]);
    for (int i = 0; i < nPixel; ++i)
        shuffle(samples + 2 * i * nSamples, nSamples, 2, rng);
    shuffle(samples, nPixel, 2 * nSamples, rng);
}



#endif // PBRT_CORE_MONTECARLO_H
