
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

#ifndef PBRT_CORE_SPECTRUM_H
#define PBRT_CORE_SPECTRUM_H

// core/spectrum.h*
#include "../core/common.h"
#include <vector>

using std::vector;
//#include "parallel.h"

namespace photonflow {

class RGBSpectrum;
using Spectrum=RGBSpectrum;

// Spectrum Utility Declarations
static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 30;
extern bool spectrumSamplesSorted(const double *lambda, const double *vals, int n);
extern void sortSpectrumSamples(double *lambda, double *vals, int n);
extern double averageSpectrumSamples(const double *lambda, const double *vals,
    int n, double lambdaStart, double lambdaEnd);
inline void XYZToRGB(const double xyz[3], double rgb[3]) {
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}


inline void RGBToXYZ(const double rgb[3], double xyz[3]) {
    xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
    xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
    xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}


enum SpectrumType { SPECTRUM_REFLECTANCE, SPECTRUM_ILLUMINANT };
extern void Blackbody(const double *wl, int n, double temp, double *vals);
extern double InterpolateSpectrumSamples(const double *lambda, const double *vals,
                                        int n, double l);

// Spectral Data Declarations
static const int nCIESamples = 471;
extern const double CIE_X[nCIESamples];
extern const double CIE_Y[nCIESamples];
extern const double CIE_Z[nCIESamples];
extern const double CIE_lambda[nCIESamples];
static const double CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;
extern const double RGB2SpectLambda[nRGB2SpectSamples];
extern const double RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const double RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const double RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const double RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const double RGBRefl2SpectRed[nRGB2SpectSamples];
extern const double RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const double RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const double RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const double RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const double RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const double RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const double RGBIllum2SpectRed[nRGB2SpectSamples];
extern const double RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const double RGBIllum2SpectBlue[nRGB2SpectSamples];

// Spectrum Declarations
template <int nSamples> class CoefficientSpectrum {
public:
    // CoefficientSpectrum Public Methods
    CoefficientSpectrum(double v = 0.0) {
        for (int i = 0; i < nSamples; ++i)
            c[i] = v;
        photonflowAssert(!hasNaNs());
    }
#ifdef DEBUG
    CoefficientSpectrum(const CoefficientSpectrum &s) {
        Assert(!s.hasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
    }
    
    CoefficientSpectrum &operator=(const CoefficientSpectrum &s) {
        Assert(!s.hasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
        return *this;
    }
#endif // DEBUG
    void print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < nSamples; ++i) {
            fprintf(f, "%f", c[i]);
            if (i != nSamples-1) fprintf(f, ", ");
        }
        fprintf(f, "]");
    }
    CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
        photonflowAssert(!s2.hasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] += s2.c[i];
        return *this;
    }
    CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
        photonflowAssert(!s2.hasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] += s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
        photonflowAssert(!s2.hasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] -= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
        photonflowAssert(!s2.hasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] /= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator*(const CoefficientSpectrum &sp) const {
        photonflowAssert(!sp.hasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] *= sp.c[i];
        return ret;
    }
    CoefficientSpectrum &operator*=(const CoefficientSpectrum &sp) {
        photonflowAssert(!sp.hasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] *= sp.c[i];
        return *this;
    }
    CoefficientSpectrum operator*(double a) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] *= a;
        photonflowAssert(!ret.hasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator*=(double a) {
        for (int i = 0; i < nSamples; ++i)
            c[i] *= a;
        photonflowAssert(!hasNaNs());
        return *this;
    }
    friend inline
    CoefficientSpectrum operator*(double a, const CoefficientSpectrum &s) {
        photonflowAssert(!std::isnan(a) && !s.hasNaNs());
        return s * a;
    }
    CoefficientSpectrum operator/(double a) const {
        photonflowAssert(!std::isnan(a));
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] /= a;
        photonflowAssert(!ret.hasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator/=(double a) {
        photonflowAssert(!std::isnan(a));
        for (int i = 0; i < nSamples; ++i)
            c[i] /= a;
        return *this;
    }
    bool operator==(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != sp.c[i]) return false;
        return true;
    }
    bool operator<(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] > sp.c[i]) return false;
        return true;
    }
    bool operator!=(const CoefficientSpectrum &sp) const {
        return !(*this == sp);
    }
    bool IsBlack() const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != 0.) return false;
        return true;
    }
    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = sqrtf(s.c[i]);
        photonflowAssert(!ret.hasNaNs());
        return ret;
    }
    template <int n> friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n> &s, double e);
    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = -c[i];
        return ret;
    }
    friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = expf(s.c[i]);
        photonflowAssert(!ret.hasNaNs());
        return ret;
    }
    CoefficientSpectrum clamp(double low = 0, double high = INFINITY) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = ::photonflow::clamp(c[i], low, high);
        photonflowAssert(!ret.hasNaNs());
        return ret;
    }
    bool hasNaNs() const {
        for (int i = 0; i < nSamples; ++i)
            if (std::isnan(c[i])) return true;
        return false;
    }
    bool Write(FILE *f) const {
        for (int i = 0; i < nSamples; ++i)
            if (fprintf(f, "%f ", c[i]) < 0) return false;
        return true;
    }
    bool Read(FILE *f) {
        for (int i = 0; i < nSamples; ++i)
            if (fscanf(f, "%f ", &c[i]) != 1) return false;
        return true;
    }
protected:
    // CoefficientSpectrum Protected Data
    double c[nSamples];
};


class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
public:
    // SampledSpectrum Public Methods
    SampledSpectrum(double v = 0.0) {
        for (int i = 0; i < nSpectralSamples; ++i) c[i] = v;
    }
    SampledSpectrum(const CoefficientSpectrum<nSpectralSamples> &v)
        : CoefficientSpectrum<nSpectralSamples>(v) { }
    static SampledSpectrum FromSampled(const double *lambda,
                                       const double *v, int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!spectrumSamplesSorted(lambda, v, n)) {
            vector<double> slambda(&lambda[0], &lambda[n]);
            vector<double> sv(&v[0], &v[n]);
            sortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        SampledSpectrum r;
        for (int i = 0; i < nSpectralSamples; ++i) {
            // Compute average value of given SPD over $i$th sample's range
            double lambda0 = lerp(double(i) / double(nSpectralSamples),
                                 sampledLambdaStart, sampledLambdaEnd);
            double lambda1 = lerp(double(i+1) / double(nSpectralSamples),
                                 sampledLambdaStart, sampledLambdaEnd);
            r.c[i] = averageSpectrumSamples(lambda, v, n, lambda0, lambda1);
        }
        return r;
    }
    static void Init() {
        // Compute XYZ matching functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            double wl0 = lerp(double(i) / double(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            double wl1 = lerp(double(i+1) / double(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            X.c[i] = averageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples,
                                            wl0, wl1);
            Y.c[i] = averageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples,
                                            wl0, wl1);
            Z.c[i] = averageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples,
                                            wl0, wl1);
        }

        // Compute RGB to spectrum functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            double wl0 = lerp(double(i) / double(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            double wl1 = lerp(double(i+1) / double(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            rgbRefl2SpectWhite.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectCyan.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectMagenta.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectYellow.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectRed.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectRed,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectGreen.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectBlue.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
                nRGB2SpectSamples, wl0, wl1);
        
            rgbIllum2SpectWhite.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectCyan.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectMagenta.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectYellow.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectRed.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectGreen.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectBlue.c[i] = averageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
                nRGB2SpectSamples, wl0, wl1);
        }
    }
    void ToXYZ(double xyz[3]) const {
        xyz[0] = xyz[1] = xyz[2] = 0.0;
        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += X.c[i] * c[i];
            xyz[1] += Y.c[i] * c[i];
            xyz[2] += Z.c[i] * c[i];
        }
        double scale = double(sampledLambdaEnd - sampledLambdaStart) /
            double(CIE_Y_integral * nSpectralSamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
    }
    double y() const {
        double yy = 0.0;
        for (int i = 0; i < nSpectralSamples; ++i)
            yy += Y.c[i] * c[i];
        return yy * double(sampledLambdaEnd - sampledLambdaStart) /
            double(CIE_Y_integral * nSpectralSamples);
    }
    void ToRGB(double rgb[3]) const {
        double xyz[3];
        ToXYZ(xyz);
        XYZToRGB(xyz, rgb);
    }
    RGBSpectrum ToRGBSpectrum() const;
    static SampledSpectrum FromRGB(const double rgb[3],
        SpectrumType type = SPECTRUM_REFLECTANCE);
    static SampledSpectrum FromXYZ(const double xyz[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        double rgb[3];
        XYZToRGB(xyz, rgb);
        return FromRGB(rgb, type);
    }
    SampledSpectrum(const RGBSpectrum &r, SpectrumType type = SPECTRUM_REFLECTANCE);
private:
    // SampledSpectrum Private Data
    static SampledSpectrum X, Y, Z;
    static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
    static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
    static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
    static SampledSpectrum rgbRefl2SpectBlue;
    static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
    static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
    static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
    static SampledSpectrum rgbIllum2SpectBlue;
};


class RGBSpectrum : public CoefficientSpectrum<3> {
    using CoefficientSpectrum<3>::c;
public:
    // RGBSpectrum Public Methods
    RGBSpectrum(double v = 0.0) : CoefficientSpectrum<3>(v) { }
    RGBSpectrum(const CoefficientSpectrum<3> &v)
        : CoefficientSpectrum<3>(v) { }
    RGBSpectrum(const RGBSpectrum &s)
        : CoefficientSpectrum<3>(s) {
        *this = s;
    }
    static RGBSpectrum fromRGB(const double rgb[3]) {
        RGBSpectrum s;
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        photonflowAssert(!s.hasNaNs());
        return s;
    }
    void toRGB(double *rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }
    const RGBSpectrum &toRGBSpectrum() const {
        return *this;
    }
    void roXYZ(double xyz[3]) const {
        RGBToXYZ(c, xyz);
    }
    static RGBSpectrum fromXYZ(const double xyz[3]) {
        RGBSpectrum r;
        XYZToRGB(xyz, r.c);
        return r;
    }
    double y() const {
        const double YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum fromSampled(const double *lambda, const double *v,
                                   int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!spectrumSamplesSorted(lambda, v, n)) {
            vector<double> slambda(&lambda[0], &lambda[n]);
            vector<double> sv(&v[0], &v[n]);
            sortSpectrumSamples(&slambda[0], &sv[0], n);
            return fromSampled(&slambda[0], &sv[0], n);
        }
        double xyz[3] = { 0, 0, 0 };
        double yint = 0.0;
        for (int i = 0; i < nCIESamples; ++i) {
            yint += CIE_Y[i];
            double val = InterpolateSpectrumSamples(lambda, v, n,
                                                   CIE_lambda[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }
        xyz[0] /= yint;
        xyz[1] /= yint;
        xyz[2] /= yint;
        return fromXYZ(xyz);
    }
};



// Spectrum Inline Functions
template <int nSamples> inline CoefficientSpectrum<nSamples>
pow(const CoefficientSpectrum<nSamples> &s, double e) {
    CoefficientSpectrum<nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i] = powf(s.c[i], e);
    photonflowAssert(!ret.hasNaNs());
    return ret;
}


inline Spectrum lerp(double t, const Spectrum &s1, const Spectrum &s2) {
    return (1.f - t) * s1 + t * s2;
}



} // namespace

#endif // PBRT_CORE_SPECTRUM_H
