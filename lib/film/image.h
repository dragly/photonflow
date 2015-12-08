
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

#ifndef PBRT_FILM_IMAGE_H
#define PBRT_FILM_IMAGE_H

// film/image.h*
#include "../core/common.h"
#include "../core/film.h"
#include "../core/sampler.h"
#include "../core/filter.h"
#include "../core/memory.h"
//#include "paramset.h"


class Pixel {
public:
    Pixel() {
        for (int i = 0; i < 3; ++i) Lxyz[i] = splatXYZ[i] = 0.0;
        weightSum = 0.0;
    }
    double Lxyz[3];
    double weightSum;
    double splatXYZ[3];
    double pad;
};

// ImageFilm Declarations
class ImageFilm : public Film {
public:
    // ImageFilm Public Methods
    ImageFilm(int xres, int yres, Filter *filt, const double crop[4]);
    ~ImageFilm() {
        delete pixels;
//        delete filter;
        delete[] filterTable;
    }
    void addSample(const CameraSample &sample, const Spectrum &L);
    void splat(const CameraSample &sample, const Spectrum &L);
    void sampleExtent(int *xstart, int *xend, int *ystart, int *yend) const;
    void pixelExtent(int *xstart, int *xend, int *ystart, int *yend) const;
    void writeImage(double splatScale);
    void updateDisplay(int x0, int y0, int x1, int y1, double splatScale);
//private:
    // ImageFilm Private Data
    Filter *filter;
    double cropWindow[4];
    int xPixelStart, yPixelStart, xPixelCount, yPixelCount;
    BlockedArray<Pixel> *pixels;
    double *filterTable;
};


//ImageFilm *CreateImageFilm(const ParamSet &params, Filter *filter);

#endif // PBRT_FILM_IMAGE_H
