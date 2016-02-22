
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


// core/sampler.cpp*
#include <iostream>
#include <numeric>

#include "stdafx.h"
#include "sampler.h"
//#include "integrator.h"
//#include "volume.h"
namespace photonflow {

// Sampler Method Definitions
Sampler::~Sampler() {
}


Sampler::Sampler(int xstart, int xend, int ystart, int yend, int spp,
                 photonflow::Time sopen, photonflow::Time sclose)
    : xPixelStart(xstart), xPixelEnd(xend), yPixelStart(ystart),
      yPixelEnd(yend), samplesPerPixel(spp), shutterOpen(sopen),
      shutterClose(sclose) { }
bool Sampler::reportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count) {
    UNUSED(samples);
    UNUSED(rays);
    UNUSED(Ls);
    UNUSED(isects);
    UNUSED(count);
    return true;
}


void Sampler::computeSubWindow(int num, int count, int *newXStart,
        int *newXEnd, int *newYStart, int *newYEnd) const {
    // Determine how many tiles to use in each dimension, _nx_ and _ny_
    int dx = xPixelEnd - xPixelStart, dy = yPixelEnd - yPixelStart;
    int nx = count, ny = 1;
    while ((nx & 0x1) == 0 && 2 * dx * ny < dy * nx) {
        nx >>= 1;
        ny <<= 1;
    }
    photonflowAssert(nx * ny == count);

    // Compute $x$ and $y$ pixel sample range for sub-window
    int xo = num % nx, yo = num / nx;
    double tx0 = double(xo) / double(nx), tx1 = double(xo+1) / double(nx);
    double ty0 = double(yo) / double(ny), ty1 = double(yo+1) / double(ny);
    *newXStart = Floor2Int(lerp(tx0, xPixelStart, xPixelEnd));
    *newXEnd   = Floor2Int(lerp(tx1, xPixelStart, xPixelEnd));
    *newYStart = Floor2Int(lerp(ty0, yPixelStart, yPixelEnd));
    *newYEnd   = Floor2Int(lerp(ty1, yPixelStart, yPixelEnd));
}



// Sample Method Definitions
Sample::Sample() {
}

uint32_t Sample::Add1D(uint32_t num) {
    SampleData d;
    d.size = num;
    if(m_metaData.size() > 0) {
        d.offset = m_metaData.back().offset + m_metaData.back().size;
    } else {
        d.offset = 0;
    }
    m_data.resize(d.offset + d.size);
    m_metaData.push_back(d);
    return d.offset;
}

Sample *Sample::Duplicate(int count) const {
    Sample *ret = new Sample[count];
    for (int i = 0; i < count; ++i) {
        ret[i].m_data = m_data;
        ret[i].m_metaData = m_metaData;
    }
    return ret;
}

size_t Sample::n1Dsize()
{
    return m_metaData.size();
}


size_t Sample::n1Dsize(int i)
{
    return m_metaData[i].size;
}

double &Sample::get(int i, int j)
{
    return m_data[m_metaData[i].offset + j];
}


} // namespace
