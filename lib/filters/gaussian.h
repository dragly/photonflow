
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

#ifndef PBRT_FILTERS_GAUSSIAN_H
#define PBRT_FILTERS_GAUSSIAN_H

// filters/gaussian.h*
#include "../core/filter.h"

// Gaussian Filter Declarations
class GaussianFilter : public Filter {
public:
    // GaussianFilter Public Methods
    GaussianFilter(double xw, double yw, double a)
        : Filter(xw, yw), alpha(a), expX(expf(-alpha * xWidth * xWidth)),
          expY(expf(-alpha * yWidth * yWidth)) { }
    double evaluate(double x, double y) const;
private:
    // GaussianFilter Private Data
    const double alpha;
    const double expX, expY;

    // GaussianFilter Utility Functions
    double Gaussian(double d, double expv) const {
        return max(0.0, double(expf(-alpha * d * d) - expv));
    }
};


//GaussianFilter *CreateGaussianFilter(const ParamSet &ps);

#endif // PBRT_FILTERS_GAUSSIAN_H
