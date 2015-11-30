
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

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

#include <functional>
#include "../core/geometry.h"
#include "../core/transform.h"
#include "../core/heyneygreenstein.h"

class VolumeGridDensity;

class Integrator
{
public:
    class iterator {
    public:
        typedef Ray value_type;
        typedef Ray& reference_type;
        typedef std::input_iterator_tag iterator_category;

        iterator(Integrator* parent, int bounce)
            : m_parent(parent)
            , m_bounce(bounce)
        {
        }

        iterator(Integrator* parent)
            : m_parent(parent)
        {
        }

        bool operator==(const iterator& other) const {
            return other.m_bounce == m_bounce;
        }

        bool operator!=(const iterator& other) const {
            return !(other == *this);
        }

        iterator& operator++() {
            m_bounce++;
            m_parent->next();
            return *this;
        }

        Ray& operator*() {
            return m_parent->m_ray;
        }

        Integrator* m_parent;

        int m_bounce = 0;
    };

    Integrator(VolumeGridDensity *volumeGridDensity, Ray startRay, int bounces, RNG &rng)
        : m_volumeGridDensity(volumeGridDensity)
        , m_ray(startRay)
        , m_bounces(bounces)
        , m_rng(&rng)
    {
    }

    iterator begin() {
        return iterator(this);
    }

    iterator end() {
        return iterator(this, m_bounces);
    }

    void next();

private:
    VolumeGridDensity *m_volumeGridDensity;
    Ray m_ray;
    int m_bounces;
    RNG *m_rng;
};

#endif // PBRT_CORE_INTEGRATOR_H
