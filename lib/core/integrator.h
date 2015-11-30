
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

class Integrator
{
public:
    class iterator {
    public:
        typedef Ray value_type;
        typedef Ray& reference_type;
        typedef std::input_iterator_tag iterator_category;

        iterator(int bounce)
        {
            m_bounce = bounce;
        }

        iterator(Ray ray, RNG *rng)
            : m_ray(ray)
            , m_rng(rng)
        {
//            next();
        }

        bool operator==(const iterator& other) const {
            return other.m_bounce == m_bounce;
        }
        bool operator!=(const iterator& other) const {
            return !(other == *this);
        }

        iterator& next() {
            double dt = 0.01;
            double g = 0.99;

            double cosTheta = Distribution::heyneyGreenstein(g, *m_rng);
            double sinTheta = sqrt(1 - cosTheta*cosTheta);
            double phi = 2.0 * M_PI * m_rng->RandomFloat();

            Vector perpendicular = m_ray.direction().perpendicular();
            Transform phiRotation = Rotate(phi, m_ray.direction());
            perpendicular = phiRotation(perpendicular);

            Transform directionRotation = Rotatec(cosTheta, sinTheta, perpendicular);

            Vector direction = directionRotation(m_ray.direction());
            direction = direction.normalized();

            Point origin = m_ray.origin() + direction * dt;
            m_ray = Ray(origin, direction);

            m_bounce += 1;

            return *this;
        }

        iterator& operator++() {
            return next();
        }

        Ray& operator*() {
            return m_ray;
        }

        int m_bounce = 0;
        Ray m_ray;
        RNG *m_rng;

        bool m_final = false;
    };

    Integrator(Ray startRay, int bounces, RNG &rng)
        : m_startRay(startRay)
        , m_bounces(bounces)
        , m_rng(&rng)
    {
    }

    iterator begin() {
        return iterator(m_startRay, m_rng);
    }

    iterator end() {
        return iterator(m_bounces);
    }

private:
    Ray m_startRay;
    int m_bounces;
    RNG *m_rng;
};

#endif // PBRT_CORE_INTEGRATOR_H
