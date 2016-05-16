
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

#include "integrator.h"
#include "../volumes/volumegrid.h"

#include <iostream>

using namespace std;

namespace photonflow {

Integrator::Integrator(VolumeGridDensity *volumeGridDensity, Ray startRay, int bounces, RNG &rng)
    : m_volumeGridDensity(volumeGridDensity)
    , m_ray(startRay)
    , m_bounces(bounces)
    , m_rng(&rng)
{
}

void Integrator::integrate(std::function<Control(const Ray& ray, photonflow::Length stepLength)> callback)
{
    for(int i = 0; i < m_bounces; i++) {
        double eps = 1e-16; // avoid -log(0) which returns inf
        photonflow::Length ds = 10.0_um * -log(m_rng->randomFloat() + eps);
        double g = m_volumeGridDensity->henyeyGreensteinFactor();

        double cosTheta = Distribution::heyneyGreenstein(g, *m_rng);
        double sinTheta = sqrt(1 - cosTheta*cosTheta);
        double phi = 2.0 * M_PI * m_rng->randomFloat();

        Vector3D perpendicular = m_ray.direction().perpendicular().normalized();
        Transform phiRotation = rotate(phi, m_ray.direction().normalized());
        perpendicular = phiRotation(perpendicular);

        Transform directionRotation = rotate(cosTheta, sinTheta, perpendicular);

        Length3D direction = directionRotation(m_ray.direction());
        direction = direction.normalized() * ds;

        Point3D origin = m_ray.origin() + direction;
        m_ray = Ray(origin, direction);
        Control result = callback(m_ray, ds);
        if(result == Control::Break) {
            break;
        }
    }
}
} // namespace
