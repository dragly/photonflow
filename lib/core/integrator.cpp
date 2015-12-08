
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

Integrator::Integrator(VolumeGridDensity *volumeGridDensity, Ray startRay, int bounces, RNG &rng)
    : m_volumeGridDensity(volumeGridDensity)
    , m_ray(startRay)
    , m_bounces(bounces)
    , m_rng(&rng)
{
}

Integrator::iterator Integrator::begin() {
    return iterator(this);
}

Integrator::iterator Integrator::end() {
    return iterator(this, m_bounces);
}

void Integrator::next() {
//    VolumeGridDensity *vr = m_volumeGridDensity;
//    double ds = 0.12  / (vr->Density(m_ray.origin()) / 40.0);
    double ds = 0.4 * -log(m_rng->RandomFloat());
//        double ds = 0.01;
    double g = 0.98;
//        double g = 1.0;

    double cosTheta = Distribution::heyneyGreenstein(g, *m_rng);
    double sinTheta = sqrt(1 - cosTheta*cosTheta);
    double phi = 2.0 * M_PI * m_rng->RandomFloat();

    Vector3D perpendicular = m_ray.direction().perpendicular();
    Transform phiRotation = Rotate(phi, m_ray.direction());
    perpendicular = phiRotation(perpendicular);

    Transform directionRotation = Rotatec(cosTheta, sinTheta, perpendicular);

    Vector3D direction = directionRotation(m_ray.direction());
    direction = direction.normalized();

    Point3D origin = m_ray.origin() + direction * ds;
    m_ray = Ray(origin, direction);
}

Integrator::iterator::iterator(Integrator *parent, int bounce)
    : m_parent(parent)
    , m_bounce(bounce)
{
}

Integrator::iterator::iterator(Integrator *parent)
    : m_parent(parent)
{
}

bool Integrator::iterator::operator==(const Integrator::iterator &other) const {
    return other.m_bounce == m_bounce;
}

bool Integrator::iterator::operator!=(const Integrator::iterator &other) const {
    return !(other == *this);
}

Integrator::iterator &Integrator::iterator::operator++() {
    m_bounce++;
    m_parent->next();
    return *this;
}

Ray &Integrator::iterator::operator*() {
    return m_parent->m_ray;
}
