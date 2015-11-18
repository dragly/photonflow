#include <catch.hpp>

#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"

#include <iostream>
#include <fstream>

using namespace std;

void vectors() {
    ofstream out("out.dat");

    RNG rng;
    rng.seed(1239);

    double dt = 0.1;

    Ray ray(Point(0, 0, 0), Vector(0.0, 0.0, 1.0), 0.0);
    for(int i = 0; i < 2000; i++) {
        double g = 0.99;
        double theta = acos(Distribution::heyneyGreenstein(g, rng));
        double phi = 2.0 * M_PI * rng.RandomFloat();

        Vector perpendicular = ray.d.perpendicular();
        Transform perpendicularRotation = Rotate(phi, ray.d);
        perpendicular = perpendicularRotation(perpendicular);

        Transform directionRotation = Rotate(theta, perpendicular);
        ray.d = directionRotation(ray.d);
        ray.d = ray.d.normalized();

        ray.o += ray.d * dt;
        out << ray.o.x << " " << ray.o.y << " " << ray.o.z << endl;
    }
}

TEST_CASE( "Heyney Greenstein", "[heyneygreenstein]" ) {
    SECTION("Vectors") {
        vectors();

        Vector left(1.0, 0.0, 0.0);
        Vector up(0.0, 0.0, 1.0);
        Transform rotation = Rotate(90, up);
        cout << rotation(left) << endl;
    }
}
