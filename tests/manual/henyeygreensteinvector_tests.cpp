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

        Vector perpendicular = ray.direction().perpendicular();
        Transform phiRotation = Rotate(phi, ray.direction());
        perpendicular = phiRotation(perpendicular);

        Transform directionRotation = Rotate(theta, perpendicular);

        Vector direction = directionRotation(ray.direction());
        direction = direction.normalized();

        Point origin = ray.origin() + direction * dt;
        ray = Ray(origin, direction);
        out << ray.origin() << endl;
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
