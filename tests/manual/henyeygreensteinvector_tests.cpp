#include <catch.hpp>

#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"
#include "core/integrator.h"

#include <iostream>
#include <fstream>

using namespace std;

//void vectors() {
//    ofstream out("out.dat");

//    RNG rng;
//    rng.seed(1239);

//    Integrator integrator;

//    Ray startRay(Point(0, 0, 0), Vector(0.0, 0.0, 1.0), 0.0);

//    integrator.integrate(startRay, 2000, rng, [&out](Ray &ray) {
//        out << 0.0 << " " << 0.0 << " " << 0.0 << endl;
//    });
//}

TEST_CASE( "Heyney Greenstein", "[heyneygreenstein]" ) {
    SECTION("Vectors") {
//        vectors();

        Vector left(1.0, 0.0, 0.0);
        Vector up(0.0, 0.0, 1.0);
        Transform rotation = Rotate(90, up);
        cout << rotation(left) << endl;
    }
}
