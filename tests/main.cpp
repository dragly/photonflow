#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "core/heyneygreenstein.h"

#include <iostream>
#include <fstream>

using namespace std;

TEST_CASE( "Proper Heyney Greenstein distribution", "[heyneygreenstein]" ) {
    RNG rng;
    rng.seed(1239);

    ofstream file("out.dat");

    for(int i = 0; i < 1e6; i++) {
        double g = 0.0;
        double theta = Phases::phaseHeyneyGreenstein(g, rng);
        file << theta << endl;
    }
}
