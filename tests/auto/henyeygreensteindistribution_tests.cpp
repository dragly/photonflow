#include <catch.hpp>

#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <iostream>
#include <fstream>

using namespace boost;
using namespace boost::accumulators;
using namespace std;

void distribution() {
    RNG rng;
    rng.seed(1239);

    int sampleCount = 1e6;
    int binCount = 20;

    accumulator_set<double, features<tag::density> > acc(tag::density::num_bins = binCount,
                                                         tag::density::cache_size = sampleCount);

    for(int i = 0; i < sampleCount; i++) {
        double g = 0.0;
        double cosTheta = Distribution::heyneyGreenstein(g, rng);
        acc(cosTheta);
    }

    auto hist = density(acc);

    double total = 0.0;

    for (unsigned int i = 1; i < hist.size() - 1; i++ )
    {
        REQUIRE(hist[i].second == Approx(1.0 / binCount).epsilon(0.01));
        total += hist[i].second;
    }

    REQUIRE(total == Approx(1.0));
}

TEST_CASE( "Heyney Greenstein Distribution", "[heyneygreensteindistribution]" ) {
    SECTION("Distribution") {
        distribution();
    }
}


