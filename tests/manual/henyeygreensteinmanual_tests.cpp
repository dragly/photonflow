#include <catch.hpp>

#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace boost;
using namespace boost::accumulators;
using namespace std;
using namespace photonflow;

using histogram_type = iterator_range<std::vector<std::pair<double, double> >::iterator >;

double writeFile(string filename, histogram_type hist, int binCount)
{
    ofstream out(filename);
    double total = 0.0;
    for (unsigned int i = 1; i < hist.size() - 1; i++ )
    {
        double theta = (hist[i].first + hist[i+1].first) * 0.5;
        out << theta << " " << hist[i].second * binCount / M_PI << endl;
        total += hist[i].second;
    }
    return total;
}

SCENARIO( "Heyney Greenstein Distribution varies with different g", "[heyneygreensteindistribution_manual]" ) {
    GIVEN("A random number generator and an accumulator") {
        RNG rng;
        rng.seed(1239);

        int sampleCount = 2e6;
        int binCount = 50;

        accumulator_set<double, features<tag::density>> acc(tag::density::num_bins = binCount,
                                                             tag::density::cache_size = sampleCount);

        WHEN("we have g = 0.0") {
            for(int i = 0; i < sampleCount; i++) {
                double g = 0.0;
                double cosTheta = Distribution::heyneyGreenstein(g, rng);
                double theta = acos(cosTheta);
                acc(theta);
            }
            THEN("the total should be 1") {
                histogram_type hist = density(acc);
                double total = writeFile("distribution_0_0.dat", hist, binCount);
                REQUIRE(total == Approx(1.0));
            }
            THEN("the distribution should be uniform (for inner points)") {
                auto hist = density(acc);
                for(unsigned int i = 5; i < hist.size() - 5; i++) {
                    double theta = (hist[i].first + hist[i+1].first) * 0.5;
                    double probability = 2 * hist[i].second * binCount / M_PI / sin(theta);
                    REQUIRE(probability == Approx(1.0).epsilon(0.01));
                }
            }
        }
        WHEN("we have g = 0.98") {
            ofstream out("distribution_0_98.dat");
            for(int i = 0; i < sampleCount; i++) {
                double g = 0.98;
                double cosTheta = Distribution::heyneyGreenstein(g, rng);
                double theta = acos(cosTheta);
                acc(theta);
//                cout << theta << " " << cosTheta << endl;
            }
            THEN("the total should be 1") {
                auto hist = density(acc);

                double total = 0.0;

                for (unsigned int i = 1; i < hist.size() - 1; i++ )
                {
                    double theta = (hist[i].first + hist[i+1].first) * 0.5;
                    out << theta << " " << hist[i].second * binCount / M_PI << endl;
                    total += hist[i].second;
                }

                REQUIRE(total == Approx(1.0));
            }
        }
        WHEN("we have g = 0.5") {
            ofstream out("distribution_0_5.dat");
            for(int i = 0; i < sampleCount; i++) {
                double g = 0.5;
                double cosTheta = Distribution::heyneyGreenstein(g, rng);
                double theta = acos(cosTheta);
                acc(theta);
            }
            THEN("the total should be 1") {
                auto hist = density(acc);

                double total = 0.0;

                for (unsigned int i = 1; i < hist.size() - 1; i++ )
                {
                    double theta = (hist[i].first + hist[i+1].first) * 0.5;
                    out << theta << " " << hist[i].second * binCount / M_PI << endl;
                    total += hist[i].second;
                }

                REQUIRE(total == Approx(1.0));
            }
        }
        WHEN("we have g = 0.98") {
            ofstream out("distribution_0_98.dat");
            for(int i = 0; i < sampleCount; i++) {
                double g = 0.98;
                double cosTheta = Distribution::heyneyGreenstein(g, rng);
                double theta = acos(cosTheta);
                acc(theta);
//                cout << theta << " " << cosTheta << endl;
            }
            THEN("the total should be 1") {
                auto hist = density(acc);

                double total = 0.0;

                for (unsigned int i = 1; i < hist.size() - 1; i++ )
                {
                    double theta = (hist[i].first + hist[i+1].first) * 0.5;
                    out << theta << " " << hist[i].second * binCount / M_PI << endl;
                    total += hist[i].second;
                }

                REQUIRE(total == Approx(1.0));
            }
        }
    }
}



