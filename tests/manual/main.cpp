#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "geometry/rectangle.h"

TEST_CASE("test", "[banana]") {
    SECTION("LOL") {
        Rectangle a(1.0, 2.0, 3.0, 4.0);

    }
}

