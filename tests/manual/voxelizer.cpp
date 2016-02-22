#include <catch.hpp>

#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"
#include "core/integrator.h"
#include "geometry/cylinderfrustum.h"
#include "io/neuromlreader.h"
#include "io/voxelizer.h"

#include <iostream>
#include <fstream>
#include <unordered_map>

#define ARMA_USE_HDF5

#include <QElapsedTimer>
#include <armadillo>

using namespace std;
using namespace photonflow;

// TODO add tests to check that voxelization works properly

void voxelizeTest() {
    std::string path("/home/svenni/Dropbox/projects/programming/neuroscience/neurona/neurona/hay_et_al_2011.nml");
    NeuroMlReader reader(path);
    vector<CylinderFrustum> cylinders = reader.cylinders();
    BoundingBox boundingBox = reader.boundingBox();

    arma::cube voxels = voxelize(cylinders, boundingBox, 2048);

    cout << "Voxels minmax: " << voxels.min() << " " << voxels.max() << endl;
    cout << "Voxels shape: " << voxels.n_slices << " " << voxels.n_rows << " " << voxels.n_cols << endl;
    cout << "Saving data to disk..." << endl;
    voxels.save("/tmp/out.raw", arma::raw_binary);
}

TEST_CASE( "Voxelizer", "[voxelizer]" ) {
    SECTION("Voxelizers") {
        voxelizeTest();
    }
}
