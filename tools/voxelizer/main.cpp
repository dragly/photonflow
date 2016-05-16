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

#include <QElapsedTimer>
#include <armadillo>

using namespace std;
using namespace photonflow;

// TODO add tests to check that voxelization works properly

int main() {
    std::string path("/home/svenni/Dropbox/projects/programming/neuroscience/neurona/neurona/hay_et_al_2011.nml");
    NeuroMlReader reader(path);
    vector<CylinderFrustum> cylinders = reader.cylinders();
    BoundingBox boundingBox = reader.boundingBox();

    arma::cube voxels = voxelize(cylinders, boundingBox, 1024);

    cout << "Voxels minmax: " << voxels.min() << " " << voxels.max() << endl;
    cout << "Voxels shape: " << voxels.n_slices << " " << voxels.n_rows << " " << voxels.n_cols << endl;
    cout << "Saving data to disk..." << endl;
//    voxels.save("/tmp/out.raw", arma::raw_binary);
    voxels.save("/tmp/out.h5", arma::hdf5_binary);
}
