#ifndef VOXELIZER_H
#define VOXELIZER_H

#include <armadillo>
#include <vector>
#include "../geometry/bbox.h"
#include "../geometry/cylinderfrustum.h"

namespace photonflow {

arma::cube voxelize(std::vector<CylinderFrustum> cylinders, const BBox &boundingBox, int maxExtent);

} // namespace

#endif // VOXELIZER_H
