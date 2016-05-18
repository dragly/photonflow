#ifndef VOXELIZER_H
#define VOXELIZER_H

#include <armadillo>
#include <vector>
#include "../core/transform.h"
#include "../geometry/bbox.h"
#include "../geometry/cylinderfrustum.h"

namespace photonflow {

arma::cube voxelize(std::vector<CylinderFrustum> cylinders, const photonflow::Transform &transform,
                    const BoundingBox &boundingBox, int maxExtent);

} // namespace

#endif // VOXELIZER_H
