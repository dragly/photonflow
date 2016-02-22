#include "voxelizer.h"

#include "../core/transform.h"

namespace photonflow {

arma::cube voxelize(std::vector<CylinderFrustum> cylinders, const BBox &boundingBox, int maxExtent)
{
    Length xSide = boundingBox.pMax[0] - boundingBox.pMin[0];
    Length ySide = boundingBox.pMax[1] - boundingBox.pMin[1];
    Length zSide = boundingBox.pMax[2] - boundingBox.pMin[2];
    Length maxLen = max(xSide, max(ySide, zSide));
    double xRatio = xSide / maxLen;
    double yRatio = ySide / maxLen;
    double zRatio = zSide / maxLen;

    // x = col, y = row, z = slice
    arma::cube voxels = arma::zeros(maxExtent * yRatio + 1, maxExtent * xRatio + 1, maxExtent * zRatio + 1);

    Length3D offset(boundingBox.pMin);
    for(CylinderFrustum& cylinder : cylinders) {
        cylinder = CylinderFrustum(cylinder.start - offset,
                                   cylinder.end - offset,
                                   cylinder.startRadius,
                                   cylinder.endRadius);
    }
    BBox offsetBox = boundingBox;
    offsetBox.pMin -= offset;
    offsetBox.pMax -= offset;

    Length step = maxLen / double(maxExtent - 1);
    Length eps = step / 2.0;
    cout << "Step: " << step.value() << endl;
    for(CylinderFrustum& cylinder : cylinders) {
        BBox localBounds(Point3D(-cylinder.h, -cylinder.startRadius, -cylinder.startRadius),
                         Point3D(cylinder.h, cylinder.startRadius, cylinder.startRadius));

        Vector3D perpendicular = cross(Vector3D(1.0, 0.0, 0.0), cylinder.direction);
        Transform rotation;
        double sinAngle = perpendicular.length();
        if(sinAngle > 0.0) {
            if(sinAngle > 1.0) {
                sinAngle -= 2.2204460492503131e-16; // typical machine precision error
            }
            double cosAngle = sqrt(1.0 - sinAngle*sinAngle);

            rotation = rotate(cosAngle, sinAngle, perpendicular);
        }
        Transform translation = translate(Length3D(cylinder.center));
        BBox bounds = translation(rotation(localBounds));
        bounds.expand(eps);

        int istart = bounds.pMin.x / step;
        int jstart = bounds.pMin.y / step;
        int kstart = bounds.pMin.z / step;

        int iend = bounds.pMax.x / step + 1;
        int jend = bounds.pMax.y / step + 1;
        int kend = bounds.pMax.z / step + 1;

        for(int i = istart; i < iend + 1; i++) {
            for(int j = jstart; j < jend + 1; j++) {
                for(int k = kstart; k < kend + 1; k++) {
                    Point3D p(step * (double(i) + 0.5), step * (double(j) + 0.5), step * (double(k) + 0.5));
                    Length3D diff = p - cylinder.center;
                    Length distance = diff.length();
                    if(distance > cylinder.h + eps && distance > cylinder.startRadius + eps) {
                        continue;
                    }
                    auto yComponent = dot(diff, cylinder.direction * 1.0_um) / 1.0_um;
                    if(fabs(yComponent) <= cylinder.h + eps) {
                        auto y2 = yComponent*yComponent;
                        auto diff2 = dot(diff, diff);
                        auto distanceToAxis = sqrt(diff2 - y2);
                        double endProportion = (yComponent + cylinder.h) / (2.0 * cylinder.h);
                        Length radius = cylinder.startRadius * (1 - endProportion) + endProportion * cylinder.endRadius;
                        if(distanceToAxis <= radius + eps) {
                            if(voxels.in_range(j, i, k)) {
                                voxels(j, i, k) = 1.0;
                            }
                        }
                    }
                }
            }
        }
    }
    return voxels;
}

} // namespace
