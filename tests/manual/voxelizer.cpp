#include <catch.hpp>

#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"
#include "core/integrator.h"
#include <pugixml.hpp>

#include <iostream>
#include <fstream>
#include <unordered_map>

#define ARMA_USE_HDF5

#include <QElapsedTimer>
#include <armadillo>

using namespace std;
using namespace photonflow;
using namespace pugi;

// TODO add tests to check that voxelization works properly

class CylinderFrustum
{
public:
    CylinderFrustum(Point3D in_start, Point3D in_end, Length in_startRadius, Length in_endRadius)
        : start(in_start)
        , end(in_end)
        , center((in_start + in_end) * 0.5)
        , startRadius(in_startRadius)
        , endRadius(in_endRadius)
        , direction((in_end - in_start).normalized())
    {
        height = (in_end - in_start).length();
        h = (in_end - in_start).length() * 0.5;
    }

    Point3D start;
    Point3D end;
    Point3D center;
    Length startRadius;
    Length endRadius;
    Vector3D direction;
    Length h;
    Length height;

    friend std::ostream& operator<< (std::ostream &out, const Point3D &point);
};

std::ostream& operator << (std::ostream &out, const CylinderFrustum &cylinder)
{
    out << "Cylinder(" << endl;
    out << "    start: " << cylinder.start.x.value() << " " << cylinder.start.y.value() << " " << cylinder.start.z.value() << " " << endl;
    out << "    end: " << cylinder.end.x.value() << " " << cylinder.end.y.value() << " " << cylinder.end.z.value() << " " << endl;
    out << "    radius: " << cylinder.startRadius.value() << endl;
    out << ");";
    return out;
}

void voxelize() {
    auto path = "/home/svenni/Dropbox/projects/programming/neuroscience/neurona/neurona/hay_et_al_2011.nml";

    vector<CylinderFrustum> cylinders;

    BBox bbox;

    xml_document doc;
    xml_parse_result result = doc.load_file(path);
    if(!result) {
        cerr << "Could not parse XML document " << path << endl;
        return;
    }
    xml_node cell = doc.child("neuroml").child("cell").child("morphology");
    unordered_map<uint, xml_node> segments;
    for(xml_node segment : cell.children("segment")) {
        uint id = segment.attribute("id").as_uint(0);
        segments.insert({id, segment});
    }
    for(xml_node segment : cell.children("segment")) {
        uint id = segment.attribute("id").as_uint(0);
        xml_node distalNode = segment.child("distal");
        Length distalRadius = distalNode.attribute("diameter").as_double() * 0.5_um;
        Length proximalRadius;
        xml_node proximal = segment.child("proximal");
        xml_node parent = segment.child("parent");
        if(!distalNode) {
            continue;
        }
        Point3D distalPoint(distalNode.attribute("x").as_double() * 1.0_um,
                      distalNode.attribute("y").as_double() * 1.0_um,
                      distalNode.attribute("z").as_double() * 1.0_um);
        Point3D proximalPoint;
        bbox = makeUnion(bbox, distalPoint);
        bool hasProximal = false;
        if(proximal) {
            hasProximal = true;
            proximalPoint = Point3D(proximal.attribute("x").as_double() * 1.0_um,
                          proximal.attribute("y").as_double() * 1.0_um,
                          proximal.attribute("z").as_double() * 1.0_um);
            proximalRadius = proximal.attribute("diameter").as_double() * 0.5_um;
        } else if(parent) {
            auto parentInMap = segments.find(parent.attribute("segment").as_uint());
            if(parentInMap != segments.end()) {
                hasProximal = true;
                xml_node parentDistal = parentInMap->second.child("distal");
                proximalPoint = Point3D(parentDistal.attribute("x").as_double() * 1.0_um,
                              parentDistal.attribute("y").as_double() * 1.0_um,
                              parentDistal.attribute("z").as_double() * 1.0_um);
                proximalRadius = parentDistal.attribute("diameter").as_double() * 0.5_um;
            }
        }
        if(hasProximal) {
            CylinderFrustum cylinder(proximalPoint, distalPoint, proximalRadius, distalRadius);
            cylinders.push_back(cylinder);
        }
    }

    cout << "Cylinders: " << cylinders.size() << endl;

    int N = 2048;
    Length xSide = bbox.pMax[0] - bbox.pMin[0];
    Length ySide = bbox.pMax[1] - bbox.pMin[1];
    Length zSide = bbox.pMax[2] - bbox.pMin[2];
    Length maxLen = max(xSide, max(ySide, zSide));
    double xRatio = xSide / maxLen;
    double yRatio = ySide / maxLen;
    double zRatio = zSide / maxLen;

    cout << "Ratios: " << xRatio << " " << yRatio << " " << zRatio << endl;

    // x = col, y = row, z = slice
    arma::cube voxels = arma::zeros(N * yRatio + 1, N * xRatio + 1, N * zRatio + 1);

    cout << "Min: " << bbox.pMin.x.value() << " " << bbox.pMin.y.value() << " " << bbox.pMin.z.value() << endl;
    cout << "Max: " << bbox.pMax.x.value() << " " << bbox.pMax.y.value() << " " << bbox.pMax.z.value() << endl;

    Length3D offset(bbox.pMin);
    for(CylinderFrustum& cylinder : cylinders) {
        cylinder = CylinderFrustum(cylinder.start - offset,
                                   cylinder.end - offset,
                                   cylinder.startRadius,
                                   cylinder.endRadius);
    }
    bbox.pMin -= offset;
    bbox.pMax -= offset;

    Length step = maxLen / double(N - 1);
    Length eps = step / 2.0;
    cout << "Step: " << step.value() << endl;
    QElapsedTimer timer;
    timer.start();
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

    cout << "Time: " << timer.elapsed() << endl;

    cout << "Voxels minmax: " << voxels.min() << " " << voxels.max() << endl;
    cout << "Voxels shape: " << voxels.n_slices << " " << voxels.n_rows << " " << voxels.n_cols << endl;
    cout << "Saving data to disk..." << endl;
    voxels.save("/tmp/out.raw", arma::raw_binary);
}

TEST_CASE( "Voxelizer", "[voxelizer]" ) {
    SECTION("Voxelizers") {
        voxelize();
    }
}
