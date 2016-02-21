#include <catch.hpp>

//#include "schema/NeuroML_v2beta3.h"
#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"
#include "core/integrator.h"
#include <pugixml.hpp>
//#include <QXmlStreamReader>

#include <iostream>
#include <fstream>
#include <unordered_map>

#define ARMA_USE_HDF5

#include <QElapsedTimer>
#include <armadillo>

using namespace std;
using namespace boost::units::photonflow;
//using namespace neurona;
using namespace pugi;

// TODO use frustum instead of cylinder (two radii)
// TODO add tests to check that voxelization works properly

class Cylinder
{
public:
    Cylinder(Point3D in_start, Point3D in_end, Length in_radius)
        : start(in_start)
        , end(in_end)
        , center((in_start + in_end) * 0.5)
        , radius(in_radius)
        , radius2(in_radius*in_radius)
        , direction((in_end - in_start).normalized())
    {
        height = (in_end - in_start).length();
        h = (in_end - in_start).length() * 0.5;
    }

    Point3D start;
    Point3D end;
    Point3D center;
    Length radius;
    Area radius2;
    Vector3D direction;
    Length h;
    Length height;

    friend std::ostream& operator<< (std::ostream &out, const Point3D &point);
};

std::ostream& operator << (std::ostream &out, const Cylinder &cylinder)
{
    out << "Cylinder(" << endl;
    out << "    start: " << cylinder.start.x.value() << " " << cylinder.start.y.value() << " " << cylinder.start.z.value() << " " << endl;
    out << "    end: " << cylinder.end.x.value() << " " << cylinder.end.y.value() << " " << cylinder.end.z.value() << " " << endl;
    out << "    radius: " << cylinder.radius.value() << endl;
    out << ");";
    return out;
}

void voxelize() {
    auto path = "/home/svenni/Dropbox/projects/programming/neuroscience/neurona/neurona/hay_et_al_2011.nml";

    vector<Cylinder> cylinders;

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
        xml_node distal = segment.child("distal");
        Length radius = distal.attribute("diameter").as_double() * 0.5_um;
        xml_node proximal = segment.child("proximal");
        xml_node parent = segment.child("parent");
        if(!distal) {
            continue;
        }
        Point3D start(distal.attribute("x").as_double() * 1.0_um,
                      distal.attribute("y").as_double() * 1.0_um,
                      distal.attribute("z").as_double() * 1.0_um);
        Point3D end;
        bbox = makeUnion(bbox, start);
        bool hasProximal = false;
        if(proximal) {
            hasProximal = true;
            end = Point3D(proximal.attribute("x").as_double() * 1.0_um,
                          proximal.attribute("y").as_double() * 1.0_um,
                          proximal.attribute("z").as_double() * 1.0_um);
        } else if(parent) {
            auto parentInMap = segments.find(parent.attribute("segment").as_uint());
            if(parentInMap != segments.end()) {
                hasProximal = true;
                xml_node parentDistal = parentInMap->second.child("distal");
                end = Point3D(parentDistal.attribute("x").as_double() * 1.0_um,
                              parentDistal.attribute("y").as_double() * 1.0_um,
                              parentDistal.attribute("z").as_double() * 1.0_um);
            }
        }
        if(hasProximal) {
            Cylinder cylinder(start, end, radius);
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
    for(Cylinder& cylinder : cylinders) {
        cylinder = Cylinder(cylinder.start - offset,
                            cylinder.end - offset,
                            cylinder.radius);
    }
    bbox.pMin -= offset;
    bbox.pMax -= offset;

    Length step = maxLen / double(N - 1);
    Length eps = step / 2.0;
    cout << "Step: " << step.value() << endl;
    QElapsedTimer timer;
    timer.start();
    for(Cylinder& cylinder : cylinders) {
        BBox localBounds(Point3D(-cylinder.h, -cylinder.radius, -cylinder.radius),
                         Point3D(cylinder.h, cylinder.radius, cylinder.radius));

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
                    if(distance > cylinder.h + eps && distance > cylinder.radius + eps) {
                        continue;
                    }
                    auto yComponent = fabs(dot(diff, cylinder.direction * 1.0_um)) / 1.0_um;
                    if(yComponent <= cylinder.h + eps) {
                        auto y2 = yComponent*yComponent;
                        auto diff2 = dot(diff, diff);
                        auto distanceToAxis = sqrt(diff2 - y2);
                        if(distanceToAxis <= cylinder.radius + eps) {
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
    voxels.save("/tmp/out.raw", arma::raw_binary);
}

TEST_CASE( "Voxelizer", "[voxelizer]" ) {
    SECTION("Voxelizers") {
        voxelize();
    }
}
