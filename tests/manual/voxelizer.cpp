#include <catch.hpp>

#include "schema/NeuroML_v2beta3.h"
#include "core/geometry.h"
#include "core/heyneygreenstein.h"
#include "core/transform.h"
#include "core/integrator.h"
#include <QXmlStreamReader>

#include <iostream>
#include <fstream>
#include <unordered_map>

#define ARMA_USE_HDF5

#include <QElapsedTimer>
#include <armadillo>

using namespace std;
using namespace boost::units::photonflow;
using namespace neurona;

class Cylinder
{
public:
    Cylinder(Vector3D in_start, Vector3D in_end, length in_radius)
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

    Vector3D start;
    Vector3D end;
    Vector3D center;
    length radius;
    area radius2;
    Vector3D direction;
    length h;
    length height;

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

BBox makeUnion(const BBox &b, const Vector3D &p) {
    BBox ret = b;
    ret.pMin.x = min(b.pMin.x, p.x);
    ret.pMin.y = min(b.pMin.y, p.y);
    ret.pMin.z = min(b.pMin.z, p.z);
    ret.pMax.x = max(b.pMax.x, p.x);
    ret.pMax.y = max(b.pMax.y, p.y);
    ret.pMax.z = max(b.pMax.z, p.z);
    return ret;
}

void voxelize() {
    auto path = "/home/svenni/Dropbox/projects/programming/neuroscience/neurona/neurona/hay_et_al_2011.nml";

    vector<Cylinder> cylinders;

    unordered_map<uint, schema::Segment> segments;

    BBox bbox;

    try
    {
        auto doc = schema::neuroml(path);

        for(schema::Cell &cell : doc->cells()) {
            if(cell.morphology()) {
                for(schema::Segment &segment : cell.morphology()->segments()) {
                    segments.insert({segment.id(), segment});
                }
                for(schema::Segment &segment : cell.morphology()->segments()) {
                    auto parentOptional = segment.parent();
                    auto proximalOptional = segment.proximal();
                    auto distal = segment.distal();
                    auto radius = segment.distal().diameter() * 0.5_um;
                    Vector3D start(distal.x() * 1.0_um, distal.y() * 1.0_um, distal.z() * 1.0_um);
                    bbox = makeUnion(bbox, start);
                    Vector3D end;
                    bool hasProximal = false;
                    if(proximalOptional) {
                        hasProximal = true;
                        end = Vector3D(proximalOptional->x() * 1.0_um, proximalOptional->y() * 1.0_um, proximalOptional->z() * 1.0_um);
                    } else if(parentOptional) {
                        hasProximal = true;
                        auto parentInMap = segments.find(segment.parent()->segment());
                        if(parentInMap != segments.end()) {
                            auto parentDistal = parentInMap->second.distal();
                            end = Vector3D(parentDistal.x() * 1.0_um, parentDistal.y() * 1.0_um, parentDistal.z() * 1.0_um);
                        }
                    }
                    if(hasProximal) {
                        Cylinder cylinder(start, end, radius);
                        cylinders.push_back(cylinder);
                    }
                }
            }
        }

    }
    catch (const xml_schema::Exception& e)
    {
        cerr << e << endl;
    }

    //    bbox = BBox(Point3D(-5.0_um, -3.0_um, -3.0_um), Point3D(3.0_um, 4.0_um, 4.0_um));
    //    cylinders.clear();
    //    cylinders.push_back(Cylinder(Vector3D(-0.0_um, -2.0_um, 0.0_um), Vector3D(0.0_um, 1.0_um, 0.0_um), 1.0_um));

    int N = 2048;
    length xSide = bbox.pMax[0] - bbox.pMin[0];
    length ySide = bbox.pMax[1] - bbox.pMin[1];
    length zSide = bbox.pMax[2] - bbox.pMin[2];
    length maxLen = max(xSide, max(ySide, zSide));
    double xRatio = xSide / maxLen;
    double yRatio = ySide / maxLen;
    double zRatio = zSide / maxLen;

    cout << "Ratios: " << xRatio << " " << yRatio << " " << zRatio << endl;

    // x = col, y = row, z = slice
    arma::cube voxels = arma::zeros(N * yRatio + 1, N * xRatio + 1, N * zRatio + 1);

    cout << "Min: " << bbox.pMin.x.value() << " " << bbox.pMin.y.value() << " " << bbox.pMin.z.value() << endl;
    cout << "Max: " << bbox.pMax.x.value() << " " << bbox.pMax.y.value() << " " << bbox.pMax.z.value() << endl;

    Vector3D offset(bbox.pMin);
    for(Cylinder& cylinder : cylinders) {
        cylinder = Cylinder(cylinder.start - offset,
                            cylinder.end - offset,
                            cylinder.radius);
    }
    bbox.pMin -= offset;
    bbox.pMax -= offset;

    length step = maxLen / double(N - 1);
    length eps = step / 2.0;
    cout << "Step: " << step.value() << endl;
    QElapsedTimer timer;
    timer.start();
    //#pragma omp parallel num_threads(8)
    //#pragma omp for
    //    for(int i = 0; i < int(voxels.n_slices); i++) {
    //        cout << "Iteration " << i << endl;
    //        for(int j = 0; j < int(voxels.n_rows); j++) {
    //            for(int k = 0; k < int(voxels.n_cols); k++) {

    //            }
    //        }
    //    }

    //    Vector3D p(step*double(k), step*double(j), step*double(i));
    //    p += Vector3D(bbox.pMin);
    for(Cylinder& cylinder : cylinders) {
        BBox localBounds(Point3D(-cylinder.h, -cylinder.radius, -cylinder.radius),
                         Point3D(cylinder.h, cylinder.radius, cylinder.radius));

        auto perpendicular2 = cross(Vector3D(1.0_um, 0.0_um, 0.0_um), cylinder.direction);
        Vector3D perpendicular = perpendicular2 / 1.0_um;
        Transform rotation;
        if(perpendicular.length() > 0.0_um) {
            double sinAngle = perpendicular.length().value();
            double cosAngle = sqrt(1 - sinAngle*sinAngle);

            rotation = rotatec(cosAngle, sinAngle, perpendicular);
        }
        Transform translation = translate(cylinder.center);
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
                    Vector3D p(step*double(i) + step / 2.0, step*double(j) + step / 2.0, step*double(k) + step / 2.0);
                    Vector3D diff = p - cylinder.center;
                    length distance = diff.length();
                    if(distance > cylinder.h + eps && distance > cylinder.radius + eps) {
                        continue;
                    }
                    auto yComponent = fabs(dot(diff, cylinder.direction)) / 1.0_um;
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
