#ifndef CYLINDERFRUSTUM_H
#define CYLINDERFRUSTUM_H

#include <ostream>
#include "point3d.h"
#include "../core/units.h"

namespace photonflow {

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

} // namespace

#endif // CYLINDERFRUSTUM_H
