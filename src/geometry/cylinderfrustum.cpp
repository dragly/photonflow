#include "cylinderfrustum.h"

namespace photonflow {

std::ostream& operator << (std::ostream &out, const CylinderFrustum &cylinder)
{
    out << "Cylinder(" << endl;
    out << "    start: " << cylinder.start.x.value() << " " << cylinder.start.y.value() << " " << cylinder.start.z.value() << " " << endl;
    out << "    end: " << cylinder.end.x.value() << " " << cylinder.end.y.value() << " " << cylinder.end.z.value() << " " << endl;
    out << "    radius: " << cylinder.startRadius.value() << endl;
    out << ");";
    return out;
}

} // namespace
