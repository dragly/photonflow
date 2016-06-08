#include "bbox.h"

namespace photonflow {

std::ostream& operator<<(ostream &out, const BoundingBox &boundingBox)
{
    out << boundingBox.pMin << boundingBox.pMax;
    return out;
}

} // namespace
