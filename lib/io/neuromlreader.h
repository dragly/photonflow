#ifndef NEUROMLREADER_H
#define NEUROMLREADER_H

#include <string>
#include <vector>
#include "../geometry/bbox.h"
#include "../geometry/cylinderfrustum.h"

namespace photonflow {

class NeuroMlReader
{
public:
    NeuroMlReader(std::string path = std::string());

    std::string path() const;
    bool load(const std::string &path);

    std::vector<CylinderFrustum> cylinders() const;
    BBox boundingBox() const;

private:
    std::string m_path;
    std::vector<CylinderFrustum> m_cylinders;
    BBox m_boundingBox;
};

} // namespace

#endif // NEUROMLREADER_H
