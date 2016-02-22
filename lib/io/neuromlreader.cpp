#include "neuromlreader.h"

#include "../geometry/bbox.h"

#include <pugixml.hpp>
#include <unordered_map>

namespace photonflow {

using pugi::xml_document;
using pugi::xml_parse_result;
using pugi::xml_node;

NeuroMlReader::NeuroMlReader(std::string path)
    : m_path(path)
{
    load(m_path);
}

std::string NeuroMlReader::path() const
{
    return m_path;
}

bool NeuroMlReader::load(const std::string &path)
{
    m_path = path;
    if(m_path.empty()) {
        cerr << "NeuroMlReader: Path is empty, cannot read file." << endl;
        return false;
    }

    m_cylinders.clear();

    BoundingBox boundingBox;

    xml_document doc;
    xml_parse_result result = doc.load_file(m_path.c_str());
    if(!result) {
        cerr << "Could not parse XML document " << m_path << endl;
        return false;
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
        boundingBox = makeUnion(boundingBox, distalPoint);
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
            boundingBox = makeUnion(boundingBox, proximalPoint);
            CylinderFrustum cylinder(proximalPoint, distalPoint, proximalRadius, distalRadius);
            m_cylinders.push_back(cylinder);
        }
    }
    m_boundingBox = boundingBox;
    return true;
}

std::vector<CylinderFrustum> NeuroMlReader::cylinders() const
{
    return m_cylinders;
}

BoundingBox NeuroMlReader::boundingBox() const
{
    return m_boundingBox;
}

}
