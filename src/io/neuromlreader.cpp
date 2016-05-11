#include "neuromlreader.h"

#include "../geometry/bbox.h"

#include <unordered_map>
#include <QXmlStreamReader>
#include <QVector3D>
#include <QFile>

namespace photonflow {

class Segment {
public:
    int id = -1;
    int parentID = -1;
    Point3D proximal;
    Point3D distal;
    Length proximalWidth = 0.0_um;
    Length distalWidth = 0.0_um;
    bool hasParentID = false;
    bool hasProximal = false;
    bool hasDistal = false;
};

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

    QString pathString = QString::fromStdString(path);
    QFile file(pathString);
    file.open(QFile::ReadOnly);
    if(!file.isOpen()) {
        cerr << "NeuroMlReader: Could not read file " << path << endl;
        return false;
    }
    QXmlStreamReader reader(&file);
    reader.readNext();
    //Reading from the file
    QVector<Segment> segments;
    bool dummySegment = true;
    Segment segment;
    while (!reader.isEndDocument())
    {
        if(reader.isStartElement()) {
            if(reader.name() == "segment") {
                if(!dummySegment) {
                    if(!segment.hasProximal && !segment.hasParentID) {
                        cerr << "WARNING: NeuroMlReader: Got segment without proximal or parent." << endl;
                    } else if(!segment.hasDistal) {
                        cerr << "WARNING: NeuroMlReader: Got segment without distal." << endl;
                    } else {
                        segments.append(segment);
                    }
                }
                dummySegment = false;
                segment = Segment();
                segment.id = reader.attributes().value("id").toInt();
                if(reader.attributes().hasAttribute("parent")) {
                    segment.parentID = reader.attributes().value("parent").toInt();
                    segment.hasParentID = true;
                }
            }
            if(reader.name() == "parent") {
                segment.parentID = reader.attributes().value("segment").toInt();
                segment.hasParentID = true;
            }
            if(reader.name() == "proximal") {
                //                qDebug() << reader.attributes().value("x") << reader.attributes().value("y") << reader.attributes().value("z");
                segment.proximal = Point3D(reader.attributes().value("x").toDouble() * 1.0_um,
                                           reader.attributes().value("y").toDouble() * 1.0_um,
                                           reader.attributes().value("z").toDouble() * 1.0_um);
                segment.proximalWidth = reader.attributes().value("diameter").toDouble() * 1.0_um;
                segment.hasProximal = true;
            }
            if(reader.name() == "distal") {
                //                qDebug() << reader.attributes().value("x") << reader.attributes().value("y") << reader.attributes().value("z") << reader.attributes().value("diameter");
                segment.distal = Point3D(reader.attributes().value("x").toDouble() * 1.0_um,
                                         reader.attributes().value("y").toDouble() * 1.0_um,
                                         reader.attributes().value("z").toDouble() * 1.0_um);
                segment.distalWidth = reader.attributes().value("diameter").toDouble() * 1.0_um;
                segment.hasDistal = true;
            }
        }
        reader.readNext();
    }
    segments.append(segment);

    for(Segment &segment : segments) {
        if(!segment.hasProximal && segment.hasParentID) {
            for(Segment &otherSegment : segments) {
                if(otherSegment.id == segment.parentID) {
                    Segment &parent = otherSegment;
                    segment.proximal = parent.distal;
                    break;
                }
            }
        }
        if(segment.proximalWidth == 0.0_um) {
            segment.proximalWidth = segment.distalWidth;
        }
    }

    m_cylinders.clear();
    for(const Segment &segment : segments) {
        CylinderFrustum cylinder(segment.proximal,
                                 segment.distal,
                                 segment.proximalWidth * 0.5,
                                 segment.distalWidth * 0.5);

        boundingBox = makeUnion(boundingBox, cylinder.start);
        boundingBox = makeUnion(boundingBox, cylinder.end);
        m_cylinders.push_back(cylinder);
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
