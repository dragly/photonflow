#include "neuronsimulator.h"
#include <io/neuromlreader.h>
#include <SimVis/CylinderData>

using namespace photonflow;

NeuronSimulator::NeuronSimulator(QNode *parent)
    : Simulator(parent)
    , m_cylinderData(new CylinderData(this))
{

}

CylinderData *NeuronSimulator::cylinderData() const
{
    return m_cylinderData;
}

const std::vector<CylinderFrustum> &NeuronSimulator::cylinders() const
{
    return m_cylinders;
}

const BoundingBox &NeuronSimulator::boundingBox() const
{
    return m_boundingBox;
}

SimulatorWorker *NeuronSimulator::createWorker() {
    return new NeuroMlWorker();
}

void NeuroMlWorker::work() {
    if(m_loaded) {
        return;
    }
    std::string path("/home/svenni/Dropbox/projects/programming/neuroscience/neurona/neurona/hay_et_al_2011.nml");
    NeuroMlReader reader(path);
    m_cylinders = reader.cylinders();
    m_boundingBox = reader.boundingBox();
    m_dirty = true;
    m_loaded = true;
}

void NeuroMlWorker::synchronizeSimulator(Simulator *simulator) {
    NeuronSimulator* neuronSimulator = qobject_cast<NeuronSimulator*>(simulator);
    QVector<CylinderVBOData> renderableData;
    renderableData.resize(m_cylinders.size());

    double scale = 10e-3;
    for(int i = 0; i < m_cylinders.size(); i++) {
        const CylinderFrustum& cylinder = m_cylinders.at(i);
        CylinderVBOData &renderableCylinder = renderableData[i];
        renderableCylinder.radius1 = scale * cylinder.startRadius.value();
        renderableCylinder.radius2 = scale * cylinder.endRadius.value();
        renderableCylinder.vertex1 = scale * QVector3D(cylinder.start.x.value(), cylinder.start.y.value(), cylinder.start.z.value());
        renderableCylinder.vertex2 = scale * QVector3D(cylinder.end.x.value(), cylinder.end.y.value(), cylinder.end.z.value());
    }
    neuronSimulator->cylinderData()->setData(renderableData);
    neuronSimulator->m_cylinders = m_cylinders;
    neuronSimulator->m_boundingBox = m_boundingBox;
}
