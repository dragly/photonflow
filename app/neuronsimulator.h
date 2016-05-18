#ifndef NEURONSIMULATOR_H
#define NEURONSIMULATOR_H

#include <QObject>
#include <SimVis/Simulator>
#include <geometry/cylinderfrustum.h>
#include <geometry/bbox.h>

class NeuroMlWorker : public SimulatorWorker
{
protected:
    virtual void work() override;
private:
    virtual void synchronizeSimulator(Simulator *simulator) override;
    bool m_loaded = false;
    bool m_dirty = true;
    std::vector<photonflow::CylinderFrustum> m_cylinders;
    photonflow::BoundingBox m_boundingBox;
};

class NeuronSimulator : public Simulator
{
    Q_OBJECT
    Q_PROPERTY(CylinderData* cylinderData READ cylinderData CONSTANT)

public:
    explicit NeuronSimulator(QNode *parent = 0);

    class CylinderData* cylinderData() const;
    const std::vector<photonflow::CylinderFrustum>& cylinders() const;
    const photonflow::BoundingBox &boundingBox() const;
    double scale() const;

signals:

public slots:

protected:
    virtual SimulatorWorker *createWorker() override;

private:
    CylinderData* m_cylinderData = nullptr;
    std::vector<photonflow::CylinderFrustum> m_cylinders;
    photonflow::BoundingBox m_boundingBox;

    friend class NeuroMlWorker;
};

#endif // NEURONSIMULATOR_H
