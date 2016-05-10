#ifndef NEURONSIMULATOR_H
#define NEURONSIMULATOR_H

#include <QObject>
#include <SimVis/Simulator>
#include <geometry/cylinderfrustum.h>

class NeuroMlWorker : public SimulatorWorker
{
protected:
    virtual void work() override;
private:
    virtual void synchronizeSimulator(Simulator *simulator) override;
    bool m_loaded = false;
    bool m_dirty = true;
    vector<photonflow::CylinderFrustum> m_cylinders;
};

class NeuronSimulator : public Simulator
{
    Q_OBJECT
    Q_PROPERTY(CylinderData* cylinderData READ cylinderData CONSTANT)

public:
    explicit NeuronSimulator(QNode *parent = 0);

    class CylinderData* cylinderData() const;

signals:

public slots:

protected:
    virtual SimulatorWorker *createWorker() override;

private:
    CylinderData* m_cylinderData = nullptr;

    QT3D_CLONEABLE(NeuronSimulator)
};

#endif // NEURONSIMULATOR_H
