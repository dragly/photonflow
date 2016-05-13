#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include "armadillo_includer.h"
#include "core/randomnumbergenerator.h"
#include "film/image.h"
#include "volumes/volumegrid.h"

#include <QQuickItem>
#include <QQuickPaintedItem>
#include <QImage>
#include <QMutex>

#include <SimVis/Simulator>

#include <memory>

namespace photonflow {

class PhotonflowWorker : public SimulatorWorker {

public:
    PhotonflowWorker();

protected:
    virtual void work() override;

private:
    virtual void synchronizeSimulator(Simulator *simulator) override;

    QImage m_image;
    bool m_isDataLoaded = false;
    int times = 1;
    std::shared_ptr<ImageFilm> film;
    int totalSampleCount = 0;
    VolumeGridDensity vr;
    vector<RNG> rngs;
};

class PhotonflowSimulator : public Simulator
{
    Q_OBJECT
    Q_PROPERTY(QImage image READ image WRITE setImage NOTIFY imageChanged)
public:
    PhotonflowSimulator(QNode *parent = 0);

    QImage image();

signals:

    void imageChanged(QImage image);

public slots:
    void requestIntegrate();

    void setImage(QImage image);

protected:
    virtual SimulatorWorker *createWorker() override;

private:
    void integrate();

    QImage m_image;
    bool m_isDataLoaded = false;
    int times = 1;
    std::shared_ptr<ImageFilm> film;
    int totalSampleCount = 0;
    VolumeGridDensity vr;
    vector<RNG> rngs;

    friend class PhotonflowWorker;
    QT3D_CLONEABLE(PhotonflowSimulator)
};

}

#endif // RENDERVIEW_H
