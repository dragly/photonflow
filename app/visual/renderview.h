#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include "armadillo_includer.h"
#include "core/randomnumbergenerator.h"
#include "geometry/cylinderfrustum.h"
#include "film/image.h"
#include "volumes/volumegrid.h"

#include <QQuickItem>
#include <QQuickPaintedItem>
#include <QImage>
#include <QMutex>

#include <SimVis/Simulator>

#include <memory>

class NeuronSimulator;

namespace photonflow {

struct Cell {
    std::vector<CylinderFrustum> m_cylinders;
};

class PhotonflowWorker : public SimulatorWorker {

public:
    PhotonflowWorker();

protected:
    virtual void work() override;

private:
    virtual void synchronizeSimulator(Simulator *simulator) override;

    QImage m_image;
    bool m_isDataLoaded = false;
    std::shared_ptr<ImageFilm> m_film;
    int m_completedSampleCount = 0;
    int m_iterationCount = 0;
    VolumeGridDensity m_volumeRegion;
    vector<RNG> m_randomNumberGenerators;
    std::vector<CylinderFrustum> m_cylinders;
    BoundingBox m_boundingBox;
    std::vector<Cell> m_cells;
    double m_renderTime = 0.0;
    Length m_lensRadius = 0.0_um;
    Length m_focalDepth = 0.0_um;
    double m_fieldOfView = 0.0;
};

class PhotonflowSimulator : public Simulator
{
    Q_OBJECT
    Q_PROPERTY(QImage image READ image WRITE setImage NOTIFY imageChanged)
    Q_PROPERTY(double emissionFactor READ emissionFactor WRITE setEmissionFactor NOTIFY emissionFactorChanged)
    Q_PROPERTY(double scatteringCoefficient READ scatteringCoefficient WRITE setScatteringCoefficient NOTIFY scatteringCoefficientChanged)
    Q_PROPERTY(double absorptionCoefficient READ absorptionCoefficient WRITE setAbsorptionCoefficient NOTIFY absorptionCoefficientChanged)
    Q_PROPERTY(double henyeyGreensteinFactor READ henyeyGreensteinFactor WRITE setHenyeyGreensteinFactor NOTIFY henyeyGreensteinFactorChanged)
    Q_PROPERTY(double renderTime READ renderTime NOTIFY renderTimeChanged)
    Q_PROPERTY(double fieldOfView READ fieldOfView WRITE setFieldOfView NOTIFY fieldOfViewChanged)
    Q_PROPERTY(double focalDepth READ focalDepth WRITE setFocalDepth NOTIFY focalDepthChanged)
    Q_PROPERTY(double lensRadius READ lensRadius WRITE setLensRadius NOTIFY lensRadiusChanged)
    Q_PROPERTY(int completedSampleCount READ completedSampleCount NOTIFY completedSampleCountChanged)
    Q_PROPERTY(Qt3DCore::QEntity* camera READ camera WRITE setCamera NOTIFY cameraChanged)

public:
    PhotonflowSimulator(QNode *parent = 0);

    QImage image();
    double emissionFactor() const;
    double absorptionCoefficient() const;
    double scatteringCoefficient() const;
    double henyeyGreensteinFactor() const;
    double renderTime() const;
    Qt3DCore::QEntity* camera() const;
    int completedSampleCount() const;
    double fieldOfView() const;
    double lensRadius() const;
    double focalDepth() const;

signals:
    void imageChanged(QImage image);
    void emissionFactorChanged(double emissionFactor);
    void absorptionCoefficientChanged(double absorptionCoefficient);
    void scatteringCoefficientChanged(double scatteringCoefficient);
    void henyeyGreensteinFactorChanged(double henyeyGreensteinFactor);
    void renderTimeChanged(double renderTime);
    void cameraChanged(Qt3DCore::QEntity* camera);
    void completedSampleCountChanged(int completedSampleCount);
    void fieldOfViewChanged(double fieldOfView);
    void lensRadiusChanged(double lensRadius);
    void focalDepthChanged(double focalDepth);

public slots:
    void clear();
    void setImage(QImage image);
    void setEmissionFactor(double emissionFactor);
    void setAbsorptionCoefficient(double absorptionCoefficient);
    void setScatteringCoefficient(double scatteringCoefficient);
    void setHenyeyGreensteinFactor(double henyeyGreensteinFactor);
    void voxelize(const QVariantList &neuronSimulators);
    void setCamera(Qt3DCore::QEntity* camera);
    void setFieldOfView(double fieldOfView);
    void setLensRadius(double lensRadius);
    void setFocalDepth(double focalDepth);

protected:
    virtual SimulatorWorker *createWorker() override;

private:
    QImage m_image;
    arma::cube m_data;
    bool m_isDataLoaded = false;
    bool m_dataDirty = false;
    bool m_clearRequested = false;
    double m_emissionCoefficient = 0.1;
    double m_absorptionCoefficient = 1.0;
    double m_scatteringCoefficient = 0.1;
    double m_henyeyGreensteinFactor = 1.0;
    double m_lensRadius = 0.0;
    double m_focalDepth = 50.0;
    double m_fieldOfView = 60.0 / 180.0;
    double m_renderTime = 0.0;
    std::vector<CylinderFrustum> m_cylinders;

    friend class PhotonflowWorker;
    Qt3DCore::QEntity* m_camera;
    int m_completedSampleCount = 0;
};

}

#endif // RENDERVIEW_H
