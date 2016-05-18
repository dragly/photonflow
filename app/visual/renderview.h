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
    std::shared_ptr<ImageFilm> m_film;
    int m_totalSampleCount = 0;
    VolumeGridDensity m_volumeRegion;
    vector<RNG> m_randomNumberGenerators;
};

class PhotonflowSimulator : public Simulator
{
    Q_OBJECT
    Q_PROPERTY(QImage image READ image WRITE setImage NOTIFY imageChanged)
    Q_PROPERTY(double emissionFactor READ emissionFactor WRITE setEmissionFactor NOTIFY emissionFactorChanged)
    Q_PROPERTY(double scatteringCoefficient READ scatteringCoefficient WRITE setScatteringCoefficient NOTIFY scatteringCoefficientChanged)
    Q_PROPERTY(double absorptionCoefficient READ absorptionCoefficient WRITE setAbsorptionCoefficient NOTIFY absorptionCoefficientChanged)
    Q_PROPERTY(double henyeyGreensteinFactor READ henyeyGreensteinFactor WRITE setHenyeyGreensteinFactor NOTIFY henyeyGreensteinFactorChanged)
public:
    PhotonflowSimulator(QNode *parent = 0);

    QImage image();
    double emissionFactor() const;
    double absorptionCoefficient() const;
    double scatteringCoefficient() const;

    double henyeyGreensteinFactor() const;

signals:
    void imageChanged(QImage image);
    void emissionFactorChanged(double emissionFactor);
    void absorptionCoefficientChanged(double absorptionCoefficient);
    void scatteringCoefficientChanged(double scatteringCoefficient);

    void henyeyGreensteinFactorChanged(double henyeyGreensteinFactor);

public slots:
    void clear();
    void setImage(QImage image);
    void setEmissionFactor(double emissionFactor);
    void setAbsorptionCoefficient(double absorptionCoefficient);
    void setScatteringCoefficient(double scatteringCoefficient);

    void setHenyeyGreensteinFactor(double henyeyGreensteinFactor);

    void voxelize(const QVariantList &neuronSimulators);

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

    friend class PhotonflowWorker;
};

}

#endif // RENDERVIEW_H
