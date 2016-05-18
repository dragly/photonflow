#include "renderview.h"

#include "armadillo_includer.h"
#include "cameras/perspective.h"
#include "core/heyneygreenstein.h"
#include "core/randomnumbergenerator.h"
#include "film/image.h"
#include "filters/box.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "samplers/random.h"
#include "volumes/volumegrid.h"
#include "neuronsimulator.h"
#include "io/voxelizer.h"

#include <Qt3DCore/QTransform>
#include <QElapsedTimer>
#include <QPainter>
#include <QMutexLocker>

#include <iostream>
#include <memory>
#include <omp.h>

using namespace std;
using namespace arma;

namespace photonflow {

static const int threadCount = 8;

PhotonflowSimulator::PhotonflowSimulator(QNode *parent)
    : Simulator(parent)
{
}

QImage PhotonflowSimulator::image()
{
    return m_image;
}

double PhotonflowSimulator::emissionFactor() const
{
    return m_emissionCoefficient;
}

double PhotonflowSimulator::absorptionCoefficient() const
{
    return m_absorptionCoefficient;
}

double PhotonflowSimulator::scatteringCoefficient() const
{
    return m_scatteringCoefficient;
}

double PhotonflowSimulator::henyeyGreensteinFactor() const
{
    return m_henyeyGreensteinFactor;
}


void PhotonflowSimulator::clear()
{
    m_clearRequested = true;
    m_image.fill(QColor("black"));
    emit imageChanged(m_image);
}

void PhotonflowSimulator::setImage(QImage image)
{
    if (m_image == image)
        return;

    m_image = image;
    emit imageChanged(image);
}

void PhotonflowSimulator::setEmissionFactor(double emissionFactor)
{
    if (m_emissionCoefficient == emissionFactor)
        return;

    m_emissionCoefficient = emissionFactor;
    emit emissionFactorChanged(emissionFactor);
}

void PhotonflowSimulator::setAbsorptionCoefficient(double absorptionCoefficient)
{
    if (m_absorptionCoefficient == absorptionCoefficient)
        return;

    m_absorptionCoefficient = absorptionCoefficient;
    emit absorptionCoefficientChanged(absorptionCoefficient);
}

void PhotonflowSimulator::setScatteringCoefficient(double scatteringCoefficient)
{
    if (m_scatteringCoefficient == scatteringCoefficient)
        return;

    m_scatteringCoefficient = scatteringCoefficient;
    emit scatteringCoefficientChanged(scatteringCoefficient);
}

void PhotonflowSimulator::setHenyeyGreensteinFactor(double henyeyGreensteinFactor)
{
    if (m_henyeyGreensteinFactor == henyeyGreensteinFactor)
        return;

    m_henyeyGreensteinFactor = henyeyGreensteinFactor;
    emit henyeyGreensteinFactorChanged(henyeyGreensteinFactor);
}

PhotonflowWorker::PhotonflowWorker()
    : SimulatorWorker()
{

    qDebug() << "PhotonflowWorker: Init";
    QElapsedTimer timer;
    timer.start();
    //    arma::Cube<short> data;
    arma::cube data;
    //    data.load("/home/svenni/Dropbox/projects/programming/neuroscience/photonflow/photonflow/notebooks/output.hdf5", hdf5_binary);
//    data.load("/home/svenni/Dropbox/projects/programming/neuroscience/photonflow/photonflow/notebooks/volume.hdf5", hdf5_binary);
    data.load("/tmp/out.h5", hdf5_binary);

    qDebug() << "PhotonflowWorker: Loaded file";
    qDebug() << "Data size:" << data.n_rows << data.n_cols << data.n_slices;
    qDebug() << "Data load time:" << timer.elapsed() << "ms";
    qDebug() << "Data max value: " << data.max();

//    data *= 255;

    BoundingBox bbox;
    photonflow::Length side = 100.0_um;
    bbox.pMin = Point3D(-side, -0.2*side, -0.2*side);
    bbox.pMax = Point3D(side, 0.2*side, 0.2*side);

    double gg = 1.0;
    double angle = 0.0;

    Transform translation = translate(Length3D(0.1 * side, 0.3 * side, 2.0 * side));
//    Transform translation = translate(Length3D(0.0 * side, 0.0 * side, 0.0 * side));
    Transform rotation = rotate(angle, Vector3D(0.0, 1.0, 0.0));

    Transform boxTransform = translation*rotation;

    Spectrum sigma_a(0.95);
    Spectrum sigma_s(0.1);
    Spectrum emita(0.1);

    qDebug() << "PhotonflowWorker: Converting file...";
//    arma::Cube<short> dataShort = arma::conv_to<arma::Cube<short>>::from(data*32000);

    qDebug() << "PhotonflowWorker: Creating volume object...";
    m_volumeRegion = VolumeGridDensity(sigma_a, sigma_s, gg, emita, bbox, boxTransform, data);

    qDebug() << "PhotonflowWorker: Init complete";
}

void PhotonflowWorker::work()
{
    if(m_randomNumberGenerators.size() < threadCount) {
        m_randomNumberGenerators.resize(threadCount);
        int i = 0;
        for(RNG& rng : m_randomNumberGenerators) {
            rng.seed(1325125 ^ i);
            i++;
        }
    }
    QElapsedTimer timer;
    timer.start();

//    const QSize size = boundingRect().size().toSize();
    QSize size(640, 480); // TODO actually get size from simulator
    if(size.width() <= 0 || size.height() <= 0 || size.width() > 1e6 || size.height() > 1e6) {
        qWarning() << "WARNING: Integrate returns due to invalid size:" << size;
        return;
    }

    const int requestedSampleCount = 1;
    const int bounces = 1000;

    const int width = size.width();
    const int height = size.height();

    const Transform cameraTransform = translate(Length3D(0.0_um, 0.0_um, 0.0_um));
    Rectangle screenWindow(-width / 2.0, -height / 2.0, width, height);
    const double crop[4] = {0.0, 1.0, 0.0, 1.0};

    double boxSize = 1.0;
    BoxFilter filter(boxSize * 0.5, boxSize * 0.5);

    if(m_image.size() != size || !m_film) {
        qDebug() << "Creating image of size" << size;
        m_image = QImage(size, QImage::Format_ARGB32);
        for(int x = 0; x < width; x++) {
            for(int y = 0; y < height; y++) {
                m_image.setPixel(x, y, QColor(0.0, 0.0, 0.0, 255.0).rgba());
            }
        }
        m_film = make_shared<ImageFilm>(width, height, &filter, crop);
    }
    const auto sopen = 0.0_us;
    const auto sclose = 1.0_us;
    const auto lensr = 0.0_um;
    const auto focald = 3.5_um;
    const double fov = 0.2;
    const PerspectiveCamera camera(cameraTransform, screenWindow, sopen, sclose, lensr, focald, fov, m_film);

#pragma omp parallel num_threads(threadCount) // OpenMP
    {
        RNG& rng = m_randomNumberGenerators.at(omp_get_thread_num());
        int actualCount = 0;
        RandomSampler sampler(0, width, 0, height, requestedSampleCount, 0.0_us, 1.0_us);
        int maxSampleCount = sampler.maximumSampleCount();

        Sample originalSample;
        originalSample.Add1D(1);
        Sample* samples = originalSample.Duplicate(maxSampleCount);

        while(true) {
            int sampleCount = sampler.moreSamples(samples, rng);
            if(sampleCount < 1) {
                break;
            }
            for(int i = 0; i < sampleCount; i++) {
                Sample sample = samples[i];
                Ray intersectRay;
                camera.generateRay(sample, &intersectRay);

                double t0;
                double t1;
                if (!m_volumeRegion.intersectP(intersectRay, &t0, &t1)) {
                    continue;
                }
                if((t1-t0) == 0.0) {
                    continue;
                }

                Spectrum Tr(1.0);
                Spectrum Lv(0.1);

                Point3D p = intersectRay(t0);
                Ray startRay(p, intersectRay.m_direction);

//                startRay = Ray(startRay.origin() + intersectRay.m_direction * 0.01, startRay.direction());

                Integrator integrator(&m_volumeRegion, startRay, bounces, rng);

                integrator.integrate([&](const Ray& ray, photonflow::Length ds) {
                    Q_UNUSED(ds);
                    if(!m_volumeRegion.fuzzyInside(ray.origin())) {
                        return Integrator::Control::Break;
                    }
//                    double factor = (1.0_um - ds).value();
                    Tr *= m_volumeRegion.absorption(ray.origin(), Length3D(), 0.0);
                    if(m_volumeRegion.Density(ray.origin()) > 60) {
                        Lv += Tr * m_volumeRegion.emission(ray.origin(), Length3D(), 0.0);
                        photonflowAssert(!Lv.hasNaNs());
                    }
                    if(Tr < Spectrum(0.01)) {
                        return Integrator::Control::Break;
                    }
                    return Integrator::Control::Continue;
                });

                Spectrum final = Lv / omp_get_num_threads();

                m_film->addSample(sample, final);

                actualCount++;
                if(!(actualCount % 100000) && omp_get_thread_num() == 0) {
                    qDebug() << "Sample count:" << actualCount;
                }
            }
        }

        delete[] samples;
    }

    m_totalSampleCount += requestedSampleCount;

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {

            Pixel& pixel = (*m_film->pixels)(x, y);
            Spectrum result = Spectrum::fromXYZ(pixel.Lxyz);

            result /= (m_totalSampleCount * boxSize);

            double rgb[3];
            result.toRGB(rgb);

//            qDebug() << rgb[0] << rgb[1] << rgb[2];

            double factor = 1.0;
            QColor color(clamp(rgb[0]*factor, 0.0, 255.0),
                    clamp(rgb[1]*factor, 0.0, 255.0),
                    clamp(rgb[2]*factor, 0.0, 255.0),
                    255.0);

            m_image.setPixel(x, y, color.rgba());
        }
    }
    qDebug() << "Done after" << timer.elapsed() << "ms";
}

void PhotonflowWorker::synchronizeSimulator(Simulator *simulator)
{
    PhotonflowSimulator *renderView = qobject_cast<PhotonflowSimulator*>(simulator);
    if(renderView->m_clearRequested) {
        m_film.reset();
        m_totalSampleCount = 0;
        renderView->m_clearRequested = false;
        m_image.fill(QColor("black"));
    }

    m_volumeRegion.setAbsorptionCoefficient(Spectrum(renderView->m_absorptionCoefficient));
    m_volumeRegion.setScatteringCoefficient(Spectrum(renderView->m_scatteringCoefficient));
    m_volumeRegion.setEmissionCoefficient(Spectrum(renderView->m_emissionCoefficient));
    m_volumeRegion.setHenyeyGreensteinFactor(renderView->m_henyeyGreensteinFactor);

    renderView->m_image = m_image;

    emit renderView->imageChanged(renderView->m_image);
}

SimulatorWorker *PhotonflowSimulator::createWorker()
{
    return new PhotonflowWorker();
}

void PhotonflowSimulator::voxelize(const QVariantList &neuronSimulators)
{
    for(auto element : neuronSimulators) {
        QVariantMap map = element.toMap();
        NeuronSimulator *neuronSimulator = map["simulator"].value<NeuronSimulator*>();
        if(!neuronSimulator) {
            qWarning() << "Could not find simulator of object";
            return;
        }
        Qt3DCore::QTransform *transform = map["transform"].value<Qt3DCore::QTransform*>();
        if(!transform) {
            qWarning() << "Could not find transform of object";
            return;
        }
        photonflow::voxelize(neuronSimulator->cylinders(), neuronSimulator->boundingBox(), 1024);
    }
}

} // namespace
