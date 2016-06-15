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
#include "neuronsimulator.h"

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

static const int threadCount = 1;

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

double PhotonflowSimulator::renderTime() const
{
    return m_renderTime;
}

Qt3DCore::QEntity *PhotonflowSimulator::camera() const
{
    return m_camera;
}

int PhotonflowSimulator::completedSampleCount() const
{
    return m_completedSampleCount;
}

double PhotonflowSimulator::fieldOfView() const
{
    return m_fieldOfView;
}

double PhotonflowSimulator::lensRadius() const
{
    return m_lensRadius;
}

double PhotonflowSimulator::focalDepth() const
{
    return m_focalDepth;
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

    const Transform cameraTransform = translate(Length3D(0.0_um, 0.0_um, -512.0_um));
    Rectangle screenWindow(-width / 2.0, -height / 2.0, width, height);
    const double crop[4] = {0.0, 1.0, 0.0, 1.0};

    double boxSize = 2.0;
    BoxFilter filter(boxSize * 0.5, boxSize * 0.5);

    if(m_image.size() != size || !m_film) {
        qDebug() << "Creating image of size" << size;
        m_image = QImage(size, QImage::Format_ARGB32);
        m_image.fill(Qt::black);
        m_film = make_shared<ImageFilm>(width, height, &filter, crop);
    }
    const auto sopen = 0.0_us;
    const auto sclose = 1.0_us;
    const PerspectiveCamera camera(cameraTransform, screenWindow, sopen, sclose, m_lensRadius, m_focalDepth, m_fieldOfView, m_film);

    int actualCount = 0;
#pragma omp parallel num_threads(threadCount) // OpenMP
    {
        RNG& rng = m_randomNumberGenerators.at(omp_get_thread_num());
        //        RandomSampler sampler(0, width, 0, height, requestedSampleCount, 0.0_us, 1.0_us);
        //        int maxSampleCount = sampler.maximumSampleCount();

        //        Sample originalSample;
        //        originalSample.Add1D(1);
        //        Sample* samples = originalSample.Duplicate(maxSampleCount);

        //        while(true) {
        //            int sampleCount = sampler.moreSamples(samples, rng);
        //            if(sampleCount < 1) {
        //                break;
        //            }
        //            for(int i = 0; i < sampleCount; i++) {
        //                Sample sample = samples[i];

        for(int i = 0; i < width * height; i++) {
            Sample sample;
            //                sample.Add1D(1);
            //                cout << sample.imageX << " " << sample.imageY << " " << sample.lensU << " " << sample.lensV << " " << sample.time.value() << endl;

            sample.imageX = rng.randomFloat() * width;
            sample.imageY = rng.randomFloat() * height;
            sample.lensU = rng.randomFloat();
            sample.lensV = rng.randomFloat();
            sample.time = 0.0_us; // TODO allow for different samplers

            Ray intersectRay;
            camera.generateRay(sample, &intersectRay);
            //                cout << intersectRay.origin().x << " "<< intersectRay.origin().y << " "<< intersectRay.origin().z << " " << intersectRay.direction().x << " " << intersectRay.direction().y << " " << intersectRay.direction().z << endl;

            double t0;
            double t1;
            if (!m_boundingBox.intersectP(intersectRay, &t0, &t1)) {
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

            Integrator integrator(startRay, bounces, rng);


            const Length cellSide = 10_um;
            const Length3D boundingBoxSize = m_boundingBox.pMax - m_boundingBox.pMin;

            size_t cellCount[3] = {
                size_t(boundingBoxSize[0] / cellSide) + 1,
                size_t(boundingBoxSize[1] / cellSide) + 1,
                size_t(boundingBoxSize[2] / cellSide) + 1
            };

            integrator.integrate([&](const Ray& ray, photonflow::Length ds) {
                if(!m_boundingBox.fuzzyInside(ray.origin())) {
                    return Integrator::Control::Break;
                }

                //                    Lv += Tr * 1.0 * ds.value();

                Length3D pointInBoundingBox = ray.origin() - m_boundingBox.pMin;
                int i = pointInBoundingBox[0] / cellSide;
                int j = pointInBoundingBox[1] / cellSide;
                int k = pointInBoundingBox[2] / cellSide;

                int index = i * cellCount[1] * cellCount[2] + j * cellCount[2] + k;


                if(index >= 0 && index < m_cells.size()) {
                    Cell &cell = m_cells[index];
                    for(const CylinderFrustum &cylinder : cell.m_cylinders) {
                        Length eps = 0.0_um;
                        Point3D p = ray.origin();
                        Length3D diff = p - cylinder.center;
                        Length distance = diff.length();
                        if(distance > cylinder.h + eps && distance > cylinder.startRadius + eps) {
                            continue;
                        }
                        auto yComponent = dot(diff, cylinder.direction * 1.0_um) / 1.0_um;
                        if(fabs(yComponent) <= cylinder.h + eps) {
                            auto y2 = yComponent*yComponent;
                            auto diff2 = dot(diff, diff);
                            auto distanceToAxis = sqrt(diff2 - y2);
                            double endProportion = (yComponent + cylinder.h) / (2.0 * cylinder.h);
                            Length radius = cylinder.startRadius * (1 - endProportion) + endProportion * cylinder.endRadius;
                            if(distanceToAxis <= radius + eps) {
                                // TODO replace with proper emission
                                Lv += Tr * 20.0 * ds.value();
                            }
                        }
                    }
                }

                if(Tr < Spectrum(0.01)) {
                    return Integrator::Control::Break;
                }
                return Integrator::Control::Continue;
            });

            Spectrum final = Lv / omp_get_num_threads();

            m_film->addSample(sample, final);

            actualCount++;
        }
        //        }

        //        delete[] samples;
    }

    m_completedSampleCount += actualCount;
    m_iterationCount += requestedSampleCount;

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {

            Pixel& pixel = (*m_film->pixels)(x, y);
            Spectrum result = Spectrum::fromXYZ(pixel.Lxyz);

            result /= pixel.weightSum;

            double rgb[3];
            result.toRGB(rgb);

            //            qDebug() << rgb[0] << rgb[1] << rgb[2];

            double factor = 1.0;
            QColor color(clamp(rgb[0]*factor, 0.0, 255.0),
                    clamp(rgb[1]*factor, 0.0, 255.0),
                    clamp(rgb[2]*factor, 0.0, 255.0),
                    255.0);

            m_image.setPixel(x, height - y - 1, color.rgba());
        }
    }
    m_renderTime = timer.elapsed();
}

const Length sideLength = 256.0_um;

void PhotonflowWorker::synchronizeSimulator(Simulator *simulator)
{
    PhotonflowSimulator *renderView = qobject_cast<PhotonflowSimulator*>(simulator);
    if(renderView->m_clearRequested) {
        m_film.reset();
        m_completedSampleCount = 0;
        m_iterationCount = 0;
        renderView->m_clearRequested = false;
        m_image.fill(QColor("black"));
    }

    if(renderView->m_dataDirty) {
        BoundingBox bbox(Point3D(-sideLength, -sideLength, -sideLength), Point3D(sideLength, sideLength, sideLength));
        double gg = 1.0;
        double angle = 0.0;
        photonflow::Length side = 100.0_um;
        Transform translation = translate(Length3D(0.0 * side, 0.0 * side, 0.0 * side));
        Transform boxTransform = translation;
        Spectrum sigma_a(0.95);
        Spectrum sigma_s(0.1);
        Spectrum emita(0.1);
        m_volumeRegion = VolumeGridDensity(sigma_a, sigma_s, gg, emita, bbox, boxTransform, renderView->m_data);

        m_cylinders = renderView->m_cylinders;
        m_focalDepth = renderView->m_focalDepth * 1.0_um;
        m_lensRadius = renderView->m_lensRadius * 1.0_um;
        m_fieldOfView = renderView->m_fieldOfView;

        m_boundingBox = BoundingBox();
        for(const auto& cylinder : m_cylinders) {
            m_boundingBox = makeUnion(m_boundingBox, cylinder.start);
            m_boundingBox = makeUnion(m_boundingBox, cylinder.end);
        }

        const Length cellSide = 10_um;
        const Length3D boundingBoxSize = m_boundingBox.pMax - m_boundingBox.pMin;

        if(boundingBoxSize[0] > 0.0_um && boundingBoxSize[1] > 0.0_um && boundingBoxSize[1] > 0.0_um) {
            size_t cellCount[3] = {
                size_t(boundingBoxSize[0] / cellSide) + 1,
                size_t(boundingBoxSize[1] / cellSide) + 1,
                size_t(boundingBoxSize[2] / cellSide) + 1
            };

            size_t cellCountTotal = cellCount[0] * cellCount[1] * cellCount[2];
            m_cells.clear();
            m_cells.resize(cellCountTotal); // TODO empty

            Length step = cellSide;
            Length eps = step; // / 2.0;

            Length3D offset(m_boundingBox.pMin);

            for(auto& originalCylinder : m_cylinders) {
                auto cylinder = CylinderFrustum(originalCylinder.start - offset,
                                                originalCylinder.end - offset,
                                                originalCylinder.startRadius,
                                                originalCylinder.endRadius);

                BoundingBox localBounds(Point3D(-cylinder.h, -cylinder.startRadius, -cylinder.startRadius),
                                        Point3D(cylinder.h, cylinder.startRadius, cylinder.startRadius));

                Vector3D perpendicular = cross(Vector3D(1.0, 0.0, 0.0), cylinder.direction);
                Transform rotation;
                double sinAngle = perpendicular.length();
                if(sinAngle > 0.0) {
                    if(sinAngle > 1.0) {
                        sinAngle -= 2.2204460492503131e-16; // typical machine precision error
                    }
                    double cosAngle = sqrt(1.0 - sinAngle*sinAngle);

                    rotation = rotate(cosAngle, sinAngle, perpendicular);
                }
                Transform translation = translate(Length3D(cylinder.center));
                BoundingBox bounds = translation(rotation(localBounds));
                bounds.expand(eps);

                int istart = int(bounds.pMin.x / step);
                int jstart = int(bounds.pMin.y / step);
                int kstart = int(bounds.pMin.z / step);

                int iend = int(bounds.pMax.x / step + 1);
                int jend = int(bounds.pMax.y / step + 1);
                int kend = int(bounds.pMax.z / step + 1);

                for(int i = istart; i < iend + 1; i++) {
                    for(int j = jstart; j < jend + 1; j++) {
                        for(int k = kstart; k < kend + 1; k++) {
                            Point3D p(step * (double(i) + 0.5), step * (double(j) + 0.5), step * (double(k) + 0.5));
                            Length3D diff = p - cylinder.center;
                            Length distance = diff.length();
                            if(distance > cylinder.h + eps && distance > cylinder.startRadius + eps) {
                                continue;
                            }
                            auto yComponent = dot(diff, cylinder.direction * 1.0_um) / 1.0_um;
                            if(fabs(yComponent) <= cylinder.h + eps) {
                                auto y2 = yComponent*yComponent;
                                auto diff2 = dot(diff, diff);
                                auto distanceToAxis = sqrt(diff2 - y2);
                                double endProportion = (yComponent + cylinder.h) / (2.0 * cylinder.h);
                                Length radius = cylinder.startRadius * (1 - endProportion) + endProportion * cylinder.endRadius;
                                if(distanceToAxis <= radius + eps) {
                                    size_t index = i * cellCount[1] * cellCount[2] + j * cellCount[2] + k;
                                    if(index < m_cells.size()) {
                                        m_cells[index].m_cylinders.push_back(originalCylinder);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        renderView->m_dataDirty = false;
    }

    //    m_volumeRegion.setAbsorptionCoefficient(Spectrum(renderView->m_absorptionCoefficient));
    //    m_volumeRegion.setScatteringCoefficient(Spectrum(renderView->m_scatteringCoefficient));
    //    m_volumeRegion.setEmissionCoefficient(Spectrum(renderView->m_emissionCoefficient));
    //    m_volumeRegion.setHenyeyGreensteinFactor(renderView->m_henyeyGreensteinFactor);

    if(renderView->m_renderTime != m_renderTime) {
        renderView->m_renderTime = m_renderTime;
        emit renderView->renderTimeChanged(m_renderTime);
    }
    if(renderView->m_completedSampleCount != m_completedSampleCount) {
        renderView->m_completedSampleCount = m_completedSampleCount;
        emit renderView->completedSampleCountChanged(m_completedSampleCount);
    }

    renderView->m_image = m_image;
    emit renderView->imageChanged(renderView->m_image);
}

SimulatorWorker *PhotonflowSimulator::createWorker()
{
    return new PhotonflowWorker();
}

void PhotonflowSimulator::voxelize(const QVariantList &neuronSimulators)
{
    m_cylinders.clear();
    for(auto element : neuronSimulators) {
        QVariantMap map = element.toMap();
        NeuronSimulator *neuronSimulator = map["simulator"].value<NeuronSimulator*>();
        if(!neuronSimulator) {
            qWarning() << "Could not find simulator of object";
            return;
        }
        Qt3DCore::QTransform *transform3d = map["transform"].value<Qt3DCore::QTransform*>();
        if(!transform3d) {
            qWarning() << "Could not find transform of object";
            return;
        }

        QMatrix4x4 t = transform3d->matrix();

        // TODO verify order
        Matrix4x4 matrix(
                    t(0, 0), t(0, 1), t(0, 2), t(0, 3)*(1.0 / neuronSimulator->scale()),
                    t(1, 0), t(1, 1), t(1, 2), t(1, 3)*(1.0 / neuronSimulator->scale()),
                    t(2, 0), t(2, 1), t(2, 2), t(2, 3)*(1.0 / neuronSimulator->scale()),
                    t(3, 0), t(3, 1), t(3, 2), t(3, 3));
        Transform transform(matrix);

        std::vector<CylinderFrustum> cylinders = neuronSimulator->cylinders();
        for(CylinderFrustum& cylinder : cylinders) {
            cylinder = CylinderFrustum(transform(cylinder.start),
                                       transform(cylinder.end),
                                       cylinder.startRadius,
                                       cylinder.endRadius);
        }
        m_cylinders.insert(m_cylinders.end(), cylinders.begin(), cylinders.end());
    }
    m_dataDirty = true;
    clear();
}

void PhotonflowSimulator::setCamera(Qt3DCore::QEntity *camera)
{
    if (m_camera == camera)
        return;

    m_camera = camera;
    emit cameraChanged(camera);
}

void PhotonflowSimulator::setFieldOfView(double fieldOfView)
{
    if (m_fieldOfView == fieldOfView)
        return;

    m_fieldOfView = fieldOfView;
    emit fieldOfViewChanged(fieldOfView);
}

void PhotonflowSimulator::setLensRadius(double lensRadius)
{
    if (m_lensRadius == lensRadius)
        return;

    m_lensRadius = lensRadius;
    emit lensRadiusChanged(lensRadius);
}

void PhotonflowSimulator::setFocalDepth(double focalDepth)
{
    if (m_focalDepth == focalDepth)
        return;

    m_focalDepth = focalDepth;
    emit focalDepthChanged(focalDepth);
}

} // namespace
