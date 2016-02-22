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

#include <QElapsedTimer>
#include <QPainter>
#include <iostream>
#include <memory>
#include <omp.h>

using namespace std;
using namespace arma;

RenderView::RenderView(QQuickItem *parent)
    : QQuickPaintedItem(parent)
{
    QElapsedTimer timer;
    timer.start();
    arma::Cube<short> data;
    data.load("/home/svenni/Dropbox/projects/programming/neuroscience/photonflow/photonflow/notebooks/output.hdf5", hdf5_binary);
//    data.load("/home/svenni/Dropbox/projects/programming/neuroscience/photonflow/photonflow/notebooks/volume.hdf5", hdf5_binary);
    qDebug() << "Data size:" << data.n_rows << data.n_cols << data.n_slices;
    qDebug() << "Data load time:" << timer.elapsed() << "ms";
    qDebug() << "Data max value: " << data.max();

//    data *= 255;

    BoundingBox bbox;
    photonflow::Length side = 100.0_um;
    bbox.pMin = Point3D(-side, -side, -0.2*side);
    bbox.pMax = Point3D(side, side, 0.2*side);

    double gg = 1.0;

    double angle = 0.5;

    Transform translation = translate(Length3D(0.1 * side, 0.3 * side, 1.0 * side));
    Transform rotation = rotate(angle, Vector3D(0.0, 1.0, 0.0));

    Transform boxTransform = translation*rotation;

//    Transform boxTransform;
    cout << "Identity: " << boxTransform.isIdentity() << endl;

    Spectrum sigma_a(0.95);
    Spectrum sigma_s(0.0);
    Spectrum emita(1.0);

    vr = VolumeGridDensity(sigma_a, sigma_s, gg, emita, bbox, boxTransform, data);
}

static const int threadCount = 4;

void RenderView::integrate()
{
    if(rngs.size() < threadCount) {
        rngs.resize(threadCount);
        int i = 0;
        for(RNG& rng : rngs) {
            rng.seed(1325125 ^ i);
            i++;
        }
    }
    QElapsedTimer timer;
    timer.start();

    const QSize size = boundingRect().size().toSize();
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

    if(m_image.size() != size || !film) {
        qDebug() << "Creating image of size" << size;
        m_image = QImage(size, QImage::Format_ARGB32);
        for(int x = 0; x < width; x++) {
            for(int y = 0; y < height; y++) {
                m_image.setPixel(x, y, QColor(0.0, 0.0, 0.0, 255.0).rgba());
            }
        }
        film = make_shared<ImageFilm>(width, height, &filter, crop);
    }
    const auto sopen = 0.0_us;
    const auto sclose = 1.0_us;
    const auto lensr = 0.0_um;
    const auto focald = 3.5_um;
    const double fov = 0.2;
    const PerspectiveCamera camera(cameraTransform, screenWindow, sopen, sclose, lensr, focald, fov, film);

#pragma omp parallel num_threads(threadCount) // OpenMP
    {
        cout << "Thread num: " << omp_get_thread_num() << " rngs size: " << rngs.size() << endl;
        RNG& rng = rngs.at(omp_get_thread_num());
        int actualCount = 0;
        RandomSampler sampler(0, width, 0, height, requestedSampleCount, 0.0_us, 1.0_us);
        int maxSampleCount = sampler.maximumSampleCount();
        qDebug() << "Max count: " << maxSampleCount;

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
                if (!vr.intersectP(intersectRay, &t0, &t1)) {
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

                Integrator integrator(&vr, startRay, bounces, rng);

                integrator.integrate([&](const Ray& ray, photonflow::Length ds) {
                    if(!vr.fuzzyInside(ray.origin())) {
                        return Integrator::Control::Break;
                    }
                    double factor = (1.0_um - ds).value();
                    Tr *= vr.sigma_a(ray.origin(), Length3D(), 0.0);
                    if(vr.Density(ray.origin()) > 60) {
                        Lv += Tr * vr.Lve(ray.origin(), Length3D(), 0.0);
                        photonflowAssert(!Lv.hasNaNs());
                    }
                    if(Tr < Spectrum(0.01)) {
                        return Integrator::Control::Break;
                    }
                    return Integrator::Control::Continue;
                });

                Spectrum final = Lv / omp_get_num_threads();

                film->addSample(sample, final);

                actualCount++;
                if(!(actualCount % 10000) && omp_get_thread_num() == 0) {
                    qDebug() << "Sample count:" << actualCount;
                }
            }
        }

        delete[] samples;
    }

    totalSampleCount += requestedSampleCount;

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {

            Pixel& pixel = (*film->pixels)(x, y);
            Spectrum result = Spectrum::fromXYZ(pixel.Lxyz);

            result /= (totalSampleCount * boxSize);

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
    update();
}

void RenderView::paint(QPainter *painter)
{
    qDebug() << "Paint!";
    painter->drawImage(0, 0, m_image);
    qDebug() << "Paint done!";
}
