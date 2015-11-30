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
    data.load("/home/svenni/Dropbox/projects/programming/neuroscience/photonflow/photonflow/notebooks/output.hdf5", hdf5_binary);
    qDebug() << "Data size:" << data.n_rows << data.n_cols << data.n_slices;
    qDebug() << "Data load time:" << timer.elapsed() << "ms";

    BBox bbox;
    bbox.pMin = Point(-1, -1, -0.3);
    bbox.pMax = Point(1, 1, 0.3);

    float gg = 1.0;

    float angle = -2.8;

    Transform translation = Translate(Vector(0.0, 0.0, 0.0));
    Transform rotation = Rotate(angle, Vector(0.0, 1.0, 0.0));

    Transform boxTransform = translation*rotation;

//    Transform boxTransform;
    cout << "Identity: " << boxTransform.IsIdentity() << endl;

    Spectrum sigma_a(0.97);
    Spectrum sigma_s(0.0);
    Spectrum emita(0.25);

    vr = VolumeGridDensity(sigma_a, sigma_s, gg, emita, bbox, boxTransform, data);
}

void RenderView::integrate()
{
    QElapsedTimer timer;
    timer.start();

    const QSize size = boundingRect().size().toSize();
    if(size.width() <= 0 || size.height() <= 0 || size.width() > 1e6 || size.height() > 1e6) {
        qWarning() << "WARNING: Integrate returns due to invalid size:" << size;
        return;
    }

    const int requestedSampleCount = 1;
    const int bounces = 200;
    const double ds = 0.01;

    const int width = size.width();
    const int height = size.height();

    const Transform cameraTransform = Translate(Vector(0, 0, -3.0));
    Rectangle screenWindow(-width / 2.0, -height / 2.0, width, height);
    const float crop[4] = {0.0, 1.0, 0.0, 1.0};

    BoxFilter filter(0.5, 0.5);

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
    const float sopen = 0.0;
    const float sclose = 1.0;
    const float lensr = 0.0;
    const float focald = 2.5;
    const float fov = 0.2;
    const PerspectiveCamera camera(cameraTransform, screenWindow, sopen, sclose, lensr, focald, fov, film);

#pragma omp parallel num_threads(8)       // OpenMP
    {
        RNG rng;
        rng.seed((1290481 ^ omp_get_thread_num()) + totalSampleCount);
        int actualCount = 0;
        RandomSampler sampler(0, width, 0, height, requestedSampleCount, 0.0, 1.0);
        int maxSampleCount = sampler.MaximumSampleCount();
        qDebug() << "Max count: " << maxSampleCount;

        Sample originalSample;
        originalSample.Add1D(1);
        Sample* samples = originalSample.Duplicate(maxSampleCount);

        while(true) {
            int sampleCount = sampler.GetMoreSamples(samples, rng);
            if(sampleCount < 1) {
                break;
            }
            for(int i = 0; i < sampleCount; i++) {
                Sample sample = samples[i];
                Ray intersectRay;
                camera.GenerateRay(sample, &intersectRay);

                float t0, t1;
                if (!vr.IntersectP(intersectRay, &t0, &t1) || (t1-t0) == 0.f) {
                    continue;
                }

                Spectrum Tr(1.0);
                Spectrum Lv(0.);

                Point p = intersectRay(t0);
                Ray startRay(p, intersectRay.m_direction);

                startRay = Ray(startRay.origin() + intersectRay.m_direction * 0.01, startRay.direction());

                Integrator integrator(&vr, startRay, bounces, rng);

                for(Ray& ray : integrator) {
                    if(!vr.inside(ray.origin())) {
                        break;
                    }
                    Tr *= vr.sigma_a(ray.origin(), Vector(), 0.0);
                    Lv += Tr * vr.Lve(ray.origin(), Vector(), 0.0);
                    if(Tr < Spectrum(0.01)) {
                        break;
                    }
                }

                Spectrum final = Lv / omp_get_num_threads();

                film->AddSample(sample, final);

                actualCount++;
                if(!(actualCount % 10000)) {
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
            Spectrum result = Spectrum::FromXYZ(pixel.Lxyz);

            result /= totalSampleCount;

            float rgb[3];
            result.ToRGB(rgb);

//            qDebug() << rgb[0] << rgb[1] << rgb[2];

            double factor = 1.0;
            QColor color(Clamp(rgb[0]*factor, 0.0, 255.0),
                    Clamp(rgb[1]*factor, 0.0, 255.0),
                    Clamp(rgb[2]*factor, 0.0, 255.0),
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
