#include "renderview.h"

#include <QElapsedTimer>
#include <QPainter>
#include <iostream>
#include <memory>
#include <omp.h>

#include "cameras/perspective.h"
#include "film/image.h"
#include "core/randomnumbergenerator.h"
#include "samplers/random.h"
#include "volumes/volumegrid.h"
#include "core/heyneygreenstein.h"

#include "filters/box.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"

using std::cout; using std::endl;
using std::unique_ptr; using std::make_unique;

RenderView::RenderView(QQuickItem *parent)
    : QQuickPaintedItem(parent)
{
}

void RenderView::integrate()
{
    QElapsedTimer timer;
    timer.start();
    QSize size = boundingRect().size().toSize();

    int requestedSampleCount = 1;
    int bounces = 500;
    double ds = 0.01;

    int width = size.width();
    int height = size.height();

    Transform identity({{1.0, 0.0, 0.0, 0.0,
                         0.0, 1.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 0.0, 0.0, 1.0}});

    Transform a({{1.0, 0.0, 0.0, 0.0,
                  0.0, 1.0, 0.0, 0.0,
                  0.0, 0.0, 1.0, 0.0,
                  0.0, 0.0, 0.0, 1.0}});
    Transform b = a;
    float screenWindow[4];
    screenWindow[0] = -width / 2.0;
    screenWindow[1] = width / 2.0;
    screenWindow[2] = - height / 2.0;
    screenWindow[3] = height / 2.0;
    float crop[4];
    crop[0] = 0.0;
    crop[1] = 1.0;
    crop[2] = 0.0;
    crop[3] = 1.0;

    BoxFilter filter(0.5, 0.5);
//    MitchellFilter filter(1.0/3.0, 1.0/3.0, 2.0, 2.0);
//    LanczosSincFilter filter(4.0, 4.0, 3.0);


    if(!film) {
        film = make_unique<ImageFilm>(width, height, &filter, crop);
    }
    float sopen = 0.0;
    float sclose = 1.0;
    float lensr = 0.0;
    float focald = 2.5;
    float fov = 0.3;
    PerspectiveCamera camera(AnimatedTransform(&a, 0.0, &b, 0.0), screenWindow,
                             sopen, sclose, lensr, focald, fov, film.get());

    BBox bbox;
    bbox.pMin = Point(-1, -1, -1);
    bbox.pMax = Point(1, 1, 1);
    Transform boxTransform = identity;
    int nx = 3;
    int ny = 3;
    int nz = 3;
    float data[27] = {1.0, 0.1, 1.0,
                      0.1, 1.0, 0.1,
                      1.0, 0.1, 1.0,
                      0.1, 1.0, 0.1,
                      1.0, 0.1, 1.0,
                      0.1, 1.0, 0.1,
                      1.0, 0.1, 1.0,
                      0.1, 1.0, 0.1,
                      1.0, 0.1, 1.0};
    float gg = 1.0;

    float angle = 0.6;
    float ca = cos(angle);
    float sa = sin(angle);
    Transform rotate({ca,   0.0,    sa,     0.0,
                      0.0,          1.0,    0.0,            0.0,
                      -sa,  0.0,    ca,     0.0,
                      0.0,          0.0,    0.0,            1.0});
    Transform translate({{1.0, 0.0, 0.0, 0.0,
                          0.0, 1.0, 0.0, 0.0,
                          0.0, 0.0, 1.0, 3.0,
                          0.0, 0.0, 0.0, 1.0}});

    boxTransform = translate*rotate*boxTransform;

    Spectrum sigma_a(0.99);
    Spectrum sigma_s(0.0);
    Spectrum emita(10.0);

    VolumeGridDensity vr(sigma_a, sigma_s, gg, emita, bbox, boxTransform, nx, ny, nz, data);

    if(m_image.size() != size) {
        m_image = QImage(size, QImage::Format_ARGB32);
        for(int x = 0; x < width; x++) {
            for(int y = 0; y < height; y++) {
                m_image.setPixel(x, y, QColor(0.0, 0.0, 0.0, 255.0).rgba());
            }
        }
    }
#pragma omp parallel num_threads(8)       // OpenMP
    {
        RNG rng;
        rng.seed((1290481 ^ omp_get_thread_num()) + totalSampleCount);
        int actualCount = 0;
        RandomSampler sampler(0, width, 0, height, requestedSampleCount, 0.0, 1.0);
        int maxSampleCount = sampler.MaximumSampleCount();
        qDebug() << "Max count: " << maxSampleCount;

        Sample origSample(&sampler);
        int scatterSampleOffset = origSample.Add1D(1);
        Sample* samples = origSample.Duplicate(maxSampleCount);
        int sampleCount = 0;

        while((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
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
                Ray ray(p, intersectRay.d, 0.0);

                float t = t0;
                for(int i = 0; i < bounces; i++) {
//                    t += ds;
                    double g = 0.98;
//                    double theta = acos(Distribution::heyneyGreenstein(g, rng));
                    double cosTheta = Distribution::heyneyGreenstein(g, rng);
                    double sinTheta = sqrt(1 - cosTheta*cosTheta);
                    double phi = 2.0 * M_PI * rng.RandomFloat();

                    Vector perpendicular = ray.d.perpendicular();
                    Transform perpendicularRotation = Rotate(phi, ray.d);
                    perpendicular = perpendicularRotation(perpendicular);

                    Transform directionRotation = Rotatec(cosTheta, sinTheta, perpendicular);
                    ray.d = directionRotation(ray.d);
                    ray.d = ray.d.normalized();

                    ray.o += ray.d * ds;

                    if(vr.Density(ray.o) <= 0.0) {
                        break;
                    }
                    Tr *= sigma_a;
                    Lv += Tr * vr.Lve(ray.o, Vector(), 0.0);
                    if(Tr < Spectrum(0.01)) {
                        break;
                    }
                }
                Spectrum final = Lv / omp_get_num_threads();
                final *= 1.0;

                film->AddSample(sample, final);

                actualCount++;
                if(!(actualCount % 10000)) {
                    qDebug() << "Sample count:" << actualCount;
                }
            }
        }
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
