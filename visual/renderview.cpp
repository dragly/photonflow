#include "renderview.h"

#include <QPainter>
#include <interpolation.h>
#include <iostream>
#include <memory>

#include "../cameras/perspective.h"
#include "../film/image.h"
#include "../core/randomnumbergenerator.h"
#include "../samplers/random.h"

using std::cout; using std::endl;
using std::unique_ptr; using std::make_unique;

RenderView::RenderView(QQuickItem *parent)
    : QQuickPaintedItem(parent)
{
}

void RenderView::integrate()
{
    RNG rng;
    rng.seed(0);
    Transform a;
    Transform b;
    float screenWindow[4];
    screenWindow[0] = 0.0;
    screenWindow[1] = 200.0;
    screenWindow[2] = 0.0;
    screenWindow[3] = 200.0;
    ImageFilm film(200, 200, nullptr, nullptr, string(), true);
    PerspectiveCamera camera(AnimatedTransform(&a, 0.0, &b, 0.0), screenWindow, 0.0, 1.0, 0.0, 1.0, 1.0, &film);

    QSize size = boundingRect().size().toSize();

    int sampleCount = 1;

    int width = size.width();
    int height = size.height();
    if(m_image.size() != size) {
        m_image = QImage(size, QImage::Format_ARGB32);
        RandomSampler sampler(0, width, 0, height, sampleCount, 0.0, 1.0);
        int maxSampleCount = sampler.MaximumSampleCount();

        Sample origSample(&sampler);
        Sample* samples = origSample.Duplicate(maxSampleCount);
        int sampleCount = 0;
        while((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
            for(int i = 0; i < sampleCount; i++) {
                Sample sample = samples[i];
                Ray ray;
                camera.GenerateRay(sample, &ray);
                double factor = 255.0;
                QColor color(Clamp(ray.d.x*factor, 0.0, 255.0),
                             Clamp(ray.d.y*factor, 0.0, 255.0),
                             Clamp(ray.d.z*factor, 0.0, 255.0),
                             255.0);
                m_image.setPixel(sample.imageX, sample.imageY, color.rgba());
            }
        }
    }
    qDebug() << "Done!";
    update();
}

void RenderView::paint(QPainter *painter)
{
    painter->drawImage(0, 0, m_image);
}
