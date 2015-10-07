#include "renderview.h"

#include <QPainter>
#include <interpolation.h>
#include <iostream>

#include "../cameras/perspective.h"
#include "../film/image.h"
#include "../core/randomnumbergenerator.h"

using std::cout; using std::endl;

RenderView::RenderView(QQuickItem *parent)
    : QQuickPaintedItem(parent)
{
}

void RenderView::integrate()
{
    RandomNumberGenerator rng;
    rng.seed(0);
    Transform a;
    Transform b;
    float screenWindow[4];
    screenWindow[0] = 0.0;
    screenWindow[1] = 200.0;
    screenWindow[2] = 0.0;
    screenWindow[3] = 200.0;
    ImageFilm film(200, 200, nullptr, nullptr, string(), true);
    PerspectiveCamera camera(AnimatedTransform(&a, 0.0, &b, 0.0), screenWindow, 0.0, 1.0, 1.0, 1.0, 0.0, &film);

    QSize size = boundingRect().size().toSize();
    int w = size.width();
    int h = size.height();
    if(m_image.size() != size) {
        m_image = QImage(size, QImage::Format_ARGB32);
        for(int x = 0; x < size.width(); x++) {
            for(int y = 0; y < size.height(); y++) {
                CameraSample sample;
                sample.imageX = x;
                sample.imageY = y;
                sample.lensU = rng.RandomFloat();
                sample.lensV = rng.RandomFloat();
                sample.time = 0.5;
                Ray ray;
                camera.GenerateRay(sample, &ray);

                m_image.setPixel(x, y, QColor(ray.d.x * 127.0, ray.d.y * 127.0, 127, 255).rgba());
            }
        }
    }
    update();
}

void RenderView::paint(QPainter *painter)
{
    painter->drawImage(0, 0, m_image);
}
