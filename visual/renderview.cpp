#include "renderview.h"

#include <QPainter>

RenderView::RenderView(QQuickItem *parent)
    : QQuickPaintedItem(parent)
{

}

void RenderView::integrate()
{
    QSize size = boundingRect().size().toSize();
    int w = size.width();
    int h = size.height();
    if(m_image.size() != size) {
        m_image = QImage(size, QImage::Format_ARGB32);
        for(int x = 0; x < size.width(); x++) {
            for(int y = 0; y < size.height(); y++) {
                m_image.setPixel(x, y, QColor(double(x) / w * 127.0, double(y) / h * 127.0, 127, 255).rgba());
            }
        }
    }
    update();
}

void RenderView::paint(QPainter *painter)
{
    painter->drawImage(0, 0, m_image);
}
