#include "imageviewer.h"

#include <QPainter>

ImageViewer::ImageViewer(QQuickItem *parent)
    : QQuickPaintedItem(parent)
{

}


void ImageViewer::paint(QPainter *painter)
{
    painter->drawImage(boundingRect(), m_image);
}
