#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include "armadillo_includer.h"
#include "core/randomnumbergenerator.h"
#include "film/image.h"

#include <QQuickItem>
#include <QQuickPaintedItem>
#include <QImage>
#include <memory>

using std::unique_ptr;

class RenderView : public QQuickPaintedItem
{
    Q_OBJECT
public:
    RenderView(QQuickItem* parent = 0);
    Q_INVOKABLE void integrate();
    virtual void paint(QPainter *painter) override;
signals:

public slots:

private:
    QImage m_image;
    bool m_isDataLoaded = false;
    int times = 1;
    unique_ptr<ImageFilm> film;
    int totalSampleCount = 0;
    arma::cube data;
};

#endif // RENDERVIEW_H
