#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include "armadillo_includer.h"
#include "core/randomnumbergenerator.h"
#include "film/image.h"
#include "volumes/volumegrid.h"

#include <QQuickItem>
#include <QQuickPaintedItem>
#include <QImage>
#include <memory>

namespace photonflow {

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
    std::shared_ptr<ImageFilm> film;
    int totalSampleCount = 0;
    VolumeGridDensity vr;
    vector<RNG> rngs;
};

}

#endif // RENDERVIEW_H
