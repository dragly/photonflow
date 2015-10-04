#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include <QQuickItem>
#include <QQuickPaintedItem>
#include <QImage>

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
};

#endif // RENDERVIEW_H
