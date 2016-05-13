#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QQuickPaintedItem>
#include <QImage>

class ImageViewer : public QQuickPaintedItem
{
    Q_OBJECT
    Q_PROPERTY(QImage image READ image WRITE setImage NOTIFY imageChanged)
public:
    ImageViewer(QQuickItem *parent = nullptr);

    QImage image()
    {
        return m_image;
    }

signals:

    void imageChanged(QImage image);

public slots:

    void setImage(QImage image)
    {
        if (m_image == image)
            return;

        m_image = image;
        emit imageChanged(image);
        update();
    }

private:
    QImage m_image;

    // QQuickPaintedItem interface
public:
    virtual void paint(QPainter *painter) override;
};

#endif // IMAGEVIEWER_H
