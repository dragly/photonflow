#include <QApplication>
#include <QQmlApplicationEngine>
#include <QSurfaceFormat>
#include <vendor.h>

#include "neuronsimulator.h"
#include "visual/renderview.h"
#include "visual/imageviewer.h"

using namespace photonflow;

int main(int argc, char *argv[])
{
//    qputenv("QSG_RENDER_LOOP", "basic"); // TODO remove this when fixed in Qt3D

    QSurfaceFormat format = QSurfaceFormat::defaultFormat();
    qDebug() << format.samples();
    format.setSamples(32);
    QSurfaceFormat::setDefaultFormat(format);

    qmlRegisterType<PhotonflowSimulator>("Photonflow", 1, 0, "PhotonflowSimulator");
    qmlRegisterType<NeuronSimulator>("Photonflow", 1, 0, "NeuronSimulator");
    qmlRegisterType<ImageViewer>("Photonflow", 1, 0, "ImageViewer");

    QApplication::setOrganizationName("Ovilab");
    QApplication::setApplicationName("Photonflow");

    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    qpm::init(app, engine);
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));

    return app.exec();
}

