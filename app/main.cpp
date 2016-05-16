#include <QApplication>
#include <QQmlApplicationEngine>
#include <vendor.h>

#include "visual/renderview.h"
#include "visual/imageviewer.h"

using namespace photonflow;

int main(int argc, char *argv[])
{
    qmlRegisterType<PhotonflowSimulator>("Photonflow", 1, 0, "PhotonflowSimulator");
    qmlRegisterType<ImageViewer>("Photonflow", 1, 0, "ImageViewer");

    QApplication::setOrganizationName("Ovilab");
    QApplication::setApplicationName("Photonflow");

    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    qpm::init(app, engine);
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));

    return app.exec();

}

