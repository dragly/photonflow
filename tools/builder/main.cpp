#include <QApplication>
#include <QQmlApplicationEngine>

#include "neuronsimulator.h"
#include <vendor.h>

int main(int argc, char *argv[])
{
    qmlRegisterType<NeuronSimulator>("Photonflow", 1, 0, "NeuronSimulator");
    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    qpm::init(app, engine);
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));

    return app.exec();

}

