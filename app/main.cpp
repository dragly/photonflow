#include <QApplication>
#include <QQmlApplicationEngine>
#include "visual/renderview.h"

using namespace photonflow;

int main(int argc, char *argv[])
{
    qmlRegisterType<RenderView>("Photonflow", 1, 0, "RenderView");
    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));

    return app.exec();
}
