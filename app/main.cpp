#include <QApplication>
#include <QQmlApplicationEngine>
#include "visual/renderview.h"

int main(int argc, char *argv[])
{
    qmlRegisterType<RenderView>("VSDS", 1, 0, "RenderView");
    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));

    return app.exec();
}

