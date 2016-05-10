TEMPLATE = app

include(../../package.pri)
include(../../vendor/SimVis/package_vendor.pri)

QT += core gui quick widgets xml 3dcore 3drender

RESOURCES += \
    qml.qrc

HEADERS += \
    neuronsimulator.h

SOURCES += \
    main.cpp \
    neuronsimulator.cpp
