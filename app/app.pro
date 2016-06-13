TEMPLATE = app

QT = core gui qml quick widgets

include(../package.pri)
include(../vendor/SimVis/package_vendor.pri)
include(../conanbuildinfo.pri)

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

RESOURCES += \
    qml.qrc \
    meshes.qrc

SOURCES += \
    main.cpp \
    visual/renderview.cpp \
    visual/imageviewer.cpp \
    neuronsimulator.cpp

HEADERS += \
    visual/renderview.h \
    visual/imageviewer.h \
    neuronsimulator.h

