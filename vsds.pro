TEMPLATE = app

QT += qml quick widgets

CONFIG += c++14

# Default rules for deployment.
include(deployment.pri)

HEADERS += \
    visual/renderview.h \
    core/randomnumbergenerator.h \
    core/transform.h \
    core/common.h \
    core/geometry.h \
    core/quaternion.h \
    core/error.h \
    core/camera.h \
    core/film.h \
    core/sampler.h \
    core/montecarlo.h \
    cameras/perspective.h \
    core/spectrum.h \
    film/image.h \
    core/filter.h \
    samplers/random.h \
    volumes/volumegrid.h \
    core/volume.h \
    core/integrator.h \
    core/renderer.h \
    core/scene.h \
    core/intersection.h \
    core/diffgeom.h \
    core/light.h \
    core/shape.h \
    core/memory.h \
    filters/box.h \
    filters/gaussian.h \
    filters/mitchell.h \
    filters/sinc.h \
    filters/triangle.h \
    stdafx.h

SOURCES += main.cpp \
    visual/renderview.cpp \
    core/randomnumbergenerator.cpp \
    core/transform.cpp \
    core/geometry.cpp \
    core/quaternion.cpp \
    core/error.cpp \
    core/camera.cpp \
    core/film.cpp \
    core/sampler.cpp \
    core/montecarlo.cpp \
    cameras/perspective.cpp \
    core/spectrum.cpp \
    film/image.cpp \
    core/filter.cpp \
    samplers/random.cpp \
    volumes/volumegrid.cpp \
    core/volume.cpp \
    core/integrator.cpp \
    core/renderer.cpp \
    core/scene.cpp \
    core/intersection.cpp \
    core/diffgeom.cpp \
    core/light.cpp \
    core/shape.cpp \
    core/memory.cpp \
    filters/box.cpp \
    filters/gaussian.cpp \
    filters/mitchell.cpp \
    filters/sinc.cpp \
    filters/triangle.cpp

RESOURCES += qml.qrc

QMAKE_CXXFLAGS += -Wno-unused-parameter

#LIBS += -lalglib

DISTFILES += \
    samplers/bestcandidate.out
