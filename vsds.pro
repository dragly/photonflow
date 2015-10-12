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
    samplers/random.h

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
    samplers/random.cpp

RESOURCES += qml.qrc

QMAKE_CXXFLAGS += -Wno-unused-parameter

LIBS += -lalglib

DISTFILES += \
    samplers/bestcandidate.out
