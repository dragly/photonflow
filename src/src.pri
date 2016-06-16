DISTFILES +=

QT += core gui opengl widgets xml

CONFIG += c++14

INCLUDEPATH = $$PWD

DEFINES += ARMA_USE_HDF5

trusty {
    LIBS += -lhdf5
} else {
    LIBS += -lhdf5_serial
    DEFINES += ARMA_HDF5_INCLUDE_DIR=/usr/include/hdf5/serial/
}

HEADERS += \
    $$PWD/cameras/perspective.h \
    $$PWD/core/camera.h \
    $$PWD/core/common.h \
    $$PWD/core/diffgeom.h \
    $$PWD/core/error.h \
    $$PWD/core/film.h \
    $$PWD/core/filter.h \
    $$PWD/core/geometry.h \
    $$PWD/core/heyneygreenstein.h \
    $$PWD/core/integrator.h \
    $$PWD/core/intersection.h \
    $$PWD/core/memory.h \
    $$PWD/core/montecarlo.h \
    $$PWD/core/quaternion.h \
    $$PWD/core/randomnumbergenerator.h \
    $$PWD/core/renderer.h \
    $$PWD/core/sampler.h \
    $$PWD/core/scene.h \
    $$PWD/core/shape.h \
    $$PWD/core/spectrum.h \
    $$PWD/core/transform.h \
    $$PWD/core/units.h \
    $$PWD/core/volume.h \
    $$PWD/film/image.h \
    $$PWD/filters/box.h \
    $$PWD/filters/gaussian.h \
    $$PWD/filters/mitchell.h \
    $$PWD/filters/sinc.h \
    $$PWD/filters/triangle.h \
    $$PWD/geometry/bbox.h \
    $$PWD/geometry/cylinderfrustum.h \
    $$PWD/geometry/normal.h \
    $$PWD/geometry/point3d.h \
    $$PWD/geometry/ray.h \
    $$PWD/geometry/raydifferential.h \
    $$PWD/geometry/rectangle.h \
    $$PWD/geometry/vector3d.h \
    $$PWD/io/neuromlreader.h \
    $$PWD/io/voxelizer.h \
    $$PWD/samplers/random.h \
    $$PWD/volumes/volumegrid.h \
    $$PWD/armadillo_includer.h \
    $$PWD/stdafx.h \
    $$PWD/cameras/orthographiccamera.h

SOURCES += \
    $$PWD/cameras/perspective.cpp \
    $$PWD/core/camera.cpp \
    $$PWD/core/diffgeom.cpp \
    $$PWD/core/error.cpp \
    $$PWD/core/film.cpp \
    $$PWD/core/filter.cpp \
    $$PWD/core/geometry.cpp \
    $$PWD/core/integrator.cpp \
    $$PWD/core/intersection.cpp \
    $$PWD/core/memory.cpp \
    $$PWD/core/montecarlo.cpp \
    $$PWD/core/quaternion.cpp \
    $$PWD/core/randomnumbergenerator.cpp \
    $$PWD/core/renderer.cpp \
    $$PWD/core/sampler.cpp \
    $$PWD/core/scene.cpp \
    $$PWD/core/shape.cpp \
    $$PWD/core/spectrum.cpp \
    $$PWD/core/transform.cpp \
    $$PWD/core/units.cpp \
    $$PWD/core/volume.cpp \
    $$PWD/film/image.cpp \
    $$PWD/filters/box.cpp \
    $$PWD/filters/gaussian.cpp \
    $$PWD/filters/mitchell.cpp \
    $$PWD/filters/sinc.cpp \
    $$PWD/filters/triangle.cpp \
    $$PWD/geometry/bbox.cpp \
    $$PWD/geometry/cylinderfrustum.cpp \
    $$PWD/geometry/normal.cpp \
    $$PWD/geometry/point3d.cpp \
    $$PWD/geometry/ray.cpp \
    $$PWD/geometry/raydifferential.cpp \
    $$PWD/geometry/rectangle.cpp \
    $$PWD/geometry/vector3d.cpp \
    $$PWD/io/neuromlreader.cpp \
    $$PWD/io/voxelizer.cpp \
    $$PWD/samplers/random.cpp \
    $$PWD/volumes/volumegrid.cpp \
    $$PWD/cameras/orthographiccamera.cpp
