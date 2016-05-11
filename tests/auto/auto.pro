TEMPLATE = app

include(../../package.pri)

SOURCES += \
    $$PWD/henyeygreenstein_tests.cpp \
    $$PWD/henyeygreensteindistribution_tests.cpp \
    $$PWD/main.cpp

DISTFILES += \
    ../../.travis.yml
