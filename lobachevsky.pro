TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    polycoeffproxy.cpp \
    lobachevskysolver.cpp \
    point.cpp \
    samplepoly.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    polycoeffproxy.h \
    lobachevskysolver.h \
    point.h \
    samplepoly.h

LIBS+=-lgmp

QMAKE_CXXFLAGS+=-std=c++11
