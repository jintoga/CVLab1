QT += core
QT += gui
CONFIG += c++14

TARGET = Lab1
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    matrix.cpp \
    sobel.cpp \
    gauss.cpp \
    pointsofinterest.cpp \
    descriptors.cpp

HEADERS += \
    matrix.h \
    sobel.h \
    gauss.h \
    pointsofinterest.h \
    descriptors.h

