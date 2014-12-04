TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES +=  Main.cpp \
	    Cortical_Column.cpp \
	    Thalamic_Column.cpp \
	    TC.cpp

HEADERS +=  macros.h \
	    ODE.h \
	    saves.h \
	    Cortical_Column.h \
	    Thalamic_Column.h \
	    Stimulation.h

QMAKE_CXXFLAGS += -std=c++11 -O3

SOURCES -= TC.cpp
