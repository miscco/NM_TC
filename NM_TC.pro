TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = TC.cpp

SOURCES +=  Cortical_Column.cpp \
	    TC.cpp		\
	    TC_mex.cpp		\
	    Thalamic_Column.cpp

HEADERS +=  Cortical_Column.h	\
	    Data_Storage.h	\
	    ODE.h		\
	    Random_Stream.h	\
	    Stimulation.h	\
	    Thalamic_Column.h

QMAKE_CXXFLAGS += -std=c++11 -O3

SOURCES -= TC_mex.cpp
