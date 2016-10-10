TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = release_binary

SOURCES +=  Cortical_Column.cpp \
			TC.cpp				\
			TC_mex.cpp			\
			Thalamic_Column.cpp

HEADERS +=  Cortical_Column.h	\
			Data_Storage.h		\
			ODE.h				\
			Random_Stream.h		\
			Stimulation.h		\
			Thalamic_Column.h

SOURCES -= TC_mex.cpp

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE *= -O3

