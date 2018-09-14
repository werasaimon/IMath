CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -std=c++11
CC=g++

OBJS =		IMath_make.o

TARGET =	IMath_bin

#----------------------- core main -------------------------#
LDFLAGS =
SOURCES   = src/main.cpp


OBJS=$(SOURCES:%.cpp=%.o)
EXECUTABLE=hello

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

#clean:
#rm -f $(OBJS) $(TARGET)
