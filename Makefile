TARGET=simulaGra
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:%.cpp=%.o)
FLAGS= -fopenmp -I/usr/include/freetype2 -std=c++11 `pkg-config --cflags opencv`
CCLIBS= -L/usr/local/share/OpenCV/3rdparty/lib -lGL -lGLU -lglfw -lftgl `pkg-config --libs opencv`

all: version.h $(OBJS)
	$(CXX) $(FLAGS) -o $(TARGET) $(OBJS) $(CCLIBS)

version.h:
	sh version.sh

.cpp.o:
	$(CXX) $(FLAGS) -c -o $*.o $<

clean:
	$(RM) $(TARGET) *.o *~ core* version.h

