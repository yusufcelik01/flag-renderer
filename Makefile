CC=gcc
CXX=g++
GLFLAGS= -g -lglfw3 -lpthread -lX11 -ldl -lXrandr -lGLEW -lGL -DGL_SILENCE_DEPRECATION -DGLM_ENABLE_EXPERIMENTAL -I. -msse2
CXXFLAGS=$(GLFLAGS)

flag: flag.cpp helpers.o
	g++ helpers.o flag.cpp -o flag $(CXXFLAGS)

bezier: bezier.cpp
	g++ bezier.cpp -o bezier $(CXXFLAGS)

all:
	g++ main.cpp -o main $(GLFLAGS)

%: %.cpp
	$(CXX) $@.cpp -o $@ $(GLFLAGS)



%.o: %.cpp %.hpp
	$(CXX) -g $< -c

.PHONY:clean
clean:
	-rm flag
	-rm *.o
