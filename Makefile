CC=gcc
CXX=g++
GLFLAGS= -g -lglfw -lpthread -lX11 -ldl -lXrandr -lGLEW -lGL -DGL_SILENCE_DEPRECATION -DGLM_ENABLE_EXPERIMENTAL -I. -msse2
CXXFLAGS=$(GLFLAGS)

main: main.cpp helpers.o
	g++ helpers.o main.cpp -o main $(CXXFLAGS)


all:
	g++ main.cpp -o main $(GLFLAGS)

%: %.cpp
	$(CXX) $@.cpp -o $@ $(GLFLAGS)


hw:
	tar -cvzf hw1.tar.gz *.cpp *.hpp *.h Makefile *.glsl metu_flag.jpg

%.o: %.cpp %.hpp
	$(CXX) -g $< -c

.PHONY:clean
clean:
	-rm main
	-rm *.o
