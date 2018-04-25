GCC=g++
GFLAGS=-lgmp -std=c++11
all: compile

compile: main.cpp
	$(GCC) $(GFLAGS) -o kry main.cpp

clean:
	rm kry