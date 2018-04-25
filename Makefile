GCC=g++
GFLAGS=-lgmp -std=c++11
all: compile

compile: main.cpp
	$(GCC) $(GFLAGS) -o rsa.o main.cpp

clean:
	rm *.o