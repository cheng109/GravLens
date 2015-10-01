include=-I/opt/local/include
lib=-L/opt/local/lib
all:
	g++ -o main.o main.cpp 
	gcc -o commons.o commons.c
	g++ -o junGL main.o commons.o
	./junGL
