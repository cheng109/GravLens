
#include=-I/Users/juncheng/work/phosim_core/source/cfitsio/include
#lib=-L/Users/juncheng/work/phosim_core/source/cfitsio/lib


INC=-I/usr/local/Cellar/boost/1.58.0/include -I/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/include -I/usr/local/include
LIB=-L/usr/local/Cellar/boost/1.58.0/lib -L/usr/local/lib -L/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/lib
CC=g++
CFLAGS=-Wall -O2
LDFLAGS=-lboost_iostreams -lboost_system -lcfitsio -larmadillo

all:
	g++ $(INC) $(LIB) $(CFLAG) $(LDFLAGS) commons.cpp main.cpp Image.cpp Model.cpp Model_test.cpp -o junGL 
	./junGL

#clean: #
#	rm *.o junGL

