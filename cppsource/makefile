#include=-I/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/include
#lib=-L/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/lib

include=-I/Users/juncheng/work/phosim_core/source/cfitsio/include
lib=-L/Users/juncheng/work/phosim_core/source/cfitsio/lib

include2=-I/usr/local/include
lib2=-L/usr/local/lib


all:
	g++ $(include) $(lib) $(include2) $(lib2)  commons.cpp main.cpp Image.cpp Model.cpp Model_test.cpp -o junGL -O2 -lcfitsio -larmadillo
	./junGL

	

