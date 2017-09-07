#CXXFLAGS=-Wall -O0 -g -Wno-deprecated -lgomp -fopenmp
CXXFLAGS=-Wall -std=c++0x -O0 -g -Wno-deprecated
LIBS=-L/usr/include/lib -lboost_program_options
LIBS2=-I /usr/include
BIN=mp_edp

 
all: $(BIN)

mp_edp: mp_edp.cpp mes.o proba.hpp
	g++ ${CXXFLAGS} ${LIBS2} mp_edp.cpp ${LIBS} mes.o -o mp_edp

clean:
	rm -f $(BIN) *.o
mes.o: mes.cpp mes.hpp proba.hpp
	g++ ${CXXFLAGS} -c ${LIBS2} mes.cpp -o mes.o


