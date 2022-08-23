#Alexander Sanfilippo
#11 June, 2022
#NEAT Project

#sequential
#CC = g++
CC = mpicxx
CFLAGS = -std=c++11 -O3 -Wall -Wreturn-type -g
OBJS = testingGrounds.o cpp_to_python.o
EXECS = testingGrounds cpp_to_python
all:${EXECS}

testingGrounds: testingGrounds.o
	${CC} ${CFLAGS} testingGrounds.o -o testingGrounds
testingGrounds.o: testingGrounds.cc gene.h Genome.h Species.h NOV.h
	${CC} ${CFLAGS} -c testingGrounds.cc

#Try to compile the simplest program to call python from c++ file:FAILS
cpp_to_python: cpp_to_python.o
	${CC} ${CFLAGS} cpp_to_python.o -o cpp_to_python -I/home/users/mschpc/2021/sanfilia/miniconda/envs/summer_project/include/python3.9
cpp_to_python.o: cpp_to_python.cc Python.h
	${CC} ${CFLAGS} -I/home/users/mschpc/2021/sanfilia/miniconda3/envs/summer_project/include/python3.9 -c cpp_to_python.cc -l python3.9
.PHONY: clean
clean:
	rm ${OBJS} ${EXECS} scripty amen cpp_to_python.o cpp_to_pythonscript.o

#reference code from CPP module for setting up makefiles
#CC = g++
#CFLAGS = -std=c++14 -fconcepts -O3 -Wall
#CFLAGS2 = -lboost_mpi -lboost_serialization -fconcepts -O3 -Wall
#OBJS = assignment5a.o assignment5b.o assignment5c.o
#EXECS = assignment5a assignment5b assignment5c
#all:${EXECS}


#assignment5a: assignment5a.o
#        ${CC} ${CFLAGS} -pthread assignment5a.o -o assignment5a
#assignment5a.o: assignment5a.cc
#        ${CC} ${CFLAGS} -c assignment5a.cc
#assignment5b: assignment5b.o
#        ${CC} ${CFLAGS} -pthread assignment5b.o -o assignment5b
#assignment5b.o: assignment5b.cc point.h ConvexHull.h
#        ${CC} ${CFLAGS} -c assignment5b.cc
#assignment5c: assignment5c.o
#        mpicxx ${CFLAGS2} assignment5c.o -o assignment5c
#assignment5c.o: assignment5c.cc point.h ConvexHull.h
#        mpicxx ${CFLAGS2} -c assignment5c.cc

