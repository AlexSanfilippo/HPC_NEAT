#Alexander Sanfilippo
#created: 11 June, 2022
#last update: 27 August, 2022
#NEAT Project

#CC = g++
CC = mpicxx
#CFLAGS = -std=c++11 -O3 -Wall -Wreturn-type -g
CFLAGS = -std=c++11 -Wall -Wreturn-type -g
OBJS = testingGrounds.o 
EXECS = testingGrounds 
all:${EXECS}

testingGrounds: testingGrounds.o
	${CC} ${CFLAGS} testingGrounds.o -o testingGrounds
testingGrounds.o: testingGrounds.cc gene.h Genome.h Species.h NOV.h
	${CC} ${CFLAGS} -c testingGrounds.cc

#Try to compile the simplest program to call python from c++ file:FAILS
#cpp_to_python: cpp_to_python.o
#	${CC} ${CFLAGS} cpp_to_python.o -o cpp_to_python -I/home/users/mschpc/2021/sanfilia/miniconda/envs/summer_project/include/python3.9
#cpp_to_python.o: cpp_to_python.cc Python.h
#	${CC} ${CFLAGS} -I/home/users/mschpc/2021/sanfilia/miniconda3/envs/summer_project/include/python3.9 -c cpp_to_python.cc -l python3.9

.PHONY: clean
clean:
	rm ${OBJS} ${EXECS}


