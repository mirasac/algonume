MAIN=main
SURNAME=casari
LIB=../../mclib/mclib

prepare: ${MAIN}.cpp
	g++ -E -DPROD ${MAIN}.cpp -o ${SURNAME}_${MAIN}.cpp

compile: main.o lib.o
	g++ main.o lib.o -o ${MAIN}.out

main.o: ${MAIN}.cpp ${LIB}.h
	g++ -c ${MAIN}.cpp -o main.o

lib.o: ${LIB}.cpp ${LIB}.h
	g++ -c ${LIB}.cpp -o lib.o

clean:
	rm *.o *.out *.dat
