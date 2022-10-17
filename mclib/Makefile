MAIN=main
LIB=../mclib/mclib

main: main.o lib.o
	g++ main.o lib.o -o main

main.o: ${MAIN}.cpp ${LIB}.h
	g++ -c ${MAIN}.cpp -o main.o

lib.o: ${LIB}.cpp ${LIB}.h
	g++ -c ${LIB}.cpp -o lib.o

clean:
	rm *.o

