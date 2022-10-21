MAIN=main
LIB=../mclib/mclib

main.out: main.o lib.o
	g++ main.o lib.o -o main.out

main.o: ${MAIN}.cpp ${LIB}.h
	g++ -c ${MAIN}.cpp -o main.o

lib.o: ${LIB}.cpp ${LIB}.h
	g++ -c ${LIB}.cpp -o lib.o

clean:
	rm *.o *.out

