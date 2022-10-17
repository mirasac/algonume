main : main.o lib.o
	g++ main.o lib.o -o main

main.o : main.cpp lib.h
	g++ -c main.cpp -o main.o

lib.o : lib.cpp lib.h
	g++ -c lib.cpp -o lib.o

clean :
	rm *.o

