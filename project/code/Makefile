MAIN=main
LIB=../../mclib/mclib

dev: main lib
	g++ ${CPPFLAGS} -o ${MAIN}.out ${MAIN}.o ${LIB}.o

prod: main lib
	g++ ${CPPFLAGS} -DPROD -o ${MAIN}.out ${MAIN}.o ${LIB}.o

main: ${MAIN}.cpp ${MAIN}.h ${LIB}.h 
	g++ ${CPPFLAGS} -c ${MAIN}.cpp -o ${MAIN}.o

lib: ${LIB}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${LIB}.cpp -o ${LIB}.o

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.dat
