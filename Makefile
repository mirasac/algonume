MAIN=main
SURNAME=casari
LIB=./mclib/mclib

dev: main lib
	g++ ${CPPFLAGS} -o ${MAIN}.out ${MAIN}.o ${LIB}.o

prod: ${MAIN}.cpp
	g++ ${CPPFLAGS} -DPROD -o ${MAIN}.cpp ${MAIN}.cpp

main: ${MAIN}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${MAIN}.cpp -o ${MAIN}.o

lib: ${LIB}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${LIB}.cpp -o ${LIB}.o

clean:
	rm *.o *.out *.dat
