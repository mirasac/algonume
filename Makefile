MAIN=main
SURNAME=casari
LIB=./mclib/mclib

dev: main lib
	g++ ${CPPFLAGS} -o ${MAIN}.out ${MAIN}.o ${LIB}.o

prod: main lib
	g++ ${CPPFLAGS} -DPROD -o ${MAIN}.out ${MAIN}.o ${LIB}.o

main: ${MAIN}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${MAIN}.cpp -o ${MAIN}.o

lib: ${LIB}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${LIB}.cpp -o ${LIB}.o

clean:
	find . -name "*.o" -type f -delete
	find . -name "*.out" -type f -delete
	find . -name "*.dat" -type f -delete
