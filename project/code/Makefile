MAIN=main
LIB=../../mclib/mclib
OBJECTS=${MAIN}.o ${LIB}.o

dev: ${OBJECTS}
	g++ ${CPPFLAGS} -o ${MAIN}.out ${OBJECTS}

prod: ${OBJECTS}
	g++ ${CPPFLAGS} -Wall -pedantic -pedantic-errors -Wextra -DPROD -o ${MAIN}.out ${OBJECTS}

${MAIN}.o: ${MAIN}.cpp constants.h ${LIB}.o
	g++ ${CPPFLAGS} -c ${MAIN}.cpp -o ${MAIN}.o

${LIB}.o: ${LIB}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${LIB}.cpp -o ${LIB}.o

radiation.o: radiation.cpp radiation.h constants.h ${LIB}.o configuration.h
	g++ ${CPPFLAGS} -c radiation.cpp -o radiation.o

utilities.o: utilities.cpp utilities.h constants.h radiation.o configuration.h
	g++ ${CPPFLAGS} -c utilities.cpp -o utilities.o

convection.o: convection.cpp convection.h
	g++ ${CPPFLAGS} -c convection.cpp -o convection.o

clean:
	rm -f *.o
	rm -f *.out
