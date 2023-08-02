MAIN=main
LIB=../../mclib/mclib
OBJECTS=${MAIN}.o ${LIB}.o functions.o radiation.o

dev: ${OBJECTS}
	g++ ${CPPFLAGS} -o ${MAIN}.out ${OBJECTS}

prod: ${OBJECTS}
	g++ ${CPPFLAGS} -DPROD -o ${MAIN}.out ${OBJECTS}

${MAIN}.o: ${MAIN}.cpp ${LIB}.h constants.h
	g++ ${CPPFLAGS} -c ${MAIN}.cpp -o ${MAIN}.o

${LIB}.o: ${LIB}.cpp ${LIB}.h
	g++ ${CPPFLAGS} -c ${LIB}.cpp -o ${LIB}.o

functions.o: functions.cpp functions.h constants.h ${LIB}.o
	g++ ${CPPFLAGS} -c functions.cpp -o functions.o

radiation.o: radiation.cpp radiation.h constants.h functions.o
	g++ ${CPPFLAGS} -c radiation.cpp -o radiation.o

clean:
	rm -f *.o
	rm -f *.out
