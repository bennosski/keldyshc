CC = icc

CCOPTS = -O2 -restrict
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

all:
	${CC} ${CCOPTS} -o a.out main.c ${LIBRARIES}

debug:
	${CC} -g -o a.out main.c ${LIBRARIES}

linalg:
	${CC} -o linalg linalg.c ${LIBRARIES}


mklmalloc: mklmalloc.o
	${CC} -o mklmalloc mklmalloc.c ${LIBRARIES}

malloc_basic: malloc_basic.o integration.o
	${CC} -o malloc_basic malloc_basic.o integration.o ${LIBRARIES}

malloc_basic.o : malloc_basic.c
	${CC} -c malloc_basic.c ${LIBRARIES}

integration.o : integration.c integration.h
	${CC} -c integration.c ${LIBRARIES}

test:
	${CC} -g -o a.out test.c ${LIBRARIES}
