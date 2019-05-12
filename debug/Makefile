CC = icc

CCOPTS = -O2 -restrict
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

all:
	${CC} ${CCOPTS} -o a.out main.c ${LIBRARIES}

debug:
	${CC} -g -o a.out main.c ${LIBRARIES}

test:
	${CC} -g -o a.out test.c ${LIBRARIES}
