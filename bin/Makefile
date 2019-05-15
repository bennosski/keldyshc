CC = icc

CCOPTS = -O3 -restrict
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

HDFROOT = /share/software/user/open/hdf5/1.10.2
LIBS =  -L$(HDFROOT)/lib -lsz -lz -lm $(HDFROOT)/lib/libhdf5.a
INCLUDE = -I$(HDFROOT)/include

MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

SRCFILES = matsubara.o integration.o params.o constants.o save.o volterra.o linalg.o
TESTFILES = test_linalg.o

IDIR = ../include
_DEPS = constants.h integration.h params.h volterra.h matsubara.h save.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


all: ${SRCFILES} main.o
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -o keldysh $^  ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

stupid: stupid.o dumb.o idiot.o
	icc -o stupid $^

stupid.o: ../test/stupid.c
	icc -I../include -c -o $@ $<

dumb.o: ../test/dumb.c
	icc -I../include -c -o $@ $<

idiot.o: ../src/idiot.c
	icc -I../include -c -o $@ $<

test: ${TESTFILES} ${SRCFILES} main_test.o
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -o ../test/keldysh_test $^  ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

%.o: ../src/%.c
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -I../include -c -o $@ $< ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

main_test.o: ../test/main_test.c
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -I../include -c -o $@ $< ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

test_linalg.o: ../test/test_linalg.c
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -I../include -c -o $@ $< ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

clean:
	rm *.o

