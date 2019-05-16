CC = icc

CCOPTS = -O3 -restrict
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

HDFROOT = /share/software/user/open/hdf5/1.10.2
LIBS =  -L$(HDFROOT)/lib -lsz -lz -lm $(HDFROOT)/lib/libhdf5.a
INCLUDE = -I$(HDFROOT)/include

MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

SRCFILES = matsubara.o integration.o params.o constants.o save.o volterra.o linalg.o util.o
TESTFILES = test_linalg.o test_MxM.o test_rcorr.o

all: ${SRCFILES} main.o
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -o keldysh $^  ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

test: ${TESTFILES} ${SRCFILES} main_test.o
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -o keldysh_test $^  ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

%.o: ../src/%.c
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -I../include -c -o $@ $< ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

%.o: ../test/%.c
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -I../include -c -o $@ $< ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

clean:
	rm *.o
