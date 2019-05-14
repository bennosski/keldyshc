CC = icc

CCOPTS = -O2 -restrict
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

HDFROOT = /share/software/user/open/hdf5/1.10.2
LIBS =  -L$(HDFROOT)/lib -lsz -lz -lm $(HDFROOT)/lib/libhdf5.a
INCLUDE = -I$(HDFROOT)/include

MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

SRCFILES = main.o matsubara.o integration.o params.o constants.o save.o volterra.o

all: ${SRCFILES}
	${CC} ${CCOPTS} $(MPI_COMPILE_FLAGS) -o keldysh $? ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

test: ${SRCFILES}
	${CC} $(MPI_COMPILE_FLAGS) -o ../test/keldysh_test $? ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

debug: ${SRCFILES}
	${CC} -g -O0 -debug all -traceback -o keldysh_debug $? ${LIBRARIES} $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)

%.o: ../src/%.c
	${CC} $(MPI_COMPILE_FLAGS) -c -I../include ${CCOPTS} ${MKL} $< $(LIBS) $(INCLUDE) $(MPI_LINK_FLAGS)


clean:
	rm *.o
