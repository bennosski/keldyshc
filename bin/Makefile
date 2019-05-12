CC = icc

CCOPTS = -O2 -restrict
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

HDFROOT = /share/software/user/open/hdf5/1.10.2
LIBS =  -L$(HDFROOT)/lib -lsz -lz -lm $(HDFROOT)/lib/libhdf5.a
INCLUDE = -I$(HDFROOT)/include

SRCFILES = main.o matsubara.o integration.o params.o constants.o save.o

all: ${SRCFILES}
	${CC} ${CCOPTS} -o keldysh $? ${LIBRARIES} $(LIBS) $(INCLUDE)
	make clean

#%.o: ../src/%.c
#	${CC} -c -I../include ${CCOPTS} ${MKL} ../include/params.h $<

%.o: ../src/%.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $< $(LIBS) $(INCLUDE)

debug:
	${CC} -g -o a.out main.c ${LIBRARIES} $(LIBS) $(INCLUDE)

test:
	${CC} -g -o a.out test.c ${LIBRARIES} $(LIBS) $(INCLUDE)


clean:
	rm *.o
