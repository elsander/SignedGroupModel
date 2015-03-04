COMP = gcc
LIB = -lgsl -lgslcblas
OPT = -Wall -O3 -DHAVE_INLINE
OBJECTS = main.o partition.o SearchAlgs.o

FindGroupsGibbs: ${OBJECTS}
	${COMP} ${OBJECTS} -o FindGroupsGibbs ${OPT} ${LIB}

partition.o: partition.c Common.h
	${COMP} partition.c -c ${OPT}

SearchAlgs.o: SearchAlgs.c partition.h Common.h
	${COMP} SearchAlgs.c -c ${OPT}

main.o: main.c Common.h SearchAlgs.h partition.h
	${COMP} main.c -c ${OPT}

clean: 
	rm -rf *.o *~ *.log FindGroupsGibbs
