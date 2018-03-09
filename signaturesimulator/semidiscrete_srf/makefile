CC = gcc
FF=gfortran

CFLAGS = -Wall -O3
LIBS = -lm 
INCLS = -I./

OUT = semiD
FOBJECTS = semiDbrf.o semiDtools.o prospect.redux.o semiDbrfe.o
COBJECTS = semiD.o prospect_interface.o price_soil.o libmatrixio.o


$(OUT):		${FOBJECTS} $(COBJECTS)
		${FF} ${CFLAGS} ${FOBJECTS} ${COBJECTS} -o ${OUT} 



#OUT = semiD_nohotspot
#FOBJECTS = semiDbrf.o semiDtools_nohotspot.o prospect.redux.o semiDbrfe.o

#$(OUT):		${FOBJECTS} $(COBJECTS)
#		${FF} ${CFLAGS} ${FOBJECTS} ${COBJECTS} -o ${OUT} 

.c.o: $<
		$(CC) ${INCLS} $(CFLAGS) -c $<

.f.o: $<
		$(FF) $(CFLAGS) -c $<

clean:
		\rm -f *.o *~ *%

#install:
#		cp prospect $(HOME)/bin/$(ARCH)/prospect.redux

