.SUFFIXES: .F .o
include ../../Makefile.include

INCDIR = ../../include
MSGLIB = ../../lib/libMSG.a

OBJS = cg3.o diststencilhalo3D.o jacobipre.o idpre.o linezpre.o

all:	cg3

cg3:	$(OBJS)
	$(MPIF77) $(LFLAGS) -o cg3 $(OBJS) $(MSGLIB)


cg3.o:		cg.h
diststencilhalo3D.o:	cg.h
jacobipre.o:	cg.h
idpre.o:	cg.h
linezpre.o:	cg.h	


.F.o :
		$(MPIF77) $(FFLAGS) -I$(INCDIR) -c $<

clean:
		$(RM) -f $(OBJS) mpif.h

veryclean:	clean
		$(RM) -f cg3