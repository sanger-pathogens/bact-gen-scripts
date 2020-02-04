CC = cc
CFLAGS = -O4

LIB = -lm

OBJS=	\
		get_excreads.o \
                fast.o \

SRCS= $(OBJS:.o=.c)

INCS=   \
                fasta.h \

PROGRAM = get_excreads

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
