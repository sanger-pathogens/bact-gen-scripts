CC = cc
CFLAGS = -O2

LIB = -lm

OBJS=	\
		fastq2fasta.o \
                fast.o \

SRCS= $(OBJS:.o=.c)

PROGRAM = fastq2fasta 

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
