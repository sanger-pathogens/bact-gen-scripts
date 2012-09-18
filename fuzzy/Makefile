#CC = gcc
CFLAGS = -O4
#CFLAGS = -g 
#CFLAGS = -O4

LIB = -lm

OBJS=	\
		main.o \
		fast.o \
		fuzzypath.o

SRCS= $(OBJS:.o=.c)

INCS=	\
		fasta.h \
		ssaha.h

PROGRAM = fuzzypath 

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
