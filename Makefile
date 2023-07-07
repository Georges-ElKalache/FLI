CC = gcc
CFLAGS = -Wall

all: my_program

my_program: basin.o FLI.o
	$(CC) $(CFLAGS) basin.o  FLI.o -o my_program -lgsl -lgslcblas -lm

basin.o: basin.c
	$(CC) $(CFLAGS) -c basin.c 

FLI.o: FLI.c
	$(CC) $(CFLAGS) -c FLI.c 
