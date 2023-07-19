CC = gcc
CFLAGS = -Wall

all: my_program

my_program: FLI.o
	$(CC) $(CFLAGS)  FLI.o -o my_program -lgsl -lgslcblas -lm


FLI.o: FLI.c
	$(CC) $(CFLAGS) -c FLI.c 
