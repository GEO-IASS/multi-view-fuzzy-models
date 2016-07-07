CC=gcc
CFLAGS=-Wall -c
EXEDIR=.

all: $(EXEDIR)/main

$(EXEDIR)/main: stex.o matrix.o util.o main.o
	$(CC) stex.o matrix.o util.o main.o -o $@ -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

util.o: util.c
	$(CC) $(CFLAGS) util.c

matrix.o: matrix.c
	$(CC) $(CFLAGS) matrix.c

stex.o: stex.c
	$(CC) $(CFLAGS) stex.c

clean:
	rm -rf *o $(EXEDIR)/main
