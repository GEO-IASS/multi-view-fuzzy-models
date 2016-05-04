CC=gcc
CFLAGS=-Wall -c
EXEDIR=.

all: $(EXEDIR)/main

$(EXEDIR)/main: util.o main.o
	$(CC) util.o main.o -o $@ -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

util.o: util.c
	$(CC) $(CFLAGS) util.c

clean:
	rm -rf *o $(EXEDIR)/main
