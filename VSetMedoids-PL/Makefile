CC=gcc
CFLAGS=-Wall -c
EXEDIR=.

all: $(EXEDIR)/main

$(EXEDIR)/main: main.o
	$(CC) main.o -o $@ -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

clean:
	rm -rf *o $(EXEDIR)/main
