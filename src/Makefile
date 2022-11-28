SRCS = $(wildcard *.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g -lhts

INSTALLDIR ?= /usr/local/bin

CC=gcc

all: $(PROGS)

%: %.c

install: $(PROGS)
	cp $(PROGS) $(INSTALLDIR)

clean: 
	rm -f $(PROGS) *.o
