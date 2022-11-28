SRCS = $(wildcard *.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g -Ihtslib-1.16 -fvisibility=hidden -fpic -c -Wall

INSTALLDIR ?= /usr/local/bin

CC=gcc

all: $(PROGS)

%: %.c

install: $(PROGS)
	cp $(PROGS) $(INSTALLDIR)

clean: 
	rm -f $(PROGS) *.o
