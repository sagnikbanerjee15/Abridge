SRCS := $(wildcard src/*.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g -fvisibility=hidden -fpic -c -Wall

LDLIBS = -lhts

INSTALLDIR ?= bin

CC=gcc

all: $(PROGS)

%: %.c

install: $(PROGS)
	cp -p $(PROGS) $(INSTALLDIR)

htslib:
	cd submodules/htslib && autoreconf -i && ./configure && make && make install

clean: 
	rm -f $(PROGS) *.o
