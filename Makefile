SRCS := $(wildcard src/*.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g -Isubmodules/htslib -fvisibility=hidden -fpic -c -Wall

LDFLAGS = -Lsubmodules/htslib -Wl,-rpath=$(PWD)/submodules/htslib

LDLIBS = -lhts

INSTALLDIR ?= bin

CC=gcc

all: $(PROGS)

%: %.c

install: $(PROGS)
	cp -p $(PROGS) $(INSTALLDIR)

htslib:
	cd submodules/htslib && autoreconf -i && ./configure && make

clean: 
	rm -f $(PROGS) *.o
