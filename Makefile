SRCS := $(wildcard src/*.c)

PROGS = $(patsubst %.c,%,$(SRCS))

CFLAGS = -Ofast -g -I$(PWD)/submodules/htslib -fvisibility=hidden -fpic -c -Wall

LDFLAGS = -L$(PWD)/submodules/htslib -Wl,-rpath=$(PWD)/submodules/htslib -lhts

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
