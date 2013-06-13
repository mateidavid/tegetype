-include Makefile.local

SHELL := /bin/bash
LD := ${CXX}

ifndef BUILD_TYPE
  BUILD_TYPE := testing
endif

ifndef CXXFLAGS
  ifeq (${BUILD_TYPE}, production)
    CXXFLAGS := -g -O2 -DNDEBUG
  else
    ifeq (${BUILD_TYPE}, testing)
      CXXFLAGS := -g -O2
    else
      CXXFLAGS := -g
    endif
  endif
endif

CXXFLAGS += -std=c++0x -mmmx -msse -msse2 -fopenmp
CXXFLAGS += -Wall -Wextra -Wno-unused-local-typedefs

# bamtools
CXXFLAGS += -I bamtools/include
LDFLAGS += -L bamtools/lib -Wl,--rpath=bamtools/lib

#GIT_VERSION=$(shell ./get_git_version)
#override CXXFLAGS+=-DGIT_VERSION=$(GIT_VERSION)


OBJS := get-frag-rates.o get-ref-gc.o globals.o Pairing.o Fasta.o
DEPS := $(OBJS:.o=.d)
TGTS := get-frag-rates get-ref-gc

.PHONY: all clean

all: ${TGTS}

clean:
	rm -f ${OBJS} ${DEPS} ${TGTS}

-include ${DEPS}

%.o : %.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $+
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -MM -MP $< >$*.d

get-frag-rates: get-frag-rates.o globals.o Pairing.o Fasta.o
	${LD} -o $@ $+ ${LDFLAGS} -lbamtools -lboost_iostreams

get-ref-gc: get-ref-gc.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams -lboost_regex
