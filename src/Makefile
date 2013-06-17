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
CXXFLAGS += -Wall -Wextra
-include Makefile.local


# bamtools
BAMTOOLS_DIR := $(shell readlink -m bamtools)
BAMTOOLS_CFLAGS := -I ${BAMTOOLS_DIR}/include
BAMTOOLS_LDFLAGS := -L ${BAMTOOLS_DIR}/lib -Wl,--rpath=${BAMTOOLS_DIR}/lib
CXXFLAGS += ${BAMTOOLS_CFLAGS}
LDFLAGS += ${BAMTOOLS_LDFLAGS}


#GIT_VERSION=$(shell ./get_git_version)
#override CXXFLAGS+=-DGIT_VERSION=$(GIT_VERSION)


OBJS := get-frag-gc.o get-ref-gc.o globals.o Pairing.o Fasta.o
DEPS := $(OBJS:.o=.d)
TGTS := get-frag-gc get-ref-gc
BIN_PATH := ../bin
TGTS_W_PATH := $(foreach tgt,${TGTS},${BIN_PATH}/${tgt})

.PHONY: all clean

all: ${TGTS_W_PATH}

clean:
	rm -f ${OBJS} ${DEPS} ${TGTS_W_PATH}

-include ${DEPS}

%.o : %.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -MMD -MP -o $@ -c $<

%/get-frag-gc: get-frag-gc.o globals.o Pairing.o Fasta.o
	${LD} -o $@ $+ ${LDFLAGS} -lbamtools -lboost_iostreams

%/get-ref-gc: get-ref-gc.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams -lboost_regex
