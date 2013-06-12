SHELL=/bin/bash

ifndef BUILD_TYPE
  BUILD_TYPE := testing
endif

ifndef CXXFLAGS
  ifeq ($(BUILD_TYPE), production)
    CXXFLAGS := -g -O2 -DNDEBUG
  else
    ifeq ($(BUILD_TYPE), testing)
      CXXFLAGS := -g -O2
    else
      CXXFLAGS := -g
    endif
  endif
endif

CXXFLAGS += -I. -mmmx -msse -msse2 -fopenmp -Wall -Wextra -Wno-deprecated

-include Makefile.local

#GIT_VERSION=$(shell ./get_git_version)
#override CXXFLAGS+=-DGIT_VERSION=$(GIT_VERSION)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  LDFLAGS := -lm -lz
else
  LDFLAGS := -lm -lz -lrt
endif

ifdef BOOST
  BOOST_INCLUDE := ${BOOST}/include
  BOOST_LIB := ${BOOST}/lib
endif

ifdef BOOST_INCLUDE
  CXXFLAGS += -I ${BOOST_INCLUDE}
endif

ifdef BOOST_LIB
  LDFLAGS += -L ${BOOST_LIB} -Wl,--rpath=${BOOST_LIB}
endif


all: get-frag-rates

Fasta.o: Fasta.?pp DNASequence.hpp globals.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
