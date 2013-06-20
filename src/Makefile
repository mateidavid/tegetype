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

CXXFLAGS += -std=c++0x -mmmx -msse -msse2 -fopenmp -Wall -Wextra
LDFLAGS += -fopenmp
-include Makefile.local


# bamtools
BAMTOOLS_DIR := $(shell readlink -m bamtools)
BAMTOOLS_CFLAGS := -I ${BAMTOOLS_DIR}/include
BAMTOOLS_LDFLAGS := -L ${BAMTOOLS_DIR}/lib -Wl,--rpath=${BAMTOOLS_DIR}/lib
CXXFLAGS += ${BAMTOOLS_CFLAGS}
LDFLAGS += ${BAMTOOLS_LDFLAGS}


#GIT_VERSION=$(shell ./get_git_version)
#override CXXFLAGS+=-DGIT_VERSION=$(GIT_VERSION)


OBJS := DNASequence.o Read.o Cigar.o Mapping.o Pairing.o Fasta.o \
	Clone.o CloneGen.o SamMapping.o SamMappingSetGen.o \
	globals.o common.o deep_size.o util.o \
	get-frag-gc.o get-ref-gc.o get-te-evidence.o combine-evidence.o \
	add-extra-sam-flags.o filter-mappings.o sam-to-fq.o \
	zc.o tee-p.o printab.o

DEPS := $(OBJS:.o=.d)

TGTS := get-frag-gc get-ref-gc get-te-evidence combine-evidence \
	add-extra-sam-flags filter-mappings sam-to-fq \
	zc tee-p printab

BIN_PATH := ../bin
TGTS_W_PATH := $(foreach tgt,${TGTS},${BIN_PATH}/${tgt})


.PHONY: all clean

all: ${TGTS_W_PATH}

clean:
	rm -f ${OBJS} ${DEPS} ${TGTS_W_PATH}

-include ${DEPS}

%.o : %.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -MMD -MP -o $@ -c $<

%.o : %.c
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -MMD -MP -o $@ -c $<


${BIN_PATH}/get-frag-gc: get-frag-gc.o globals.o Pairing.o Fasta.o
	${LD} -o $@ $+ ${LDFLAGS} -lbamtools -lboost_iostreams

${BIN_PATH}/get-ref-gc: get-ref-gc.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams -lboost_regex

${BIN_PATH}/get-te-evidence: get-te-evidence.o globals.o Clone.o CloneGen.o Mapping.o \
	SamMapping.o SamMappingSetGen.o Pairing.o common.o Read.o Cigar.o \
	DNASequence.o deep_size.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams

${BIN_PATH}/combine-evidence: combine-evidence.o globals.o Pairing.o Fasta.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams

${BIN_PATH}/add-extra-sam-flags: add-extra-sam-flags.o globals.o util.o deep_size.o \
	Pairing.o DNASequence.o Read.o Mapping.o CloneGen.o \
	SamMapping.o SamMappingSetGen.o common.o Cigar.o Clone.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams

${BIN_PATH}/filter-mappings: filter-mappings.o globals.o util.o deep_size.o Pairing.o \
	DNASequence.o Read.o CloneGen.o SamMapping.o SamMappingSetGen.o common.o \
	Clone.o Mapping.o Cigar.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams

${BIN_PATH}/sam-to-fq: sam-to-fq.o globals.o util.o deep_size.o Pairing.o \
	DNASequence.o Read.o CloneGen.o SamMapping.o SamMappingSetGen.o common.o \
	Clone.o Mapping.o Cigar.o
	${LD} -o $@ $+ ${LDFLAGS} -lboost_iostreams

${BIN_PATH}/zc: zc.o
	${LD} -o $@ $+ ${LDFLAGS} -lz

${BIN_PATH}/tee-p: tee-p.o
	${LD} -o $@ $+ ${LDFLAGS}

${BIN_PATH}/printab: printab.o
	${LD} -o $@ $+ ${LDFLAGS}
