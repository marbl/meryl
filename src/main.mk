
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.

ifeq "$(strip ${DESTDIR})" ""
  DESTDIR      := 
endif

ifeq "$(strip ${PREFIX})" ""
  ifeq "$(strip ${DESTDIR})" ""
    PREFIX     := $(realpath ..)
  else
    PREFIX     := /meryl
  endif
endif

ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := $(DESTDIR)$(PREFIX)/$(OSTYPE)-$(MACHINETYPE)/obj
endif

ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := $(DESTDIR)$(PREFIX)/$(OSTYPE)-$(MACHINETYPE)
endif

TARGET       := libmeryl.a

SOURCES      := AS_global.C \
                \
                utility/mt19937ar.C \
                utility/AS_UTL_stackTrace.C \
                utility/AS_UTL_alloc.C \
                utility/AS_UTL_fileIO.C \
                utility/AS_UTL_fasta.C \
                utility/readBuffer.C \
                utility/memoryMappedFile.C \
                \
                libsequence.C \
                libkmer.C \
                libkmer-reader.C \
                libkmer-writer.C \
                libkmer-statistics.C \
                libkmer-exact.C \
                libbits.C

ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += utility/libbacktrace/atomic.c \
                utility/libbacktrace/backtrace.c \
                utility/libbacktrace/dwarf.c \
                utility/libbacktrace/elf.c \
                utility/libbacktrace/fileline.c \
                utility/libbacktrace/mmap.c \
                utility/libbacktrace/mmapio.c \
                utility/libbacktrace/posix.c \
                utility/libbacktrace/print.c \
                utility/libbacktrace/simple.c \
                utility/libbacktrace/sort.c \
                utility/libbacktrace/state.c \
                utility/libbacktrace/unknown.c
endif

SRC_INCDIRS  := . \
                utility

SUBMAKEFILES := meryl.mk \
                lookup.mk \
                elias-fano.mk \
                sequence.mk \
                libbitsTest.mk
