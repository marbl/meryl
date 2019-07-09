
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)
endif

TARGET   := meryl-lookup
SOURCES  := meryl-lookup.C \

SRC_INCDIRS  := . .. ../utility ../stores

#  If we're part of Canu, build with canu support.
#  Otherwise, don't.

ifneq ($(wildcard stores/sqStore.H), )

SRC_CXXFLAGS := -DCANU

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

else

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lmeryl
TGT_PREREQS := libmeryl.a

endif

SUBMAKEFILES :=
