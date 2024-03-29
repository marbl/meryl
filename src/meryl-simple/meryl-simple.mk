TARGET   := meryl-simple
SOURCES  := meryl-simple.C

SRC_INCDIRS := .

#  If we're part of Canu, build with canu support and use Canu's copy of
#  meryl-utility.  Otherwise, don't.
ifneq ($(wildcard stores/sqStore.H), )
  SRC_CXXFLAGS := -DCANU
  SRC_INCDIRS  := ../../../utility/src ../../../stores

#  If we're part of something else, include the something else's
#  utility directory.
else ifneq ($(wildcard meryl/src/meryl/meryl.C), )
  SRC_INCDIRS  := ../../../utility/src

#  Otherwise, we're building directly in the meryl repo.
else
  SRC_INCDIRS  := ../utility/src

endif


TGT_LDFLAGS  := -L${TARGET_DIR}/lib
TGT_LDLIBS   := -l${MODULE}
TGT_PREREQS  := lib${MODULE}.a
