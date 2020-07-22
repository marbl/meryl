TARGET   := meryl-check
SOURCES  := meryl-check.C \

SRC_INCDIRS  := . ../utility/src/utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lmeryl
TGT_PREREQS := libmeryl.a

SUBMAKEFILES :=
