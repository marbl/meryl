TARGET   := meryl-import
SOURCES  := meryl-import.C \
            ../meryl/merylCountArray.C

SRC_INCDIRS  := . ../utility/src/utility ../meryl

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lmeryl
TGT_PREREQS := libmeryl.a

SUBMAKEFILES :=
