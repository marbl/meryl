TARGET   := meryl
SOURCES  := meryl.C \
            merylCommandBuilder-isAlias.C \
            merylCommandBuilder-isFilter.C \
            merylCommandBuilder-isOption.C \
            merylCommandBuilder-labelSelect.C \
            merylCommandBuilder-processWord.C \
            merylCommandBuilder-valueSelect.C \
            merylCommandBuilder.C \
            merylCountArray.C \
            merylFilter.C \
            merylInput.C \
            merylOp-count-memorySize.C \
            merylOp-count.C \
            merylOp-countSequential.C \
            merylOp-countSimple.C \
            merylOp-countThreads.C \
            merylOp-histogram.C \
            merylOp-nextMer.C \
            merylOp.C \
            merylOpCompute.C \
            merylOpTemplate.C

SRC_INCDIRS := .

#  If we're part of Canu, build with canu support and use Canu's copy of
#  meryl-utility.  Otherwise, don't.
ifneq ($(wildcard stores/sqStore.H), )
  SRC_CXXFLAGS := -DCANU
  SRC_INCDIRS  := ../../../utility/src/utility ../../../stores

#  If we're part of something else, include the something else's
#  utility directory.
else ifneq ($(wildcard meryl/src/meryl/meryl.C), )
  SRC_INCDIRS  := ../../../utility/src/utility

#  Otherwise, we're building directly in the meryl repo.
else
  SRC_INCDIRS  := ../utility/src/utility

endif


TGT_LDFLAGS  := -L${TARGET_DIR}/lib
TGT_LDLIBS   := -l${MODULE}
TGT_PREREQS  := lib${MODULE}.a
