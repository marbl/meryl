
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "meryl.H"
#include <cmath>



bool
merylOpCompute::nextMer(void) {
  bool   isEmpty = false;

 nextMerAgain:

  //  Get some logging out of the way.

  if (verbosity.showEverything() == true) {
    char  kmerString[256];

    fprintf(stderr, "\n");
    fprintf(stderr, "merylOp::nextMer()-- STARTING for operation #%u\n", _ot->_ident);

    for (uint32 ii=0; ii<_inputs.size(); ii++)
      if (_inputs[ii]->_valid == true)
        fprintf(stderr, "merylOp::nextMer()--   CURRENT STATE: input %s kmer %s count %u\n",
                _inputs[ii]->name(),
                _inputs[ii]->_kmer.toString(kmerString),
                _inputs[ii]->_kmer._val);
      else
        fprintf(stderr, "merylOp::nextMer()--   CURRENT STATE: input %s <empty>\n",
                _inputs[ii]->name());
  }

  //  Grab the next mer for every input that was active in the last
  //  iteration.  (on the first call, all inputs were 'active' last time)

  //if ((verbosity.showDetails()    == true) &&
  //    (verbosity.showEverything() == false)) {
  //  fprintf(stderr, "\n");
  //  fprintf(stderr, "merylOp::nextMer()-- STARTING for operation #%u\n", _ot->_ident);
  //}

  for (uint32 ii=0; ii<_actLen; ii++) {
    if (verbosity.showDetails() == true)
      fprintf(stderr, "merylOp::nextMer()-- CALL NEXTMER on input actIdx " F_U32 "\n", _actIdx[ii]);
    _inputs[_actIdx[ii]]->nextMer();
  }

  //  Find the smallest kmer in any input, and remember the values and labels
  //  of the kmer in each input file.

  findOutputKmer();

  //  If no active kmers, we're done.  There are no kmers left in any input.

  if (_actLen == 0) {
    //if (verbosity.showDetails() == true) {
    //  fprintf(stderr, "merylOp::nextMer()-- No inputs found, all done here.\n");
    //}
    return(false);
  }

  //  The findOutputKmer() functions will decide, based on present-in, if we
  //  should skip this kmer.  If we're told to skip it, do it right now;
  //  there is no need to do any processing on values or labels.

  if (verbosity.showDetails() == true) {
    char  kmerString[256];
    fprintf(stderr, "merylOp::nextMer()-- op #%u activeLen " F_U32 " kmer %s\n", _ot->_ident, _actLen, _kmer.toString(kmerString));
  }

  //  Figure out the value/label of the output kmer.

  findOutputValue();
  findOutputLabel();

  //  Decide if we really want to output.

  if (isKmerFilteredOut() == true) {
    if (verbosity.showDetails() == true) {
      char  kmerString[256];
      fprintf(stderr, "merylOp::nextMer()-- FILTERED for operation #%u with kmer %s count %u label %016lx%s\n",
              _ot->_ident, _kmer.toString(kmerString), _kmer._val, _kmer._lab, ((_writerSlice != nullptr) && (_kmer._val != 0)) ? " OUTPUT" : "");
      fprintf(stderr, "\n");
    }
    goto nextMerAgain;
  }

  //  Output the fully processed kmer to whatever outputs exist.

  if (verbosity.showDetails() == true) {
    char  kmerString[256];
    fprintf(stderr, "merylOp::nextMer()-- FINISHED for operation #%u with kmer %s count %u label %016lx%s\n",
            _ot->_ident, _kmer.toString(kmerString), _kmer._val, _kmer._lab, ((_writerSlice != nullptr) && (_kmer._val != 0)) ? " OUTPUT" : "");
    fprintf(stderr, "\n");
  }

  outputKmer();
  printKmer();

  if (_stats)
    _stats->addValue(_kmer._val);

  //  We have now loaded the nextMer; return and let our clients query us.

  return(true);
}
