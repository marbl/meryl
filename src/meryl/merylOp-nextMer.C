
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



void
printKmer(FILE *f, kmer pk, bool printACGTorder) {
  char   outstr[256];
  char  *outptr = outstr;

  if (printACGTorder == true)              //  If requested, recompute the canonical
    pk.recanonicalizeACGTorder();          //  mer in ACGT order.  Yuck.

  pk.toString(outptr);                     //  Convert the kmer to ASCII, then
  while (*outptr)                          //  advance to the end of the string.
    outptr++;

  *outptr++ = '\t';                        //  Add the value.  There is always
  outptr = toDec(pk._val, outptr);         //  a value to add.

#warning need to print label as binary or hex, user supplied
  if (kmer::labelSize() > 0) {             //  If a label exists, add it too.
    *outptr++ = '\t';
    outptr = toBin(pk._lab, outptr, kmer::labelSize());
  }

  *outptr++ = '\n';                        //  Terminate the string and
  *outptr++ = 0;                           //  emit it.

#pragma omp critical (printLock)           //  fputs() is not thread safe and will
  fputs(outstr, f);                        //  happily intermix on e.g. Linux.
}



bool
merylOpCompute::isKmerFilteredOut(void) {
  bool   r = false;

  if (_ot->_filter.size() == 0)
    return(false);

  for (uint32 ii=0; ii<_ot->_filter.size(); ii++) {
    bool t = true;

    for (uint32 tt=0; tt<_ot->_filter[ii].size(); tt++)
      t &= _ot->_filter[ii][tt].isTrue(_kmer, _actLen, _act, _actIdx, _actRdx);

    r |= t;
  }

  return(r == false);
}



bool
merylOpCompute::nextMer(void) {
  bool   isEmpty = false;

 nextMerAgain:

  //  Grab the next mer for every input that was active in the last
  //  iteration.  (on the first call, all inputs were 'active' last time)

  for (uint32 ii=0; ii<_actLen; ii++)
    _inputs[_actIdx[ii]]->nextMer();

  //  Find the smallest kmer in any input, and remember the values and labels
  //  of the kmer in each input file.

  findOutputKmer();

  //  If no active kmers, we're done.  There are no kmers left in any input.

  if (_actLen == 0)
    return(false);

  //  Figure out the value/label of the output kmer.

  findOutputValue();
  findOutputLabel();

  //  If this kmer is filtered, go back and get another kmer from the inputs.

  if (isKmerFilteredOut() == true)
    goto nextMerAgain;

  //  Output the fully processed kmer to whatever outputs exist.

  if (_writerSlice != nullptr)
    _writerSlice->addMer(_kmer);

  if (_printer != nullptr)
    printKmer(_printer->file(), _kmer, _ot->_printACGTorder);

  if (_stats)
    _stats->addValue(_kmer._val);

  //  We have now loaded the nextMer; return and let our clients query us.

  return(true);
}
