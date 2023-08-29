
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



//  Process options for any current command.
//
//  Some of these require a command to exist, some can
//  be deferred until the command is encountered.
//
bool
merylCommandBuilder::isOptionWord(void) {
  merylOpTemplate  *op = getCurrent();

#if 0
  bool  isNum = isDecNumber(_optString, 0);

  //  Compatibility mode might need this.
  //     greater-than ### ----> add a new selector action 'select value:>#'
  if ((op->needsThreshold() == true) && (isNum == true)) {
    op->setThreshold(strtouint64(_optString));
    return(true);
  }

  if ((op->needsConstant() == true) && (isNum == true)) {
    op->setConstant(strtouint64(_optString));
    return(true);
  }
#endif

  //  Some options have no value, unfortunately.

  if (strcmp(_optString, "compress") == 0) {
    _doCompression = true;
    return(true);
  }

  //  The rest should be key=value options.

  KeyAndValue   kv(_optString);
  char const   *key = kv.key();
  char const   *val = kv.value();

  if ((key == nullptr) ||
      (val == nullptr))
    return(false);

  uint32 valLen = strlen(val);
  uint32 val32  = strtouint32(val);
  uint64 val64  = strtouint64(val);
  double valDB  = strtodouble(val);

  //  Number of kmers expected for counting.
  if (strcmp(key, "n") == 0) {
    if (op->_isCounting == true)
      op->_counting->setExpectedNumberOfKmers(val64);
    else
      sprintf(_errors, "option '%s' encountered for non-counting operation.", _optString);
    return(true);
  }

  //  A suffix to select kmers by when counting.
  if (strcmp(key, "count-suffix") == 0) {
    if (op->_isCounting == true)
      op->_counting->setCountSuffix(val);
    else
      sprintf(_errors, "option '%s' encountered for non-counting operation.", _optString);
    return(true);
  }

  //  Threshold values for less-than, etc, specifed as a fraction of the
  //  total distinct kmers, or as a word-frequency, or as an absolute count.
#if 0
  if ((strcmp(key, "d") == 0) ||
      (strcmp(key, "distinct") == 0)) {
    op->setFractionDistinct(valDB);
    return(true);
  }

  if ((strcmp(key, "f") == 0) ||
      (strcmp(key, "word-frequency") == 0)) {
    op->setWordFrequency(valDB);
    return(true);
  }

  if ((strcmp(key, "t") == 0) ||            //  See above for special case of this!
      (strcmp(key, "threshold") == 0)) {
    op->setThreshold(val64);
    return(true);
  }
#endif

  //  Segment of input, for counting from seqStore.  Useless otherwise.
  if ((strcmp(key, "segment") == 0) &&
      (isDecNumber(val, '/'))) {
    decodeRange(val, _segment, _segmentMax);
#ifndef CANU
    sprintf(_errors, "option '%s' available only with Canu support.", _optString);
#endif
    return(true);
  }

  //  If nothing triggered, we don't recognize it as an option.  Maybe the
  //  filename had an '=' in it?

  return(false);
}
