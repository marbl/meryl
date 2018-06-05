
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#include "meryl.H"



void
merylOperation::findMinCount(void) {
  _count = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_actCount[ii] < _count)
      _count = _actCount[ii];
}



void
merylOperation::findMaxCount(void) {
  _count = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_count < _actCount[ii])
      _count = _actCount[ii];
}



void
merylOperation::findSumCount(void) {
  _count = 0;
  for (uint32 ii=0; ii<_actLen; ii++)
    _count += _actCount[ii];
}



bool
merylOperation::nextMer(bool isRoot) {

  char  kmerString[256];

  //  Find the smallest kmer in the _inputs, and save their counts in _actCount.
  //  Mark which input was used in _actIndex.

  if (_verbosity >= sayDetails) {
    fprintf(stderr, "\n");
    fprintf(stderr, "merylOp::nextMer()-- STARTING for operation %s\n",
            toString(_operation));
  }

  //  If we're counting, do that entirely right now, then reset so we
  //  can iterate through the kmers we just counted.
  //
  //  It's a bit more complicated than I like.  We need to close the output
  //  (so all the data gets written and an index created) before opening the input,
  //  but to close t

  if ((_operation == opCount) ||
      (_operation == opCountForward) ||
      (_operation == opCountReverse)) {
    if (_kmer.merSize() <= 16)
      countSimple();
    else
      count();

    //  Done with the inputs, so forget about them.

    clearInputs();

    //  Remember the name of the data we just created.

    char  dataName[FILENAME_MAX+1];

    strncpy(dataName, _output->filename(), FILENAME_MAX);

    //  Close the output and forget about it.

    delete _output;
    _output = NULL;

    //  If we're the root node, nobody is going to read our kmers,
    //  and we can just return that there are no kmers.

    if (isRoot)
      return(false);

    //  Otherwise, add the output we just made as an input.

    if (_verbosity >= sayConstruction)
      fprintf(stderr, "merylOp::nextMer()-- CONVERTING '%s' to '%s'.\n",
              toString(opCount), toString(opPassThrough));

    _operation = opPassThrough;

    addInput(new kmerCountFileReader(dataName));
  }

  if (_verbosity >= sayEverything)
    for (uint32 ii=0; ii<_inputs.size(); ii++)
      fprintf(stderr, "merylOp::nextMer()--   CURRENT STATE: input %s kmer %s count " F_U64 " %s\n",
              _inputs[ii]->_name,
              _inputs[ii]->_kmer.toString(kmerString),
              _inputs[ii]->_count,
              _inputs[ii]->_valid ? "valid" : "INVALID");

  //  Grab the next mer for every input that was active in the last iteration.
  //  (on the first call, all inputs were 'active' last time)
  //
  for (uint32 ii=0; ii<_actLen; ii++) {
    if (_verbosity >= sayDetails)
      fprintf(stderr, "merylOp::nextMer()-- CALL NEXTMER on input actIndex " F_U32 "\n", _actIndex[ii]);
    _inputs[_actIndex[ii]]->nextMer();
  }

  _actLen = 0;

  //  Log.

  if (_verbosity >= sayEverything)
    for (uint32 ii=0; ii<_inputs.size(); ii++)
      fprintf(stderr, "merylOp::nextMer()--   BEFORE OPERATION: input %s kmer %s count " F_U64 " %s\n",
              _inputs[ii]->_name,
              _inputs[ii]->_kmer.toString(kmerString),
              _inputs[ii]->_count,
              _inputs[ii]->_valid ? "valid" : "INVALID");

  //  Build a list of the inputs that have the smallest kmer.

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->_valid == false)
      continue;

    //  If we have no active kmer, or the input kmer is smaller than the one we
    //  have, reset the list.

    if ((_actLen == 0) ||
        (_inputs[ii]->_kmer < _kmer)) {
      _actLen = 0;
      _kmer              = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      if (_verbosity >= sayDetails)
        fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s. reset\n", _kmer.toString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, if the input kmer is the one we have, save the count to the list.

    else if (_inputs[ii]->_kmer == _kmer) {
      //_kmer             = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      if (_verbosity >= sayDetails)
        fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s\n", _kmer.toString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, the input kmer comes after the one we're examining, ignore it.

    else {
    }
  }

  //  If no active kmers, we're done.

  if (_actLen == 0) {
    if (_verbosity >= sayDetails) {
      fprintf(stderr, "merylOp::nextMer()-- No inputs found, all done here.\n");
      fprintf(stderr, "\n");
    }
    _valid = false;
    return(false);
  }

  //  Otherwise, active kmers!  Figure out what the count should be.

  if (_verbosity >= sayDetails)
    fprintf(stderr, "merylOp::nextMer()-- op %s activeLen " F_U32 " kmer %s\n", toString(_operation), _actLen, _kmer.toString(kmerString));

  //  If math-subtract gets implemented, use negative-zero to mean "don't output" and positive-zero
  //  to mean zero.  For now, count=0 means don't output.

  //  Set the count to zero, meaning "don't output the kmer".  Intersect depends on this,
  //  skipping most of it's work if all files don't have the kmer.
  _count = 0;

  switch (_operation) {
    case opCount:
    case opCountForward:
    case opCountReverse:
      fprintf(stderr, "ERROR: got %s, but shouldn't have.\n", toString(_operation));
      assert(0);
      break;

    case opPassThrough:
    case opUnion:                           //  Union
      _count = 1;
      break;

    case opUnionMin:                        //  Union, retain smallest count
      findMinCount();
      break;

    case opUnionMax:                        //  Union, retain largest count
      findMaxCount();
      break;

    case opUnionSum:                        //  Union, sum all counts
      findSumCount();
      break;

    case opIntersect:                       //  Intersect
      if (_actLen == _inputs.size())
        _count = 1;
      break;

    case opIntersectMin:                    //  Intersect, retain smallest count
      if (_actLen == _inputs.size())
        findMinCount();
      break;

    case opIntersectMax:                    //  Intersect, retain largest count
      if (_actLen == _inputs.size())
        findMaxCount();
      break;

    case opIntersectSum:                    //  Intersect, sum all counts
      if (_actLen == _inputs.size())
        findSumCount();
      break;

    case opDifference:
      break;

    case opSymmetricDifference:
      break;

    case opPrint:
      if (_inputs.size() != 1)
        fprintf(stderr, "merylOp::nextMer()-- ERROR: 'print' can operate on one input only; this has " F_SIZE_T " inputs.\n",
                _inputs.size()), exit(1);

      fprintf(stdout, "%s\t" F_U64 "\n", _kmer.toString(kmerString), _actCount[0]);
      break;

    case opNothing:
      break;
  }  

  //  If flagged for output, output!

  if (_verbosity >= sayDetails) {
    fprintf(stderr, "merylOp::nextMer()-- FINISHED for operation %s with kmer %s count " F_U64 "%s\n",
            toString(_operation), _kmer.toString(kmerString), _count, ((_output != NULL) && (_count != 0)) ? " OUTPUT" : "");
    fprintf(stderr, "\n");
  }

  if ((_output != NULL) &&
      (_count  != 0))
    _output->addMer(_kmer, _count);

  return(true);
}
