
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


merylOperation::merylOperation(merylOp op, uint32 k, uint64 numMers, uint32 threads, uint64 memory) {
  fprintf(stderr, "Create operation '%s'\n", toString(op));

  _operation     = op;
  _k             = k;
  _numMers       = numMers;
  _threads       = threads;
  _memory        = memory;

  _outputName[0] = 0;
  _output        = NULL;

  _actLen        = 0;
  _actCount      = new uint64 [1024];
  _actIndex      = new uint32 [1024];

  _count         = 0;
  _valid         = true;
}


merylOperation::~merylOperation() {
  fprintf(stderr, "Destroy op %s\n", toString(_operation));

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  delete    _output;
  delete [] _actCount;
  delete [] _actIndex;
}




void
merylOperation::addInput(merylOperation *operation) {
  fprintf(stderr, "Adding input from operation '%s' to operation '%s'\n",
          toString(operation->_operation), toString(_operation));

  _inputs.push_back(new merylInput(operation));
  _actIndex[_actLen++] = _inputs.size() - 1;
}


void
merylOperation::addInput(char *name, kmerCountFileReader *reader) {
  fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
          name, toString(_operation));

  _inputs.push_back(new merylInput(name, reader));
  _actIndex[_actLen++] = _inputs.size() - 1;
}


void
merylOperation::addInput(char *name, dnaSeqFile *sequence) {
  fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
          name, toString(_operation));

  _inputs.push_back(new merylInput(name, sequence));
  _actIndex[_actLen++] = _inputs.size() - 1;
}



void
merylOperation::addOutput(char *name, kmerCountFileWriter *writer) {
  fprintf(stderr, "Adding output to file '%s' from operation '%s'\n",
          name, toString(_operation));

  strncpy(_outputName, name, FILENAME_MAX);
  _output = writer;
}




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



char const *
toString(merylOp op) {
  switch (op) {
    case opCount:                return("opCount");                break;
    case opCountForward:         return("opCountForward");         break;
    case opCountReverse:         return("opCountReverse");         break;
    case opUnion:                return("opUnion");                break;
    case opUnionMin:             return("opUnionMin");             break;
    case opUnionMax:             return("opUnionMax");             break;
    case opUnionSum:             return("opUnionSum");             break;
    case opIntersect:            return("opIntersect");            break;
    case opIntersectMin:         return("opIntersectMin");         break;
    case opIntersectMax:         return("opIntersectMax");         break;
    case opIntersectSum:         return("opIntersectSum");         break;
    case opDifference:           return("opDifference");           break;
    case opSymmetricDifference:  return("opSymmetricDifference");  break;
    case opComplement:           return("opComplement");           break;
    case opNothing:              return("opNothing");              break; 
  }

  assert(0);
  return(NULL);
}

