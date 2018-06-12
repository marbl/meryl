
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


bool            merylOperation::_showProgress = false;
merylVerbosity  merylOperation::_verbosity    = sayStandard;


merylOperation::merylOperation(merylOp op, uint32 threads, uint64 memory) {
  _operation     = op;

  _parameter     = UINT64_MAX;
  _expNumKmers   = 0;

  _maxThreads    = threads;
  _maxMemory     = memory;

  _output        = NULL;

  _actLen        = 0;
  _actCount      = new uint64 [1024];
  _actIndex      = new uint32 [1024];

  _count         = 0;
  _valid         = true;
}


merylOperation::~merylOperation() {

  clearInputs();

  delete    _output;
  delete [] _actCount;
  delete [] _actIndex;
}




void
merylOperation::clearInputs(void) {

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  _inputs.clear();

  _actLen = 0;
}



void
merylOperation::checkInputs(const char *name) {

  if ((_actLen > 1) && ((_operation == opLessThan)    ||
                        (_operation == opGreaterThan) ||
                        (_operation == opEqualTo)     ||
                        (_operation == opPrint))) {
    fprintf(stderr, "merylOp::addInput()-- ERROR: Can't add input '%s' to operation '%s': only one input supported.\n",
            name, toString(_operation));
    exit(1);
  }
}



void
merylOperation::addInput(merylOperation *operation) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from operation '%s' to operation '%s'\n",
            toString(operation->_operation), toString(_operation));

  if ((_inputs.size() > 0) && (_operation == opPrint))
    fprintf(stderr, "ERROR: 'print' operation can have exactly one input.\n"), exit(1);

  _inputs.push_back(new merylInput(operation));
  _actIndex[_actLen++] = _inputs.size() - 1;

  checkInputs(toString(operation->getOperation()));
}


void
merylOperation::addInput(kmerCountFileReader *reader) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
            reader->filename(), toString(_operation));

  if ((_inputs.size() > 0) && (_operation == opPrint))
    fprintf(stderr, "ERROR: 'print' operation can have exactly one input.\n"), exit(1);

  _inputs.push_back(new merylInput(reader->filename(), reader));
  _actIndex[_actLen++] = _inputs.size() - 1;

  checkInputs(reader->filename());
}


void
merylOperation::addInput(dnaSeqFile *sequence) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
            sequence->filename(), toString(_operation));

  if ((_inputs.size() > 0) && (_operation == opPrint))
    fprintf(stderr, "ERROR: 'print' operation can have exactly one input.\n"), exit(1);

  _inputs.push_back(new merylInput(sequence->filename(), sequence));
  _actIndex[_actLen++] = _inputs.size() - 1;

  checkInputs(sequence->filename());
}



void
merylOperation::addOutput(kmerCountFileWriter *writer) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding output to file '%s' from operation '%s'\n",
            writer->filename(), toString(_operation));

  if (_output)
    fprintf(stderr, "ERROR: already have an output set!\n"), exit(1);

  _output = writer;
}




char const *
toString(merylOp op) {
  switch (op) {
    case opCount:                return("opCount");                break;
    case opCountForward:         return("opCountForward");         break;
    case opCountReverse:         return("opCountReverse");         break;
    case opPassThrough:          return("opPassThrough");          break;
    case opLessThan:             return("opLessThan");             break;
    case opGreaterThan:          return("opGreaterThan");          break;
    case opEqualTo:              return("opEqualTo");              break;
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
    case opPrint:                return("opPrint");                break;
    case opNothing:              return("opNothing");              break; 
  }

  assert(0);
  return(NULL);
}

