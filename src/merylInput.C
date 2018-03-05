
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



merylInput::merylInput(merylOperation *o) {
  _operation   = o;
  _stream      = NULL;
  _sequence    = NULL;
  _count       = 0;
  _valid       = false;

  strncpy(_name, toString(_operation->getOperation()), FILENAME_MAX);
}



merylInput::merylInput(char *n, kmerStreamReader *s) {
  _operation   = NULL;
  _stream      = s;
  _sequence    = NULL;
  _count       = 0;
  _valid       = false;

  strncpy(_name, n, FILENAME_MAX);
}



merylInput::merylInput(char *n, dnaSeqFile *f) {
  _operation   = NULL;
  _stream      = NULL;
  _sequence    = f;
  _count       = 0;
  _valid       = true;    //  Trick nextMer into doing something without a valid mer.

  strncpy(_name, n, FILENAME_MAX);
}



merylInput::~merylInput() {
  fprintf(stderr, "Destroy input %s\n", _name);
  delete _stream;
  delete _operation;
}



void
merylInput::nextMer(void) {
  char kmerString[256];

  if (_stream) {
    fprintf(stderr, "merylIn::nextMer('%s')--\n", _name);

    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
    _count = _stream->theCount();

    fprintf(stderr, "merylIn::nextMer('%s')-- now have valid=" F_U32 " kmer %s count " F_U64 "\n",
            _name, _valid, _kmer.toString(kmerString), _count);
    fprintf(stderr, "\n");
  }

  if (_operation) {
    fprintf(stderr, "merylIn::nextMer(%s)--\n", _name);

    _valid = _operation->nextMer();
    _kmer  = _operation->theFMer();
    _count = _operation->theCount();

    fprintf(stderr, "merylIn::nextMer(%s)-- now have valid=" F_U32 " kmer %s count " F_U64 "\n",
            _name, _valid, _kmer.toString(kmerString), _count);
    fprintf(stderr, "\n");
  }

  if (_sequence) {
    fprintf(stderr, "merylIn::nextMer(%s)--\n", _name);
  }
}

