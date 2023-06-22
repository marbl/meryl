
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


merylOpTemplate::merylOpTemplate(uint32 ident) {
  _ident = ident;
}

merylOpTemplate::~merylOpTemplate() {

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  _inputs.clear();

  for (uint32 ii=0; ii<64; ii++)
    delete _computes[ii];

  delete    _writer;
  delete    _printer;
  delete [] _printerName;

  delete    _statsFile;
  delete    _histoFile;
}



void
merylOpTemplate::addInputFromOp(merylOpTemplate *otin, std::vector<char const *> &err) {

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addInputFromOp()-- action #%u <- input from action #%u\n",
            _ident, otin->_ident);

  _inputs.push_back(new merylInput(otin));
}


void
merylOpTemplate::addInputFromDB(char const *dbName, std::vector<char const *> &err) {

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addInputFromDB()-- action #%u <- input from file '%s'\n",
            _ident, dbName);

  _inputs.push_back(new merylInput(new merylFileReader(dbName)));
}



//  Probably should go away if we can get rid of the err's on the other addInput methods.
void
merylOpTemplate::addInputFromDB(char const *dbName) {
  _inputs.push_back(new merylInput(new merylFileReader(dbName)));
}


void
merylOpTemplate::addInputFromSeq(char const *sqName, bool doCompression, std::vector<char const *> &err) {

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addInputFromSeq()-- action #%u <- input from file '%s'\n",
            _ident, sqName);

  _inputs.push_back(new merylInput(new dnaSeqFile(sqName), doCompression));
}



void
merylOpTemplate::addInputFromCanu(char const *sqName, uint32 segment, uint32 segmentMax, std::vector<char const *> &err) {

#ifdef CANU

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addInputFromCanu()-- action #%u <- input from sqStore '%s'\n",
            _ident, sqName);

  _inputs.push_back(new merylInput(new sqStore(sqName), segment, segmentMax));

#endif
}



void
merylOpTemplate::addOutput(char const *wrName, std::vector<char const *> &err) {

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addOutput()-- action #%u -> output to database '%s'\n",
            _ident, wrName);

  if (_writer) {
    sprintf(err, "Operation #%u already writes output to '%s',", _ident, _writer->filename());
    sprintf(err, "  can't add another output to '%s'!", wrName);
    sprintf(err, "");
    return;
  }

  _writer = new merylFileWriter(wrName);
}



//  Add a printer to the template command.
//
//  Three outcomes:
//
//   - The name is '-' and output will go to stdout.
//
//   - The name does not contain two or more '#' symbols, and output will
//     go to a single file with that name.
//
//   - The name does contain two or more '#' symbols, and output will
//     go to 64 files, one per slice, replacing ##'s with digits.
//
//  The first two cases both define _printer, while the third does not.
//  merylOpCompute::addPrinter() uses this to decide if it should open
//  per-slice output or use the global output.
void
merylOpTemplate::addPrinter(char const *prName, bool ACGTorder, std::vector<char const *> &err) {

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addOutput()-- action #%u -> print to '%s'\n",
            _ident, (prName == nullptr) ? "(stdout)" : prName);

  //  Fail if we've already got a printer assigned.

  if (_printerName) {
    sprintf(err, "Operation #%u is already printing to '%s',", _ident, _printerName);
    sprintf(err, "  can't add another output to '%s'.", prName);
    sprintf(err, "");
    return;
  }

  //  Decide if this is to stdout.

  if ((prName == nullptr) || (strcmp("-", prName) == 0)) {
    _printer        = new compressedFileWriter(prName);
    _printerName    = duplicateString("(stdout)");
    _printACGTorder = ACGTorder;
    return;
  }

  //  Decide if this is a normal file.

  uint32  len = 0;

  for (char const *suf = strchr(prName, '#'); ((suf) && (*suf == '#')); suf++)
    len++;

  if (len < 2) {
    _printer        = new compressedFileWriter(prName);
    _printerName    = duplicateString(prName);
    _printACGTorder = ACGTorder;
  }
  else {
    _printer        = nullptr;
    _printerName    = duplicateString(prName);
    _printACGTorder = ACGTorder;
  }
}



//
//  Add 'histogram' or 'statistics' output to this operation.
//  Fails if 

void
merylOpTemplate::addHistogram(char const *hiName, bool asStats, std::vector<char const *> &err) {

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "addOutput()-- action #%u -> %s to '%s'\n",
            _ident, (asStats == true) ? "statistics" : "histogram", hiName);

  if (asStats == true) {
    if (_statsFile != nullptr) {
      sprintf(err, "Operation #%u already has 'statistics' output to file '%s',", _ident, _statsFile->filename());
      sprintf(err, "  can't add another output to file '%s'.", hiName);
      sprintf(err, "");
      return;
    }
    _statsFile = new compressedFileWriter(hiName);
  }

  if (asStats == false) {
    if (_histoFile != nullptr) {
      sprintf(err, "Operation #%u already has 'histogram' output to file '%s',", _ident, _histoFile->filename());
      sprintf(err, "  can't add another output to file '%s'.", hiName);
      sprintf(err, "");
      return;
    }
    _histoFile = new compressedFileWriter(hiName);
  }
}



//  This is called by merylCommandBuilder after the entire tree has been
//  built.
//
//  The purpose is to allow the filters a chance to figure out what
//  'all' means and to fail if there are any errors (like asking for input
//  4 when there are only 3 available).
//
//  It is NOT intended to do any processing, like examining histograms to
//  choose thresholds.  That is done in initializeTemplate().
//
void
merylOpTemplate::finalizeTemplateInputs(std::vector<char const *> &err) {

  for (uint32 f1=0; f1<_filter    .size(); f1++)
    for (uint32 f2=0; f2<_filter[f1].size(); f2++)
      _filter[f1][f2].finalizeFilterInputs(this, err);
}



//  This is called by spawnThreads() just before the template is copied
//  into (64) merylOpCompute objects.
//
//  The purpose is to allow inputs and outputs to be opened before threads
//  are spawned, and to allow the flters a chance to examine histograms or
//  anything else on disk.
//
void
merylOpTemplate::finalizeTemplateParameters(void) {

  for (uint32 ii=0; ii<_inputs.size(); ii++)       //  Forward the request to any inputs
    if (_inputs[ii]->_template != nullptr)         //  that are actions.
      _inputs[ii]->_template->finalizeTemplateParameters();

  if (_writer)                                   //  Create the master output object.  We'll later
    _writer->initialize(0, false);               //  request per-thread writer objects from this.

  //  Any filters that need to query meryl databases for parameters
  //  should do so now.

  for (uint32 f1=0; f1<_filter    .size(); f1++)
  for (uint32 f2=0; f2<_filter[f1].size(); f2++)
    _filter[f1][f2].finalizeFilterParameters(this);
}



void
merylOpTemplate::finishAction(void) {

  //  Forward the request to any inputs that are actions.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    if (_inputs[ii]->_template != nullptr)
      _inputs[ii]->_template->finishAction();

  //  Gather stats from the compute threads then output.

  if ((_statsFile != nullptr) ||
      (_histoFile != nullptr)) {
    merylHistogram  stats;

    for (uint32 ss=0; ss<64; ss++)
      stats.insert( _computes[ss]->_stats );

    if (_statsFile != nullptr)
      stats.reportStatistics(_statsFile->file());

    if (_histoFile != nullptr)
      stats.reportHistogram(_histoFile->file());
  }

  //  Delete the compute objets.

  //  Anything else should be done in the various destructors.
}





void
merylOpTemplate::doCounting(uint64 allowedMemory,
                            uint32 allowedThreads) {
  if (_counting == nullptr)
    return;

  if (_writer == nullptr)
    return;

  _counting->_lConstant = _labelConstant;

  //  Call the counting method.
  _counting->doCounting(_inputs, allowedMemory, allowedThreads, _writer);

  if (_onlyConfig == true)   //  If only configuring, stop now.
    return;

  //  Convert this op into a pass through
#warning NEED TO RESET COUNTING OPERATION TO PASS THROUGH
  //  Fiddle with the operation.
  //   - remove the output; it's already been written.
  //   - remove all the inputs
  //   - convert the operation to a simple 'pass through'
  //   - add the counted output as an input

  char name[FILENAME_MAX + 1];

  strncpy(name, _writer->filename(), FILENAME_MAX + 1);   //  know which input to open later.

  //  Close the output and forget about it.
  delete _writer;
  _writer = nullptr;

  //  Close the inputs and forget about them too.
  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];
  _inputs.clear();

  //  But now remember what that output was and make it an input.
  addInputFromDB(name);
};
