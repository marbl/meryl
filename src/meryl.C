
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



int
main(int argc, char **argv) {
  stack<merylOperation *>   opStack;

  uint32                    merSize    = 0;
  uint64                    numMers    = 0;
  uint32                    maxThreads = (uint32)1;
  uint64                    maxMem     = (uint64)2 * 1024 * 1024 * 1024;

  uint32                    outputArg  = UINT32_MAX;

  vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    char    mcidx[FILENAME_MAX];
    char    mcdat[FILENAME_MAX];
    char    opt[FILENAME_MAX];
    uint32  optLen = strlen(argv[arg]);

    //uint32   creating    = 0;
    uint32   terminating = 0;

    //  Scan for options.

    if      (strcmp(argv[arg], "-h") == 0) {
      err.push_back(NULL);
      continue;
    }

    else if (strcmp(argv[arg], "-k") == 0) {
      merSize = strtouint32(argv[++arg]);
      kmerTiny::setSize(merSize);
      continue;
    }

    else if (strcmp(argv[arg], "-n") == 0) {
      numMers = strtouint64(argv[++arg]);
      continue;
    }

    else if (strcmp(argv[arg], "-m") == 0) {
      maxMem = strtouint64(argv[++arg]) * 1024 * 1024 * 1024;
      continue;
    }

    else if (strcmp(argv[arg], "-t") == 0) {
      maxThreads = strtouint32(argv[++arg]);
      continue;
    }


    //  Save a copy of the argument, just in case it's a filename, before
    //  we munge out any brackets.

    strncpy(opt, argv[arg], FILENAME_MAX);

    //  If we have [ as the first character, make a new operation (with no
    //  operation set yet) and add it to the inputs, then push on the stack.
    //
    //  We can get 0 or 1 open bracket at a time.  Seeing two in a row is an error,
    //  but it isn't caught.

    if (opStack.empty() == true) {
      //creating = true;
    }

    if (opt[0] == '[') {
      strncpy(opt, argv[arg]+1, FILENAME_MAX);
      optLen--;

      fprintf(stderr, "CREATENEW operation\n");

      //creating = true;
    }

    //  If we have a ] as the last character, strip it off and remember.
    //
    //  We can get any number of closing brackets.

    while (opt[optLen-1] == ']') {
      opt[optLen-1] = 0;
      optLen--;

      fprintf(stderr, "TERMINATE operation\n");

      terminating++;
    }

    //  Now that brackets are stripped, make meryl database names for the arg.

    snprintf(mcidx, FILENAME_MAX, "%s.mcidx", opt);
    snprintf(mcdat, FILENAME_MAX, "%s.mcdat", opt);

    merylOp             op       = opNothing;
    kmerStreamWriter  *writer   = NULL;
    kmerStreamReader  *reader   = NULL;
    dnaSeqFile         *sequence = NULL;


    if      (opt[0] == 0)
      ;  //  Got a single bracket, nothing to do here except make it not be an error.

    else if (strcmp(opt, "union") == 0)                  op = opUnion;
    else if (strcmp(opt, "union-min") == 0)              op = opUnionMin;
    else if (strcmp(opt, "union-max") == 0)              op = opUnionMax;
    else if (strcmp(opt, "union-sum") == 0)              op = opUnionSum;
    else if (strcmp(opt, "intersect") == 0)              op = opIntersect;
    else if (strcmp(opt, "intersect-min") == 0)          op = opIntersectMin;
    else if (strcmp(opt, "intersect-max") == 0)          op = opIntersectMax;
    else if (strcmp(opt, "intersect-sum") == 0)          op = opIntersectSum;
    else if (strcmp(opt, "difference") == 0)             op = opDifference;
    else if (strcmp(opt, "symmetric-difference") == 0)   op = opSymmetricDifference;
    else if (strcmp(opt, "complement") == 0)             op = opComplement;
    else if (strcmp(opt, "count") == 0)                  op = opCount;
    else if (strcmp(opt, "count-forward") == 0)          op = opCountForward;
    else if (strcmp(opt, "count-reverse") == 0)          op = opCountReverse;

    //  If we see 'output', flag the next arg as being the output name.
    //  If this arg is flagged as output, add an output using the bracket-stripped name.
    //  If this arg is a valid meryl file, make it an input.

    else if (strcmp(opt, "output") == 0)
      outputArg = arg+1;

    else if (arg == outputArg)
      writer = new kmerStreamWriter(opt, merSize, 0, merSize/2, false);

    else if ((opStack.size() > 0) &&
             (opStack.top()->isCounting() == false) &&
             (AS_UTL_fileExists(mcidx)    == true) &&
             (AS_UTL_fileExists(mcdat)    == true))
      reader = new kmerStreamReader(opt);

    else if ((opStack.size() > 0) &&
             (opStack.top()->isCounting() == true) &&
             (AS_UTL_fileExists(opt)      == true))
      sequence = new dnaSeqFile(opt);

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Don't know what to do with '%s'.\n", opt);
      err.push_back(s);
    }


    //  Create a new operation or set inputs or output.  Or do nothing.

    if (op != opNothing) {
      merylOperation *newOp = new merylOperation(op, merSize, numMers, maxThreads, maxMem);

      if (opStack.empty() == false)
        opStack.top()->addInput(newOp);

      opStack.push(newOp);
      op = opNothing;
    }

    if ((writer != NULL) &&                    //  If nothing on the stack, wait for an
        (opStack.size() > 0)) {                //  operation to show up.
      opStack.top()->addOutput(opt, writer);
      writer = NULL;
    }

    if (reader != NULL) {                      //  A reader exists only if an
      opStack.top()->addInput(opt, reader);    //  operation is on the stack.
      reader = NULL;
    }

    if (sequence != NULL) {                    //  A sequence input exists only if
      opStack.top()->addInput(opt, sequence);  //  an operation is on the stack.
      sequence = NULL;
    }

    //  Now that we're done with the operation, if it was terminating, pop it off the stack.

    for (; terminating > 0; terminating--) {
      fprintf(stderr, "TERMINATE op '%s'\n", toString(opStack.top()->getOperation()));

      opStack.pop();
    }
  }

  //  If any errors, fail.

  if ((argc == 1) ||        //  No commands
      (err.size() > 0)) {   //  Errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -k <K>       work with mers of size K bases\n");
    fprintf(stderr, "  -n <N>       expect N mers in the input\n");
    fprintf(stderr, "  -m <M>       use no more than M gigabytes of memory\n");
    fprintf(stderr, "  -t <T>       use no more than T compute threads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  //  Now just walk through the kmers until nothing is left.

  merylOperation *op = opStack.top();

  fprintf(stderr, "START\n");
  fprintf(stderr, "START operation %s\n", toString(op->getOperation()));
  fprintf(stderr, "START\n");

  //  The counting operations need to be special cased.

  //if ((op->getOperation() == opCount) ||
  //    (op->getOperation() == opCountForward) ||
  //    (op->getOperation() == opCountReverse))
  //  op->count();

  //  Now just walk through the kmers until nothing is left.

  while (op->nextMer() == true)
    ;

  //  Done!

  fprintf(stderr, "DONE\n");
  fprintf(stderr, "DONE operation %s\n", toString(op->getOperation()));
  fprintf(stderr, "DONE\n");

  delete op;

  return(0);
}
