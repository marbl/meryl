
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


bool
isDigit(char c) {
  return(('0' <= c) && (c <= '9'));
}

bool
isNumber(char *s) {

  if (s == NULL)
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if ((isDigit(s[ii]) == false) &&
        (s[ii] != '.'))
      return(false);

  return(true);
}




int
main(int argc, char **argv) {
  uint32                    optStringLen = 0;
  char                      optString[FILENAME_MAX+1];
  char                      inoutName[FILENAME_MAX+1];
  char                      indexName[FILENAME_MAX+1];

  stack<merylOperation *>   opStack;
  merylOp                   opName         = opNothing;
  uint32                    outputArg      = UINT32_MAX;
  kmerCountFileWriter      *writer         = NULL;
  kmerCountFileReader      *reader         = NULL;
  dnaSeqFile               *sequence       = NULL;

  uint32                    terminating    = 0;

  uint32                    physThreads    = omp_get_max_threads();     //  Absolute maximum limits on
  uint64                    physMemory     = getPhysicalMemorySize();   //  memory= and threads= values.

  uint32                    allowedThreads = physThreads;               //  Global limits, if memory= or
  uint64                    allowedMemory  = physMemory;                //  threads= is set before any operation.


  vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {

    //  Save a few copies of the command line word.

    optStringLen = strlen(argv[arg]);

    strncpy(optString, argv[arg], FILENAME_MAX);
    strncpy(inoutName, argv[arg], FILENAME_MAX);
    strncpy(indexName, argv[arg], FILENAME_MAX);
    strncat(indexName, "/merylIndex", FILENAME_MAX - optStringLen - 1);

    //  Scan for options.  If any trigger, set the option and move on to the next word on the command line.

    if     ((strcmp(optString, "-h")   == 0) ||
            (strcmp(optString, "help") == 0)) {
      err.push_back(NULL);
      continue;
    }

    //  Kmer size.
    else if ((optStringLen > 2) &&
             (strncmp(optString, "k=", 2) == 0) &&
             (isNumber(optString + 2) == true)) {
      kmerTiny::setSize(strtouint32(optString + 2));
      continue;
    }

    //  Number of kmers expected for counting.
    else if ((opStack.size() > 0) &&
             (optStringLen > 2) &&
             (strncmp(optString, "n=", 2) == 0) &&
             (isNumber(optString + 2) == true)) {
      opStack.top()->setExpectedNumberOfKmers(strtouint64(optString + 2));
      continue;
    }

    //  Threshold values for less-than, greater-than and equal-to are just a number.
    else if ((opStack.size() > 0) &&
             (isNumber(optString))) {
      opStack.top()->setParameter(strtouint64(optString));
      continue;
    }

    //  Memory limit, either global or per-task.
    else if ((optStringLen > 7) &&
             (strncmp(optString, "memory=", 7) == 0) &&
             (isNumber(optString + 7) == true)) {
      uint64 memory = (uint64)(strtodouble(optString + 7) * 1024 * 1024 * 1024);

      if (memory > physMemory) {
        char *s = new char [1024];
        snprintf(s, 1024, "Requested memory '%s' (GB) is more than physical memory %.2f GB.",
                 optString, physMemory / 1024.0 / 1024.0 / 1024.0);
        err.push_back(s);
      }

      if (opStack.size() == 0)
        allowedMemory = memory;
      else
        opStack.top()->setMemoryLimit(memory);
             
      continue;
    }

    //  Thread limit, either global or per-task.
    else if ((optStringLen > 8) &&
             (strncmp(optString, "threads=", 8) == 0) &&
             (isNumber(optString + 8) == true)) {
      uint32 threads = strtouint32(optString + 8);

      if (opStack.size() == 0)
        allowedThreads = threads;
      else
        opStack.top()->setThreadLimit(threads);

      //omp_set_num_threads(allowedThreads);

      continue;
    }

    //  Global thread limit
    else if ((opStack.size() == 0) &&
             (optStringLen > 8) &&
             (strncmp(optString, "threads=", 8) == 0) &&
             (isNumber(optString + 8) == true)) {
      continue;
    }



    else if (strncmp(optString, "-V", 2) == 0) {      //  Anything that starts with -V
      for (uint32 vv=1; vv<strlen(optString); vv++)   //  increases verbosity by the
        merylOperation::increaseVerbosity();          //  number of letters.
      continue;
    }

    else if (strcmp(optString, "-Q") == 0) {
      merylOperation::beQuiet();
      continue;
    }

    else if (strcmp(optString, "-P") == 0) {
      merylOperation::showProgress();
      continue;
    }



    //  Ignore '[' at the start of the string.  Their purpose is in the matching ']' which
    //  tells us to stop adding inputs to the current command.
    //
    //  There should only be one opening bracket.

    if (optString[0] == '[') {
      for (uint32 ii=0; ii<optStringLen; ii++)
        optString[ii] = optString[ii+1];

      optStringLen--;
    }

    //  If we have a ] as the last character, strip it off and remember that we need to
    //  close the command on the stack after we process this arg.
    //
    //  We can get any number of closing brackets.

    while ((optStringLen > 0) &&
           (optString[optStringLen-1] == ']')) {
      optString[optStringLen-1] = 0;
      optStringLen--;

      terminating++;
    }

    //  Now, parse this word.  Decide if it's a new operation, or an output name, or an input file.

    if      (0 == optStringLen)
      ;  //  Got a single bracket, nothing to do here except make it not be an error.

    else if (0 == strcmp(optString, "count"))                  opName = opCount;
    else if (0 == strcmp(optString, "count-forward"))          opName = opCountForward;
    else if (0 == strcmp(optString, "count-reverse"))          opName = opCountReverse;
    else if (0 == strcmp(optString, "less-than"))              opName = opLessThan;
    else if (0 == strcmp(optString, "greater-than"))           opName = opGreaterThan;
    else if (0 == strcmp(optString, "equal-to"))               opName = opEqualTo;
    else if (0 == strcmp(optString, "union"))                  opName = opUnion;
    else if (0 == strcmp(optString, "union-min"))              opName = opUnionMin;
    else if (0 == strcmp(optString, "union-max"))              opName = opUnionMax;
    else if (0 == strcmp(optString, "union-sum"))              opName = opUnionSum;
    else if (0 == strcmp(optString, "intersect"))              opName = opIntersect;
    else if (0 == strcmp(optString, "intersect-min"))          opName = opIntersectMin;
    else if (0 == strcmp(optString, "intersect-max"))          opName = opIntersectMax;
    else if (0 == strcmp(optString, "intersect-sum"))          opName = opIntersectSum;
    else if (0 == strcmp(optString, "difference"))             opName = opDifference;
    else if (0 == strcmp(optString, "symmetric-difference"))   opName = opSymmetricDifference;
    else if (0 == strcmp(optString, "print"))                  opName = opPrint;

    else if (0 == strcmp(optString, "output"))            //  Flag the next arg as the output name
      outputArg = arg + 1;                                //  if we see 'output'.

    else if (arg == outputArg)                            //  If this is the output nane, make a new
      writer = new kmerCountFileWriter(inoutName);        //  output writer.

    else if ((opStack.size() > 0) &&                      //  If a command exists,
             (opStack.top()->isCounting()  == false) &&   //  and it isn't for counting,
             (AS_UTL_fileExists(indexName) == true))      //  and the meryl index file exists,
      reader = new kmerCountFileReader(inoutName);        //  add a kmerCountFile as input to the current command.

    else if ((opStack.size() > 0) &&                      //  If a command exists,
             (opStack.top()->isCounting()  == true) &&    //  and it IS for counting,
             (AS_UTL_fileExists(inoutName)  == true))     //  and the file exists,
      sequence = new dnaSeqFile(inoutName);               //  add a sequence file as input to the current command.

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Don't know what to do with '%s'.", optString);
      err.push_back(s);
    }

    //  Now, do something with the parsed word.
    //    If 'op' is set, make a new command.
    //    If 'writer' exists, set the output of the top most command to that.
    //    If 'reader' or 'sequence' exist, add it to the inputs of the top most command.

    if (opName != opNothing) {
      merylOperation *newOp = new merylOperation(opName, allowedThreads, allowedMemory);

      if (opStack.empty() == false)        //  If a command exists, the new command
        opStack.top()->addInput(newOp);    //  supplies input to the existing command.

      opStack.push(newOp);                 //  Make the new command the current command.
      opName = opNothing;
    }

    if ((writer != NULL) &&                //  Add the writer to the top most command, if
        (opStack.size() > 0)) {            //  one actually exists.  If not, wait until the
      opStack.top()->addOutput(writer);    //  command is created.
      writer = NULL;
    }

    if (reader != NULL) {                  //  Add the reader to the top most command.  The
      opStack.top()->addInput(reader);     //  top most command always exists (else we'd error out
      reader = NULL;                       //  when creating the reader object above).
    }

    if (sequence != NULL) {                //  Same story, different object.
      opStack.top()->addInput(sequence);
      sequence = NULL;
    }

    //  Finally, if we've been told to terminate the command, do so.

    for (; terminating > 0; terminating--)
      opStack.pop();
  }

  //  If any errors, fail.

  if ((argc == 1) ||        //  No commands
      (err.size() > 0)) {   //  Errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  k=<K>        work with mers of size K bases\n");
    fprintf(stderr, "  n=<N>        expect N mers in the input\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m <M>       use no more than M gigabytes of memory\n");
    fprintf(stderr, "  -t <T>       use no more than T compute threads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  A meryl command line is formed as a series of commands and files, possibly\n");
    fprintf(stderr, "  grouped using square brackets.  Each command operates on the file(s) that\n");
    fprintf(stderr, "  are listed after it.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  COMMANDS:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    count\n");
    fprintf(stderr, "    count-forward\n");
    fprintf(stderr, "    count-reverse\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    less-than\n");
    fprintf(stderr, "    greater-than\n");
    fprintf(stderr, "    equal-to\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    union\n");
    fprintf(stderr, "    union-min\n");
    fprintf(stderr, "    union-max\n");
    fprintf(stderr, "    union-sum\n");
    fprintf(stderr, "    intersect\n");
    fprintf(stderr, "    intersect-min\n");
    fprintf(stderr, "    intersect-max\n");
    fprintf(stderr, "    intersect-sum\n");
    fprintf(stderr, "    difference\n");
    fprintf(stderr, "    symmetric-difference\n");
    fprintf(stderr, "    print\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  EXAMPLES:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    meryl -k 16 count input.fasta output counted-kmers \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      Count the k=16 mers in 'input.fasta' and write results to the meryl directory 'counted-kmers'.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  //  Pop the stack until we get back to the root operation.

  while (opStack.size() > 1)
    opStack.pop();

  merylOperation *op = opStack.top();

  //  Dump the tree that we're going to process.

  fprintf(stderr, "Detected %u available threads and %.3f GB memory.\n",
          physThreads, physMemory / 1024.0 / 1024.0 / 1024.0);

  //  Then keep calling nextMer() on the root node until all kmers have been processed.

  while (op->nextMer(true) == true)
    ;

  //  Cleanup by deleting the root operation, which will then delete all its inputs, etc.

  delete op;

  fprintf(stderr, "Bye.\n");

  return(0);
}
