
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

//
//  Process each word on the command line:
//
//    Debug options and requests for help.
//     - dumpIndex
//     - dumpFile
//     - -h, -help, --help, help
//
//    Global options.  These are easier to do outside the command builder
//    (and besides, they're global options).
//      -k <kmer-size>
//      -l <label-size>
//      -m <mem-in-gb>    and -memory, --memory
//      -t <thread-limit> and -threads, --threads
//      -V, -VV, etc.
//  
//    Legacy options are k=<kmer-size>, memory=<mem-in-gb> and
//    threads=<thread-limit>.
//
//    If none of the above match, the usual case, toss the word to the
//    command builder and let it take care of it.
//
//  Once the command line has been scanned, finish building the command
//  trees, linking outputs to inputs, computing computable parameters, etc.
//  Then check for errors and display/fail if any are found.
//
//  Counting operations are a big headache.  They don't fit into the
//  tree nicely:
//   - they do their own threading, so only one thread can start the operation
//   - when done, they transform themselves into a pass-through operation that
//     simply reads the (just created) database and passes kmers through.
//
//  So, we special case them here.  This steps through each action in the tree,
//  counting kmers, writing to a new database, and finally converting the action
//  to a null pass-through.
//
//  Once counting is done, the action tree is expanded into 64 copies, one
//  for each database slice, that are run in parallel.  After the run,
//  deleting the commandBuilder will (recursively) delete all the action
//  nodes which will close and open files, etc.
//

int
main(int argc, char **argv) {
  merylCommandBuilder  *B = new merylCommandBuilder;

  argc = globals.initialize(argc, argv);  //  Handles --version, initializes max memory and threads.

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if ((globals.processDebugOption (arg, argv, err) == true) ||
        (globals.processGlobalOption(arg, argv, err) == true) ||
        (globals.processLegacyOption(arg, argv, err) == true))
      ;  //  Do nothing, process() has already done the work,
    else //  Otherwise, let the command builder do the work.
      B->processWord(argv[arg]);
  }

  B->buildTrees();   //  Finalize the processing tree and check for errors.

  if ((argc == 1) ||                //  No commands
      (B->numTrees() == 0) ||       //  No actions
      (err.size() > 0)) {           //  Command line errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);

#include "meryl-usage.H"

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  if (B->displayTreesAndErrors() == true)  //  Halt on errors.
    return 1;

  if (globals.stopAfterConfigure())  //  Print message and return success.
    return fprintf(stderr, "Stopping after configuration.\n"), 0;

  B->performCounting(globals.allowedMemory(), globals.allowedThreads());
  B->spawnThreads(globals.allowedThreads());
  B->runThreads(globals.allowedThreads());

  if (globals.showStandard() == true)
    fprintf(stderr, "\n"
                    "Cleaning up.\n");

  delete B;

  if (globals.showStandard() == true)
    fprintf(stderr, "\n"
                    "Bye.\n");

  return 0;
}
