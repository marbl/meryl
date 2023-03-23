
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

merylVerbosity  verbosity;



int
main(int argc, char **argv) {
  merylCommandBuilder  *B = new merylCommandBuilder;
  bool                  stopAfterConfigure = false;

  argc = AS_configure(argc, argv);

  uint64  allowedMemory  = getMaxMemoryAllowed();   //  Remember any system-imposed
  uint32  allowedThreads = getMaxThreadsAllowed();  //  limit on cpu/threads.

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {

    //
    //  Scan for debug options and requests for help.
    //

    if (strcmp(argv[arg], "dumpIndex") == 0) {         //  Report the index for the dataset.
      arg++;                                           //  It's just the parameters used for encoding.
      delete new merylFileReader(argv[arg++], true);   //  Expects a meryl db directory as a parameter.
      exit(0);
    }

    if (strcmp(argv[arg], "dumpFile") == 0) {          //  Dump the index for a single data file.
      arg++;                                           //  Expects a meryl file prefix as a parameter.
      dumpMerylDataFile(argv[arg++]);                  //  (e.g., db.meryl/0x000000)
      exit(0);
    }

    if ((strcmp(argv[arg], "-h")     == 0) ||
        (strcmp(argv[arg], "-help")  == 0) ||
        (strcmp(argv[arg], "--help") == 0) ||
        (strcmp(argv[arg], "help")   == 0)) {
      err.push_back(NULL);
      continue;
    }

    //
    //  Handle global options.  These are easier to do outside the
    //  command builder (and besides, they're global options).
    //

    if (strcmp(argv[arg], "-k") == 0) {
      kmerTiny::setSize(strtouint32(argv[++arg]));
      continue;
    }
#ifdef LEGACY_OPTIONS
    if (strncmp(argv[arg], "k=", 2) == 0) {
      fprintf(stderr, "WARNING: obsolete '%s' supplied; use '-k %s' instead.\n",
              argv[arg], argv[arg]+2);
      kmerTiny::setSize(strtouint32(argv[arg]+2));
      continue;
    }
#endif

    if (strcmp(argv[arg], "-l") == 0) {
      kmerTiny::setLabelSize(strtouint32(argv[++arg]));
      continue;
    }

    if ((strcmp(argv[arg],  "-m")      == 0) ||
        (strcmp(argv[arg],  "-memory") == 0) ||
        (strcmp(argv[arg], "--memory") == 0)) {
      allowedMemory = getAllowedMemory(argv[++arg], err);
      continue;
    }
#ifdef LEGACY_OPTIONS
    if (strncmp(argv[arg], "memory=", 7) == 0) {
      fprintf(stderr, "WARNING: obsolete '%s' supplied; use '-m %s' instead.\n",
              argv[arg], argv[arg]+7);
      allowedMemory = getAllowedMemory(argv[arg]+7, err);
      continue;
    }
#endif

    if ((strcmp(argv[arg],  "-t")       == 0) ||
        (strcmp(argv[arg],  "-threads") == 0) ||
        (strcmp(argv[arg], "--threads") == 0)) {
      allowedThreads = getAllowedThreads(argv[++arg], err);
      continue;
    }
#ifdef LEGACY_OPTIONS
    if (strncmp(argv[arg], "threads=", 8) == 0) {
      fprintf(stderr, "WARNING: obsolete '%s' supplied; use '-t %s' instead.\n",
              argv[arg], argv[arg]+8);
      allowedThreads = getAllowedThreads(argv[arg]+8, err);
      continue;
    }
#endif

    if (strncmp(argv[arg], "-V", 2) == 0) {          //  Anything that starts with -V
      for (uint32 vv=1; vv<strlen(argv[arg]); vv++)  //  increases verbosity by the
        verbosity.increaseVerbosity();               //  number of letters.
      continue;
    }

    if (strcmp(argv[arg], "-Q") == 0) {
      verbosity.beQuiet();
      continue;
    }

    if (strcmp(argv[arg], "-P") == 0) {
      verbosity.enableProgressReport();
      continue;
    }

    if (strcmp(argv[arg], "-C") == 0) {
      stopAfterConfigure = true;
      continue;
    }

    //
    //  Throw the option to merylCommandBuilder and let it figure it out.
    //

    B->processWord(argv[arg]);
  }

  //
  //  All done parsing the command line.  Finalize the trees and check for
  //  errors.
  //

  B->buildTrees();

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

  //
  //  Display the action tree.
  //

  if (verbosity.showStandard() == true) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Found %u command tree%s.\n", B->numTrees(), (B->numTrees() == 1) ? "" : "s");

    for (uint32 ii=0; ii<B->numTrees(); ii++)
      B->printTree(B->getTree(ii), 0, 0);
  }

  //
  //  Display any errors in that tree.
  //

  if (B->numErrors() > 0) {         //  Errors from parsing actions.
    fprintf(stderr, "Errors detected in the action tree%s:\n", (B->numTrees() == 1) ? "" : "s");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<B->numErrors(); ii++)
      if (B->getErrors()[ii] != nullptr)
        fprintf(stderr, "  %s\n", B->getErrors()[ii]);

    exit(1);
  }

  //
  //  Stop if we're only configuring.
  //

  if (stopAfterConfigure) {
    fprintf(stderr, "Stopping after configuration.\n");
    return(0);
  }

  //
  //  Counting operations are a big headache.  They don't fit into the
  //  tree nicely:
  //   - they do their own threading, so only one thread can start the operation
  //   - when done, they transform themselves into a pass-through operation that
  //     simply reads the (just counted) input and passes kmers through.
  //
  //  So, we special case them here.  This steps through each action in the tree,
  //  counting kmers, writing to a new database, and finally converting the action
  //  to a null pass-through.
  //

  B->performCounting(allowedMemory, allowedThreads);

  //  Initialize nodes for all the threads.  All the root nodes need to be
  //  initialized before we spawn, so we get thresholds set correctly.

  B->spawnThreads(allowedThreads);

  //  Process each file, in parallel.  Just keep getting the next mer and let
  //  each op do their work.

  omp_set_num_threads(allowedThreads);

  for (uint32 rr=0; rr<B->numTrees(); rr++) {
    merylOpTemplate *tpl = B->getTree(rr);

    //  Actions that were previously a count or a histo/stats on a database
    //  are all done and do not need to process their inputs.
    //
    //  Better - just remove the inputs from them?  But then we need to know
    //  if it supplies input to anything else.
    //
    //if (tpl->_type == merylOpType::opNothing) { //  Was previously a count, histo or stats
    //  continue;                                 //  operation, but it's all done now.
    //}

    if (verbosity.showStandard() == true) {
      fprintf(stderr, "\n");
      fprintf(stderr, "PROCESSING TREE #%u using %u thread%s.\n", rr+1, getMaxThreadsAllowed(), getMaxThreadsAllowed() == 1 ? "" : "s");
    }
    //B->printTree(tpl, 0, 11);

#pragma omp parallel for schedule(dynamic, 1)
    for (uint32 ff=0; ff<64; ff++) {
      merylOpCompute *cpu = B->getTree(rr, ff);

      while (cpu->nextMer() == true)
        ;
    }

    //  Signal that we're done processing.  This will (recursively) collect
    //  any statistics the user has requested be generated.

    tpl->finishAction();

    //  Simply deleting the root template node is enough to delete the entire
    //  tree, including the compute nodes.
    //
    //  NOTE that the pointers in merylCommandBuilder are all dangling.
    //
    delete tpl;
  }

  //  Now that everything is done, delete!

  if (verbosity.showStandard() == true) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Cleaning up.\n");
  }

  delete B;

  if (verbosity.showStandard() == true) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Bye.\n");
  }
  return(0);
}
