
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

#include "meryl-lookup.H"


void
lookupGlobal::initialize(void) {

  //  Compute the max length of any single label.  Used during output.

  lookupDBlabelLen = 0;

  for (uint32 ll=0; ll<lookupDBlabel.size(); ll++)
    lookupDBlabelLen = std::max(lookupDBlabelLen, (uint32)strlen(lookupDBlabel[ll]));
}



void
lookupGlobal::loadLookupTables(void) {
  std::vector<merylFileReader *>    merylDBs;    //  Input meryl database.
  std::vector<double>               minMem;      //  Estimated min memory for lookup table.
  std::vector<double>               optMem;      //  Estimated max memory for lookup table.

  //  Open input meryl databases, initialize lookup.

  for (uint32 ii=0; ii<lookupDBname.size(); ii++) {
    merylDBs .push_back(new merylFileReader(lookupDBname[ii]));
    minMem   .push_back(0.0);
    optMem   .push_back(0.0);
    lookupDBs.push_back(new merylExactLookup());
  }

  //  Estimate memory needed for each lookup table.

  double   minMemTotal = 0.0;
  double   optMemTotal = 0.0;

  for (uint32 ii=0; ii<lookupDBname.size(); ii++) {
    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Estimating memory usage for '%s'.\n", lookupDBname[ii]);
    fprintf(stderr, "--\n");

    double  minm, optm;
    lookupDBs[ii]->estimateMemoryUsage(merylDBs[ii], maxMemory, minm, optm, minV, maxV);

    minMemTotal += minm;
    optMemTotal += optm;
  }

  //  Use either the smallest or 'fastest' table, or fail, depending on how
  //  much memory the use lets us use.

  bool  useOpt = (optMemTotal <= maxMemory);
  bool  useMin = (minMemTotal <= maxMemory) && (useOpt == false);

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Minimal memory needed: %.3f GB%s\n", minMemTotal, (useMin) ? "  enabled" : "");
  fprintf(stderr, "-- Optimal memory needed: %.3f GB%s\n", optMemTotal, (useOpt) ? "  enabled" : "");
  fprintf(stderr, "-- Memory limit           %.3f GB\n",   maxMemory);
  fprintf(stderr, "--\n");

  if ((useMin == false) &&
      (useOpt == false)) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not enough memory to load databases.  Increase -memory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (doEstimate == true) {
    fprintf(stderr, "-- Stopping after memory estimated reported; -estimate option enabled.\n");
    exit(0);
  }

  //  Now load the data and forget about the input databases.

  for (uint32 ii=0; ii<lookupDBname.size(); ii++) {
    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", lookupDBname[ii]);
    fprintf(stderr, "--\n");

    if (lookupDBs[ii]->load(merylDBs[ii], maxMemory, useMin, useOpt, minV, maxV) == false)
      exit(1);

    delete merylDBs[ii];
  }
}



//  Open input sequences.
void
lookupGlobal::openInputs(void) {

  if (seqName1) {
    fprintf(stderr, "-- Opening input sequences '%s'.\n", seqName1);
    seqFile1 = openSequenceFile(seqName1);
  }

  if (seqName2) {
    fprintf(stderr, "-- Opening input sequences '%s'.\n", seqName2);
    seqFile2 = openSequenceFile(seqName2);
  }
}



//  Open output writers.
void
lookupGlobal::openOutputs(void) {

  if (outName1) {
    fprintf(stderr, "-- Opening output file '%s'.\n", outName1);
    outFile1 = new compressedFileWriter(outName1);
  }

  if (outName2) {
    fprintf(stderr, "-- Opening output file '%s'.\n", outName1);
    outFile2 = new compressedFileWriter(outName2);
  }
}






int
main(int argc, char **argv) {
  lookupGlobal  *G = new lookupGlobal;

  argc = AS_configure(argc, argv);

  uint32  lThreads = 0;   //  Threads for loading kmer databases
  uint32  nThreads = 0;   //  Threads for computing stuff

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      G->seqName1 = argv[++arg];

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->seqName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-mers") == 0) {
      while ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->lookupDBname.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-labels") == 0) {
      while ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->lookupDBlabel.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-output") == 0) {
      G->outName1 = argv[++arg];

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->outName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      G->minV = (kmvalu)strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      G->maxV = (kmvalu)strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      nThreads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-loadthreads") == 0) {
      lThreads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      G->maxMemory = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-bed") == 0) {
      G->reportType = lookupOp::opBED;

    } else if (strcmp(argv[arg], "-bed-runs") == 0) {
      G->reportType = lookupOp::opBED;
      G->mergeBedRuns   = true;

    } else if (strcmp(argv[arg], "-wig-count") == 0) {
      G->reportType = lookupOp::opWIGcount;

    } else if (strcmp(argv[arg], "-wig-depth") == 0) {
      G->reportType = lookupOp::opWIGdepth;

    } else if (strcmp(argv[arg], "-existence") == 0) {
      G->reportType = lookupOp::opExistence;

    } else if (strcmp(argv[arg], "-include") == 0) {
      G->reportType = lookupOp::opInclude;

    } else if (strcmp(argv[arg], "-exclude") == 0) {
      G->reportType = lookupOp::opExclude;

    } else if (strcmp(argv[arg], "-10x") == 0) {
      G->is10x = true;

    } else if (strcmp(argv[arg], "-estimate") == 0) {
      G->doEstimate = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      G->showProgress = true;

    } else if (strcmp(argv[arg], "-help") == 0) {
      err.push_back(nullptr);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  //  Check for invalid usage.

  G->checkInvalid(err);

  if (err.size() > 0) {
    switch (G->reportType) {
      case lookupOp::opNone:          help(argv[0]);                 break;
      case lookupOp::opBED:           helpBED(argv[0]);              break;
      case lookupOp::opWIGcount:      helpWIGcount(argv[0]);         break;
      case lookupOp::opWIGdepth:      helpWIGdepth(argv[0]);         break;
      case lookupOp::opExistence:     helpExistence(argv[0]);        break;
      case lookupOp::opInclude:       helpIncludeExclude(argv[0]);   break;
      case lookupOp::opExclude:       helpIncludeExclude(argv[0]);   break;
    }

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    return(1);
  }

  if (nThreads == 0)   nThreads = getMaxThreadsAllowed();
  if (lThreads == 0)   lThreads = nThreads;

  omp_set_num_threads(lThreads);   //  Enable threads for loading data.

  G->initialize();
  G->loadLookupTables();
  G->openInputs();
  G->openOutputs();

  omp_set_num_threads(nThreads);   //  Enable threads for computing results.

  switch (G->reportType) {
    case lookupOp::opNone:                                break;
    case lookupOp::opBED:           dumpExistence(G);     break;
    case lookupOp::opWIGcount:      dumpExistence(G);     break;
    case lookupOp::opWIGdepth:      dumpExistence(G);     break;
    case lookupOp::opExistence:     reportExistence(G);   break;
    case lookupOp::opInclude:       filter(G);            break;
    case lookupOp::opExclude:       filter(G);            break;
  }

  delete G;
  fprintf(stderr, "Bye!\n");

  return(0);
}


char const *
makeString(char const *t, char const *l) {
  char   *s = new char [1024];
  sprintf(s, t, l);
  return(s);
}

void
lookupGlobal::checkInvalid(std::vector<char const *> &err) {

  //  If there is no report type, we can skip all the other checks.

  if (reportType == lookupOp::opNone) {
    err.push_back("No report-type (-bed, -wig-count, -wig-depth, -existence, -include, -exclude) supplied.\n");
    return;
  }

  //  Everybody needs at least one input, one database, and one output.

  if (seqName1 == nullptr)
    err.push_back("No input sequences (-sequence) supplied.\n");

  if (lookupDBname.size() == 0)
    err.push_back("No meryl database (-mers) supplied.\n");

  if (outName1 == nullptr)
    err.push_back("No output file (-output) supplied.\n");

  //  Only include/exclude can take a second input.

  if ((reportType != lookupOp::opInclude) &&
      (reportType != lookupOp::opExclude)) {
    if (seqName2 != nullptr)
      err.push_back(makeString("Only one input sequence (-sequence) supported for %s.\n", toString(reportType)));

    if (outName2 != nullptr)
      err.push_back(makeString("Only one output file (-output) supported for %s.\n", toString(reportType)));
  }

  //  If include/exclude, be sure there is an output for the second input,
  //  and that there is only one input database.

  if ((reportType == lookupOp::opInclude) ||
      (reportType == lookupOp::opExclude)) {
    if ((seqName2 != nullptr) &&
        (outName2 == nullptr))
      err.push_back("No second output file (-output) supplied for second input (-input) file.\n");

    if ((seqName2 == nullptr) &&
        (outName2 != nullptr))
      err.push_back("No second input file (-input) supplied for second output (-output) file.\n");

    if (lookupDBname.size()  > 1)
      err.push_back(makeString("Only one meryl database (-mers) supported for %s.\n", toString(reportType)));
  }

  //  Reject labels for things that don't use them.

  if (((reportType == lookupOp::opWIGcount) ||
       (reportType == lookupOp::opWIGdepth) ||
       (reportType == lookupOp::opExistence) ||
       (reportType == lookupOp::opInclude) ||
       (reportType == lookupOp::opExclude)) && (lookupDBlabel.size() > 0))
    err.push_back(makeString("Labels (-labels) not supported for %s.\n", toString(reportType)));
}
