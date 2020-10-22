
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

#include "runtime.H"

#include "kmers.H"
#include "system.H"
#include "sequence.H"
#include "bits.H"


#define OP_NONE       0
#define OP_DUMP       1
#define OP_EXISTENCE  2
#define OP_INCLUDE    3
#define OP_EXCLUDE    4



void
dumpExistence(dnaSeqFile                  *sfile,
              compressedFileWriter        *ofile,
              vector<merylExactLookup *>  &klookup,
              vector<const char *>        &klabel) {

  //  Build a list of labels for each database.  If no labels are provided,
  //  this is just an empty string.

  char   **labels = new char * [klookup.size()];

  for (uint32 ll=0; ll<klookup.size(); ll++) {

    //  If we don't have the ll'th input label, make an empty string.

    if (klabel.size() <= ll) {
      labels[ll]    = new char [1];
      labels[ll][0] = 0;
      continue;
    }

    //  Otherwise, we have a label, so allocate space for a tab, a copy of
    //  the label, and a NUL byte, then create the string we'll output.

    labels[ll] = new char [strlen(klabel[ll]) + 2];

    labels[ll][0] = '\t';
    strcpy(labels[ll] + 1, klabel[ll]);
  }

  //  Scan each sequence against each database.

  char     fString[65];
  char     rString[65];
  dnaSeq   seq;

  for (uint32 seqId=0; sfile->loadSequence(seq); seqId++) {
    kmerIterator  kiter(seq.bases(), seq.length());

    while (kiter.nextBase()) {
      if (kiter.isValid() == false) {
        fprintf(ofile->file(), "%s\t%u\t%lu\t%c\n",
                seq.name(),
                seqId,
                kiter.position(),
                kiter.isACGTbgn() ? 'n' : 'N');
      }

      else {
        for (uint32 dd=0; dd<klookup.size(); dd++) {
          uint64  fValue = 0;
          uint64  rValue = 0;
          bool    fExists = klookup[dd]->exists(kiter.fmer(), fValue);
          bool    rExists = klookup[dd]->exists(kiter.rmer(), rValue);

          fprintf(ofile->file(), "%s\t%u\t%lu\t%c\t%s\t%lu\t%s\t%lu\t%s\n",
                  seq.name(),
                  seqId,
                  kiter.position(),
                  (fExists || rExists) ? 'T' : 'F',
                  kiter.fmer().toString(fString), fValue,
                  kiter.rmer().toString(rString), rValue,
                  labels[dd]);
        }
      }
    }
  }
}



void
reportExistence(dnaSeqFile                  *sfile,
                compressedFileWriter        *ofile,
                vector<merylExactLookup *>  &klookup,
                vector<const char *>        &klabel) {
  dnaSeq   seq;

  while (sfile->loadSequence(seq)) {
    kmerIterator  kiter(seq.bases(), seq.length());

    uint64   nKmer      = 0;
    uint64   nKmerFound = 0;

    while (kiter.nextMer()) {
      nKmer++;

      if ((klookup[0]->value(kiter.fmer()) > 0) ||
          (klookup[0]->value(kiter.rmer()) > 0))
        nKmerFound++;
    }

    fprintf(ofile->file(), "%s\t%lu\t%lu\t%lu\n", seq.name(), nKmer, klookup[0]->nKmers(), nKmerFound);
  }
}



void
filter(dnaSeqFile                  *sfile1,
       dnaSeqFile                  *sfile2,
       compressedFileWriter        *ofile1,
       compressedFileWriter        *ofile2,
       vector<merylExactLookup *>  &klookup,
       bool                         outputIfFound) {

  //  Do nothing if there are no sequences.

  if ((sfile1 == NULL) && (sfile2 == NULL))
    return;

  //  While we load sequences from all files supplied...

  dnaSeq  seq1;
  dnaSeq  seq2;

  uint64   nReads      = 0;
  uint64   nReadsFound = 0;

  while (((sfile1 == NULL) || (sfile1->loadSequence(seq1))) &&
         ((sfile2 == NULL) || (sfile2->loadSequence(seq2)))) {
    uint32 nKmerFound = 0;

    nReads++;

    if (seq1.length() > 0) {
      kmerIterator  kiter(seq1.bases(), seq1.length());

      while (kiter.nextMer())
        if ((klookup[0]->value(kiter.fmer()) > 0) ||
            (klookup[0]->value(kiter.rmer()) > 0))
          nKmerFound++;
    }

    if (seq2.length() > 0) {
      kmerIterator  kiter(seq2.bases(), seq2.length());

      while (kiter.nextMer())
        if ((klookup[0]->value(kiter.fmer()) > 0) ||
            (klookup[0]->value(kiter.rmer()) > 0))
          nKmerFound++;
    }

    //  Report the sequence if:
    //    any kmers are found and     ifFound
    //    no  kmers are found and not ifFound

    if ((nKmerFound > 0) == outputIfFound) {
      nReadsFound++;

      if (sfile1 != NULL) {
        if (seq1.quals()[0] == 0)   fprintf(ofile1->file(), ">%s nKmers=%u\n%s\n",        seq1.name(), nKmerFound, seq1.bases());
        else                        fprintf(ofile1->file(), "@%s nKmers=%u\n%s\n+\n%s\n", seq1.name(), nKmerFound, seq1.bases(), seq1.quals());
      }

      if (sfile2 != NULL) {
        if (seq2.quals()[0] == 0)   fprintf(ofile2->file(), ">%s nKmers=%u\n%s\n",        seq2.name(), nKmerFound, seq2.bases());
        else                        fprintf(ofile2->file(), "@%s nKmers=%u\n%s\n+\n%s\n", seq2.name(), nKmerFound, seq2.bases(), seq2.quals());
      }
    }
  }

  fprintf(stderr, "\nIncluding %lu reads (or read pairs) out of %lu.\n", nReadsFound, nReads);
}



int
main(int argc, char **argv) {
  char const     *seqName1 = nullptr;
  char const     *seqName2 = nullptr;

  char const     *outName1 = "-";
  char const     *outName2 = nullptr;

  vector<const char *>  inputDBname;
  vector<const char *>  inputDBlabel;

  uint64          minV       = 0;
  uint64          maxV       = UINT64_MAX;
  uint32          threads    = getMaxThreadsAllowed();
  double          memory     = getMaxMemoryAllowed() / 1024.0 / 1024.0 / 1024.0;
  uint32          reportType = OP_NONE;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      seqName1 = argv[++arg];

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        seqName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-mers") == 0) {
      while ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        inputDBname.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-labels") == 0) {
      while ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        inputDBlabel.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-output") == 0) {
      outName1 = argv[++arg];

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        outName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      minV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      maxV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      threads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memory = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-dump") == 0) {
      reportType = OP_DUMP;

    } else if (strcmp(argv[arg], "-existence") == 0) {
      reportType = OP_EXISTENCE;

    } else if (strcmp(argv[arg], "-include") == 0) {
      reportType = OP_INCLUDE;

    } else if (strcmp(argv[arg], "-exclude") == 0) {
      reportType = OP_EXCLUDE;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((seqName1 == NULL) && (seqName2 == NULL))
    err.push_back("No input sequences (-sequence) supplied.\n");
  if (inputDBname.size() == 0)
    err.push_back("No query meryl database (-mers) supplied.\n");
  if (reportType == OP_NONE)
    err.push_back("No report-type (-existence, etc) supplied.\n");
  if ((seqName1 != nullptr) && (outName1 == nullptr))
    err.push_back("No output file (-output) supplied.\n");
  if ((seqName2 != nullptr) && (outName2 == nullptr))
    err.push_back("No second output file (-output) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s <report-type> \\\n", argv[0]);
    fprintf(stderr, "         -sequence <input1.fasta> [<input2.fasta>] \\\n");
    fprintf(stderr, "         -output   <output1>      [<output2>]\n");
    fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] \\\n");
    fprintf(stderr, "         -labels   <input1name>   [<input2name>]   [...]\n");
    fprintf(stderr, "  Query the kmers in meryl database(s) <input.meryl> with the sequences\n");
    fprintf(stderr, "  in <input.fasta>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Multiple databases are supported.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Up to two inptu sequences are supported (only for -include / -exclude).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Input files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Output from each input is sent to the associated output file.  Files will be\n");
    fprintf(stderr, "  compressed if the appropriate extension is supplied (gz, bz2 or xz).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Each input database can be filtered by value.  More advanced filtering\n");
    fprintf(stderr, "  requires a new database to be constructed using meryl.\n");
    fprintf(stderr, "    -min   m    Ignore kmers with value below m\n");
    fprintf(stderr, "    -max   m    Ignore kmers with value above m\n");
    fprintf(stderr, "    -threads t  Number of threads to use when constructing lookup table.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory m   Don't use more than m GB memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Exactly one report type must be specified.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -existence");
    fprintf(stderr, "    Report a tab-delimited line for each sequence showing the number of kmers\n");
    fprintf(stderr, "    in the sequence, in the database, and in both.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output:  seqName <tab> mersInSeq <tab> mersInDB <tab> mersInBoth\n");
    fprintf(stderr, "      seqName    - name of the sequence\n");
    fprintf(stderr, "      mersInSeq  - number of mers in the sequence\n");
    fprintf(stderr, "      mersInDB   - number of mers in the meryl database\n");
    fprintf(stderr, "      mersInBoth - number of mers in the sequence that are\n");
    fprintf(stderr, "                   also in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -dump\n");
    fprintf(stderr, "    Report a tab-delimited line reporting each kmer in the input sequences, in\n");
    fprintf(stderr, "    order, annotated with the value of the kmer in the input database.  If the kmer\n");
    fprintf(stderr, "    does not exist in the database its value will be reported as zero.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output:  seqName <tab> seqId <tab> seqPos <tab> exists <tab> fwd-mer <tab> fwd-val <tab> rev-mer <tab> rev-val\n");
    fprintf(stderr, "      seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "      seqId      - numeric version of the seqName (0-based)\n");
    fprintf(stderr, "      seqPos     - start position (0-based) of the kmer in the sequence\n");
    fprintf(stderr, "      exists     - 'T' if the kmer exists in the database, 'F' if it does not\n");
    fprintf(stderr, "      fwd-mer    - forward mer sequence\n");
    fprintf(stderr, "      fwd-val    - value of the forward mer in the database\n");
    fprintf(stderr, "      rev-mer    - reverse mer sequence\n");
    fprintf(stderr, "      rev-val    - value of the reverse mer in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -include / -exclude\n");
    fprintf(stderr, "    Extract sequences containing (-include) or not containing (-exclude) kmers in\n");
    fprintf(stderr, "    any input database.  Output sequences are written in the same format as the input\n");
    fprintf(stderr, "    sequences, with the number of kmers found added to the name.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    If two input files are supplied, the corresponding sequences are treated as a pair.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output:  sequence given format (fasta or fastq) with the number of overlapping kmers appended\n");
    fprintf(stderr, "             if pairs of sequences are given, R1 will be stdout and R2 be named as <output.r2>\n");
    fprintf(stderr, "              <output.r2> will be automatically compressed if ends with .gz, .bz2, or xs\n");
    fprintf(stderr, "      seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "      mersInBoth - number of mers in both sequence and in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -exclude       Extract sequences *NOT containing* kmers in <input.meryl>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "     output:  sequence given format (fasta or fastq) without reads containing kmers\n");
    fprintf(stderr, "              if pairs of sequences are given, R1 will be stdout and R2 be named as <output.r2>\n");
    fprintf(stderr, "              <output.r2> will be automatically compressed if ends with .gz, .bz2, or xs\n");
    fprintf(stderr, "         seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(threads);

  //  Open the kmers, build a lookup table.

  vector<merylFileReader *>   merylDBs;
  vector<merylExactLookup *>  kmerLookups;
  vector<double>              minMem;
  vector<double>              optMem;

  double                      minMemTotal = 0.0;
  double                      optMemTotal = 0.0;

  bool                        useMin = false;
  bool                        useOpt = false;

  for (uint32 ii=0; ii<inputDBname.size(); ii++) {
    merylFileReader   *merylDB    = new merylFileReader(inputDBname[ii]);
    merylExactLookup  *kmerLookup = new merylExactLookup();

    merylDBs   .push_back(merylDB);
    kmerLookups.push_back(kmerLookup);
    minMem     .push_back(0.0);
    optMem     .push_back(0.0);
  }

  for (uint32 ii=0; ii<inputDBname.size(); ii++) {
    double  minMem = 0.0;
    double  optMem = 0.0;

    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Estimating memory usage for '%s'.\n", inputDBname[ii]);
    fprintf(stderr, "--\n");

    kmerLookups[ii]->estimateMemoryUsage(merylDBs[ii], memory, minMem, optMem, minV, maxV);

    minMemTotal += minMem;
    optMemTotal += optMem;
  }

  if      (optMemTotal <= memory)
    useOpt = true;
  else if (minMemTotal <= memory)
    useMin = true;

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Minimal memory needed: %.3f GB%s\n", minMemTotal, (useMin) ? "  enabled" : "");
  fprintf(stderr, "-- Optimal memory needed: %.3f GB%s\n", optMemTotal, (useOpt) ? "  enabled" : "");
  fprintf(stderr, "-- Memory limit           %.3f GB\n",   memory);
  fprintf(stderr, "--\n");

  if ((useMin == false) &&
      (useOpt == false)) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not enough memory to load databases.  Increase -memory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  for (uint32 ii=0; ii<inputDBname.size(); ii++) {
    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", inputDBname[ii]);
    fprintf(stderr, "--\n");

    if (kmerLookups[ii]->load(merylDBs[ii], memory, useMin, useOpt, minV, maxV) == false)
      exit(1);
  }


  for (uint32 ii=0; ii<inputDBname.size(); ii++) {
    delete merylDBs[ii];
  }

  //  Open input sequences.

  dnaSeqFile  *seqFile1 = NULL;
  dnaSeqFile  *seqFile2 = NULL;

  if (seqName1 != NULL) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName1);

    seqFile1 = new dnaSeqFile(seqName1);
  }

  if (seqName2 != NULL) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName2);

    seqFile2 = new dnaSeqFile(seqName2);
  }

  //  Open output writers.

  compressedFileWriter  *outFile1 = (outName1 == nullptr) ? nullptr : new compressedFileWriter(outName1);
  compressedFileWriter  *outFile2 = (outName2 == nullptr) ? nullptr : new compressedFileWriter(outName2);

  //  Do something.

  if (reportType == OP_DUMP)
    dumpExistence(seqFile1, outFile1, kmerLookups, inputDBlabel);

  if (reportType == OP_EXISTENCE)
    reportExistence(seqFile1, outFile1, kmerLookups, inputDBlabel);

  if (reportType == OP_INCLUDE)
    filter(seqFile1, seqFile2, outFile1, outFile2, kmerLookups, true);

  if (reportType == OP_EXCLUDE)
    filter(seqFile1, seqFile2, outFile1, outFile2, kmerLookups, false);

  //  Done!

  delete seqFile1;
  delete seqFile2;

  delete outFile1;
  delete outFile2;

  for (uint32 ii=0; ii<kmerLookups.size(); ii++)
    delete kmerLookups[ii];

  fprintf(stderr, "Bye!\n");

  exit(0);
}
