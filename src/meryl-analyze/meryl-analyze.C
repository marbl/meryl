
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
#include "sequence.H"
#include "bits.H"




int
main(int argc, char **argv) {
  char   *inputDBname = NULL;
  bool    verbose     = false;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-mers") == 0) {
      inputDBname = argv[++arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputDBname == NULL)
    err.push_back("No query meryl database (-mers) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }



  fprintf(stderr, "Open meryl database '%s'.\n", inputDBname);
  merylFileReader   *merylDB = new merylFileReader(inputDBname);
  uint64             nKmers  = 0;

  char               fstr[65];

  uint32             maxCount = 16 * 1024 * 1024;
  uint64            *AGhist[65];
  uint64            *TChist[65];

  for (uint32 ii=0; ii<65; ii++) {
    AGhist[ii] = new uint64 [maxCount];
    TChist[ii] = new uint64 [maxCount];

    memset(AGhist[ii], 0, sizeof(uint64) * maxCount);
    memset(TChist[ii], 0, sizeof(uint64) * maxCount);
  }


  while (merylDB->nextMer() == true) {
    uint32  value = merylDB->theValue();
    kmer    fmer  = merylDB->theFMer();
    kmdata  fbits = fmer;

    uint32  fscore = 0,  fa = 0, fg = 0;
    uint32  rscore = 0,  rt = 0, rc = 0;

    for (uint32 ii=0; ii<kmer::merSize(); ii++) {
      kmdata fbase = fbits & 0x03;

      switch (fbase) {
        case 0x00:  //  A
          if ((rt > 0) && (rc > 0))   rscore += rt + rc;
          rt = rc = 0;

          fa++;
          break;
        case 0x01:  //  C
          rc++;

          if ((fa > 0) && (fg > 0))   fscore += fa + fg;
          fa = fg = 0;
          break;
        case 0x02:  //  T
          rt++;

          if ((fa > 0) && (fg > 0))   fscore += fa + fg;
          fa = fg = 0;
          break;
        case 0x03:  //  G
          if ((rt > 0) && (rc > 0))   rscore += rt + rc;
          rt = rc = 0;

          fg++;
          break;
      }

      fbits >>= 2;
    }

    if ((fa > 0) && (fg > 0))   fscore += fa + fg;
    if ((rt > 0) && (rc > 0))   rscore += rt + rc;

    if (verbose)
      fprintf(stderr, "%s  %8u  AG= %2u TC= %2u\n",
              fmer.toString(fstr), value,
              fscore, rscore);

    if (fscore < maxCount) {
      AGhist[fscore][value]++;
    }

    if (rscore < maxCount) {
      TChist[rscore][value]++;
    }

    if ((++nKmers % 100000) == 0)
      fprintf(stderr, "Loaded %li kmers.\n", nKmers);
  }

  delete merylDB;



  for (uint32 ll=0; ll<65; ll++) {
    char    outName[FILENAME_MAX+1];

    sprintf(outName, "hist-runlen=%02u-ag", ll);

    FILE *F = AS_UTL_openOutputFile(outName);

    for (uint32 cc=0; cc<maxCount; cc++) {
      if (AGhist[ll][cc] > 0)
        fprintf(F, "%2u\t%9u\t%lu\n", ll, cc, AGhist[ll][cc]);
    }

    fclose(F);
  }




  for (uint32 ll=0; ll<65; ll++) {
    char    outName[FILENAME_MAX+1];

    sprintf(outName, "hist-runlen=%02u-tc", ll);

    FILE *F = AS_UTL_openOutputFile(outName);

    for (uint32 cc=0; cc<maxCount; cc++) {
      if (TChist[ll][cc] > 0)
        fprintf(F, "%2u\t%9u\t%lu\n", ll, cc, TChist[ll][cc]);
    }

    fclose(F);
  }





  exit(0);
}
