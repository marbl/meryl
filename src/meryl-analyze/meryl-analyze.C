
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

#include "kmers.H"
#include "sequence.H"
#include "bits.H"

using namespace merylutil;
using namespace merylutil::kmers::v1;

#define OP_NONE   0
#define OP_GA     1
#define OP_GC     2
#define OP_GT     3

//     T - type of the histogram counters
//     V - type of the thing we're counting (must be integral)
template<typename T, typename V>
class denseHistogram {
public:
  denseHistogram(V minValue, V maxValue) {
    _minValue = minValue;
    _maxValue = maxValue;

    _smallestV = _maxValue;
    _largestV  = _minValue;

    _histoLen = 0;
    _histoMax = maxValue - minValue + 1;
    _histo    = new T [_histoMax];
  };

  ~denseHistogram() {
  };

  uint32     minValue(void)    { return(_smallestV); };
  uint32     maxValue(void)    { return(_largestV);  };

  void       insert(V value) {
    if ((_minValue <= value) &&
        (value     <= _maxValue)) {
      _smallestV = std::min(_smallestV, value);
      _largestV  = std::max(_largestV,  value);

      _histo[value - _minValue]++;
    }
  };

  T          report(V value) {
    if ((_minValue <= value) &&
        (value     <= _maxValue))
      return(_histo[value - _minValue]);
    else
      return(0);
  };

private:
  V         _minValue;      //  Minimum value we'll accept into the histogram
  V         _maxValue;

  V         _smallestV;    //  Minimum value we have seen in the input data
  V         _largestV;

  uint32    _histoLen;      //  Allocated histogram.  It doesn't grow.
  uint32    _histoMax;
  T        *_histo;
};


  
template<typename T, typename V>
class sparseHistogram {
public:
  sparseHistogram() {
  };
  sparseHistogram(V minValue, V maxValue) {
    initialize(minValue, maxValue);
  };

  ~sparseHistogram() {
  };

  void       initialize(V minValue, V maxValue) {
    _minValue = minValue;
    _maxValue = maxValue;

    _smallestV = _maxValue;
    _largestV  = _minValue;
  };

  uint32     minValue(void)    { return(_smallestV); };
  uint32     maxValue(void)    { return(_largestV);  };

  void       insert(V value) {
    if ((_minValue <= value) &&
        (value     <= _maxValue)) {
      _smallestV = std::min(_smallestV, value);
      _largestV  = std::max(_largestV,  value);

      _histo[value]++;
    }
  };

  T          report(V value) {
    if ((_minValue <= value) &&
        (value     <= _maxValue) &&
        (_histo.count(value) > 0))
      return(_histo[value]);
    else
      return(0);
  };

private:
  V              _minValue;      //  Minimum value we'll accept into the histogram
  V              _maxValue;

  V              _smallestV;     //  Minimum value we have seen in the input data
  V              _largestV;

  std::map<V,T>  _histo;         //  Histogram data.
};


void
printHist(char* outName, sparseHistogram<uint64,uint32> hist[]) {
  FILE *F = merylutil::openOutputFile(outName);

  for (uint32 ll=0; ll<=kmer::merSize(); ll++) {
    if (hist[ll].minValue() <= hist[ll].maxValue()) {
      for (uint32 cc=hist[ll].minValue(); cc<=hist[ll].maxValue(); cc++) {
        if (hist[ll].report(cc) > 0)
          fprintf(F, "%u\t%u\t%lu\n", ll, cc, hist[ll].report(cc));
      }
    }
  }
  fclose(F);
}

void
histGC(merylFileReader* merylDB,
       char*            outPrefix,
       bool             verbose ) {

  uint64             nKmers  = 0;
  char               fstr[65];
  uint32             maxCount = UINT32_MAX;

  sparseHistogram<uint64,uint32>    GCHist[65];
  sparseHistogram<uint64,uint32>    ATHist[65];

  for (uint32 ii=0; ii<=kmer::merSize(); ii++) {
    GCHist[ii].initialize(0llu, UINT32_MAX);
    ATHist[ii].initialize(0llu, UINT32_MAX);
  }

  while (merylDB->nextMer() == true) {
    uint32  value = merylDB->theValue();
    kmer    fmer  = merylDB->theFMer();
    kmdata  fbits = fmer;

    uint32  fscore = 0, g = 0, c = 0;
    uint32  rscore = 0, a = 0, t = 0;

    for (uint32 ii=0; ii<kmer::merSize(); ii++) {
      kmdata fbase = fbits & 0x03;

      switch (fbase) {
        case 0x01:  //  C
          c++;
          break;
        case 0x03:  //  G
          g++;
          break;
        
        case 0x00:  //  A
          a++;
          break;
        case 0x02:  //  T
          t++;
          break;
      }

      fbits >>= 2;
    }

    fscore = c + g;
    rscore = a + t;

    if (verbose) {
      fprintf(stderr, "%s  %8u  C= %2u G= %2u\n",
              fmer.toString(fstr), value,
              c, g);
      fprintf(stderr, "%s  %8u  A= %2u T= %2u\n",
              fmer.toString(fstr), value,
              a, t);
    }

    if (fscore < maxCount) {
      GCHist[fscore].insert(value);
    }

    if (rscore < maxCount) {
      ATHist[rscore].insert(value);
    }

    if ((++nKmers % 100000000) == 0)
      fprintf(stderr, "Processed %li kmers.\n", nKmers);
  }
  fprintf(stderr, "Processed %li kmers in total.\n\n", nKmers);

  fprintf(stderr, "Output histogram\n");

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.GC.hist", outPrefix);
  printHist(outName, GCHist);

  sprintf(outName, "%s.AT.hist", outPrefix);
  printHist(outName, ATHist);

}

void
histGA(merylFileReader* merylDB,
       char*            outPrefix,
       bool             verbose ) {

  uint64             nKmers  = 0;
  char               fstr[65];
  uint32             maxCount = UINT32_MAX;

  sparseHistogram<uint64,uint32>    CombinedHist[65];
  sparseHistogram<uint64,uint32>    AGhist[65];
  sparseHistogram<uint64,uint32>    TChist[65];

  for (uint32 ii=0; ii<=kmer::merSize(); ii++) {
    CombinedHist[ii].initialize(0llu, UINT32_MAX);
    AGhist[ii].initialize(0llu, UINT32_MAX);
    TChist[ii].initialize(0llu, UINT32_MAX);
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
      AGhist[fscore].insert(value);
    }

    if (rscore < maxCount) {
      TChist[rscore].insert(value);
    }

    if (fscore > rscore) {
      CombinedHist[fscore].insert(value);
    } else {
      CombinedHist[rscore].insert(value);
    }

    if ((++nKmers % 100000000) == 0)
      fprintf(stderr, "Processed %li kmers.\n", nKmers);
  }
  fprintf(stderr, "Processed %li kmers in total.\n\n", nKmers);

  fprintf(stderr, "Output histogram\n");

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.GA_TC.hist", outPrefix);
  printHist(outName, CombinedHist);

  sprintf(outName, "%s.GA.hist", outPrefix);
  printHist(outName, AGhist);

  sprintf(outName, "%s.TC.hist", outPrefix);
  printHist(outName, TChist);

}

void
histGT(merylFileReader* merylDB,
       char*            outPrefix,
       bool             verbose ) {

  uint64             nKmers  = 0;
  char               fstr[65];
  uint32             maxCount = UINT32_MAX;

  sparseHistogram<uint64,uint32>    CombinedHist[65];
  sparseHistogram<uint64,uint32>    GThist[65];
  sparseHistogram<uint64,uint32>    AChist[65];

  for (uint32 ii=0; ii<=kmer::merSize(); ii++) {
    CombinedHist[ii].initialize(0llu, UINT32_MAX);
    GThist[ii].initialize(0llu, UINT32_MAX);
    AChist[ii].initialize(0llu, UINT32_MAX);
  }


  while (merylDB->nextMer() == true) {
    uint32  value = merylDB->theValue();
    kmer    fmer  = merylDB->theFMer();
    kmdata  fbits = fmer;

    uint32  fscore = 0,  fg = 0, ft = 0;
    uint32  rscore = 0,  ra = 0, rc = 0;

    for (uint32 ii=0; ii<kmer::merSize(); ii++) {
      kmdata fbase = fbits & 0x03;

      switch (fbase) {
        case 0x00:  //  A
          if ((ft > 0) && (fg > 0))   fscore += ft + fg;
          ft = fg = 0;

          ra++;
          break;
        case 0x01:  //  C
          rc++;

          if ((ft > 0) && (fg > 0))   fscore += ft + fg;
          ft = fg = 0;
          break;
        case 0x02:  //  T
          ft++;

          if ((ra > 0) && (rc > 0))   rscore += ra + rc;
          ra = rc = 0;
          break;
        case 0x03:  //  G
          if ((ra > 0) && (rc > 0))   rscore += ra + rc;
          ra = rc = 0;

          fg++;
          break;
      }

      fbits >>= 2;
    }

    if ((ft > 0) && (fg > 0))   fscore += ft + fg;
    if ((ra > 0) && (rc > 0))   rscore += ra + rc;

    if (verbose)
      fprintf(stderr, "%s  %8u  GT= %2u AC= %2u\n",
              fmer.toString(fstr), value,
              fscore, rscore);

    if (fscore < maxCount) {
      GThist[fscore].insert(value);
    }

    if (rscore < maxCount) {
      AChist[rscore].insert(value);
    }

    if (fscore > rscore) {
      CombinedHist[fscore].insert(value);
    } else {
      CombinedHist[rscore].insert(value);
    }

    if ((++nKmers % 100000000) == 0)
      fprintf(stderr, "Processed %li kmers.\n", nKmers);
  }
  fprintf(stderr, "Processed %li kmers in total.\n\n", nKmers);

  fprintf(stderr, "Output histogram\n");

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.GT_AC.hist", outPrefix);
  printHist(outName, CombinedHist);

  sprintf(outName, "%s.GT.hist", outPrefix);
  printHist(outName, GThist);

  sprintf(outName, "%s.AC.hist", outPrefix);
  printHist(outName, AChist);

}




int
main(int argc, char **argv) {
  char   *inputDBname = NULL;
  char   *outPrefix   = NULL;
  bool    verbose     = false;
  uint32  reportType  = OP_NONE;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  int                        arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-mers") == 0) {
      inputDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-prefix") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-ga") == 0) {
      reportType = OP_GA;

    } else if (strcmp(argv[arg], "-gc") == 0) {
      reportType = OP_GC;

    } else if (strcmp(argv[arg], "-gt") == 0) {
      reportType = OP_GT;

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
    fprintf(stderr, "usage: %s -mers <meryldb> -prefix <prefix> (-gc | -ga | -gt) \n", argv[0]);
    fprintf(stderr, "  -mers <meryldb>   : meryl database to analyze.\n");
    fprintf(stderr, "  -prefix <prefix>  : prefix for output file(s).\n");
    fprintf(stderr, "  -gc               : generate a histogram of G and C counts. A and T counts will be output as well.\n");
    fprintf(stderr, "  -ga               : generate a histogram of 2-mer microsatellite GA counts. TC counts and combined GA_TC will be output as well.\n");
    fprintf(stderr, "                      combined output reports the GA or TC count, whichever is larger, for a given kmer.\n");
    fprintf(stderr, "  -gt               : generate a histogram of 2-mer microsatellite GT counts. AC counts and combined GT_AC will be output as well.\n");
    fprintf(stderr, "                      combined output reports the GT or AC count, whichever is larger, for a given kmer.\n");
    fprintf(stderr, "  -verbose          : print verbose output.\n");
    fprintf(stderr, "output              : score (out of \"k\" mers) <tab> k-mer multiplicity <tab> count\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  fprintf(stderr, "Open meryl database '%s'.\n", inputDBname);
  merylFileReader   *merylDB = new merylFileReader(inputDBname);

  if (reportType == OP_GC)
    histGC(merylDB, outPrefix, verbose);

  if (reportType == OP_GA)
    histGA(merylDB, outPrefix, verbose);

  if (reportType == OP_GT)
    histGT(merylDB, outPrefix, verbose);

  fprintf(stderr, "Clean up..\n\n");

  delete merylDB;

  fprintf(stderr, "Bye!\n");

  exit(0);
}
