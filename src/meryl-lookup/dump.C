
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
#include "strings.H"
#include <fstream>
#include <iostream>

using namespace std;  //  For ifstream and string.


void
dumpExistence(dnaSeqFile                       *sfile,
              compressedFileWriter             *ofile,
              std::vector<merylExactLookup *>  &klookup,
              std::vector<const char *>        &klabel,
              std::vector<const char *>        &kname) {

  //  Build a list of labels for each database.  If no labels are provided,
  //  this is just an empty string.
  //
  const  uint32 numLookupDB = klookup.size();
  char   **labels = new char * [numLookupDB];

  for (uint32 ll=0; ll<numLookupDB; ll++) {

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

  uint32 ctgn = sfile->numberOfSequences();

  const int threads = omp_get_max_threads();
  uint32 chunks = ctgn;
  
  //  Run chunks at maximum of MAX_FILES. Setting this too high will break network storage
  const int64 openMax = sysconf(_SC_OPEN_MAX) / 2; // prevent too many open files
  if (chunks > openMax) { chunks = openMax;  }
  uint32 chunkLeft = chunks;
  string tmpPrefix = sfile->filename() + string("_tmp_");
  for (uint32 dd = 0; dd < numLookupDB; dd++) {
    splitToWords fp(kname.at(dd), splitPaths);
    tmpPrefix += fp.last() + string("_");
  }

  string tmpFilename;

  fprintf(stderr, "\nTotal of %u sequences found. Will be processed over %d threads, with maximum %u intermediate %s#.dump files\n", ctgn, threads, chunks, tmpPrefix.c_str());

  for (uint32 ii = 0; ii < ctgn; ii += chunks) {
    if ( ii + chunks > ctgn ) chunkLeft = ctgn - ii;
    fprintf(stderr, "Reading sequence %u - %u  ... \n", ii, ii+chunkLeft);

#pragma omp parallel for private(seq, tmpFilename)
    for (uint32 seqId = ii; seqId < ii + chunkLeft; seqId++) {

#pragma omp critical
      {
        //  fprintf(stderr, "Load\t%u\n", seqId);
        sfile->findSequence(seqId);
        sfile->loadSequence(seq);
      }
      //  fprintf(stderr, "Iterate\t%u\n", seqId);
      kmerIterator  kiter(seq.bases(), seq.length());
      tmpFilename = tmpPrefix + to_string(seqId - ii) + ".dump";
      compressedFileWriter  *tmpFile = new compressedFileWriter(tmpFilename.c_str());

      while (kiter.nextBase()) {
        if (kiter.isValid() == false) {
          fprintf(tmpFile->file(), "%s\t%u\t%lu\t%c\n",
              seq.ident(),
              seqId,
              kiter.position(),
              kiter.isACGTbgn() ? 'n' : 'N');
        }

        else {
          for (uint32 dd=0; dd<numLookupDB; dd++) {
            kmvalu  fValue = 0;
            kmvalu  rValue = 0;
            bool    fExists = klookup[dd]->exists(kiter.fmer(), fValue);
            bool    rExists = klookup[dd]->exists(kiter.rmer(), rValue);

            fprintf(tmpFile->file(), "%s\t%u\t%lu\t%c\t%s\t%u\t%s\t%u\t%s\n",
                seq.ident(),
                seqId,
                kiter.position(),
                (fExists || rExists) ? 'T' : 'F',
                kiter.fmer().toString(fString), fValue,
                kiter.rmer().toString(rString), rValue,
                labels[dd]);
          }
        }
      }
      delete tmpFile;
      //  fprintf(stderr, "Done\t%u\n", seqId);

      //  write output, sequence order is unsorted
#pragma omp critical
      {
        //  fprintf(stderr, "Write\t%u\n", seqId);
        ifstream tmpFileI(tmpFilename, std::ios_base::binary);
        if (tmpFileI.is_open()) {
          cout << tmpFileI.rdbuf();
          tmpFileI.close();
          remove(tmpFilename.c_str());
        } else {
          fprintf(stderr, "Unable to open file %s\n", tmpFilename.c_str());
        }
        //  fprintf(stderr, "Done\t%u\n", seqId);
      }
    }

    /*
    fprintf(stderr, "Pushing outputs\n");
    for (uint32 tt=0; tt < chunkLeft; tt++) {
      tmpFilename = tmpPrefix + "tmp." + to_string(tt) + ".dump";
      ifstream tmpFile(tmpFilename, std::ios_base::binary);

      if (tmpFile.is_open()) {
        cout << tmpFile.rdbuf();
        tmpFile.close();
        remove(tmpFilename.c_str());
      } else {
        fprintf(stderr, "Unable to open file %s\n", tmpFilename.c_str());
      }
    }
    */
  }
}

