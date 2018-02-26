
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

#ifndef LIBSEQUENCE_H
#define LIBSEQUENCE_H

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "readBuffer.H"

//  Saves the file offset of the first byte in the record:
//    for FASTA, the '>'
//    for FASTQ, the '@'.

class dnaSeqIndexEntry {
public:
  dnaSeqIndexEntry() {
    _fileOffset     = UINT64_MAX;
    _sequenceLength = 0;
  };
  ~dnaSeqIndexEntry() {
  };

  uint64   _fileOffset;
  uint64   _sequenceLength;
};




class dnaSeqFile {
public:
  dnaSeqFile(const char *filename, bool indexed=false);
  ~dnaSeqFile();

  compressedFileReader  *_file;
  readBuffer            *_buffer;

  dnaSeqIndexEntry      *_index;
  uint64                 _indexLen;
  uint64                 _indexMax;

  //  If indexed, searches the index for the proper sequence.
  //
  //  If not indexed, searches forward in the file for the sequence.  If not found,
  //  the file will be at the end.
  //
  //  In both cases, the file is left positioned at the start of the sequence header.
  //
  //  Returns true if found, false if not.
  //
  bool     findSequence(uint64 i);
  bool     findSequence(const char *name);

private:
  uint64
  loadFASTA(char   *&name,     uint32   nameMax,
            char   *&seq,
            uint8  *&qlt,      uint64   seqMax);

  uint64
  loadFASTQ(char   *&name,     uint32   nameMax,
            char   *&seq,
            uint8  *&qlt,      uint64   seqMax);


public:
  //  Return the next sequence in the file.
  //  Returns the number of bases loaded; zero if EOF.
  //
  bool   loadSequence(char   *&name,     uint32   nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64   seqMax,
                      uint64  &seqLen);

  //  Returns a chunk of sequence from the file, up to 'length' bases or
  //  the end of the current sequence.
  //
  //  Returns the number of bases loaded; zero if EOF.
  //
  bool   loadBases(char    *seq,
                   uint64   maxLength,
                   uint64  &seqLength);
};


#endif  //  LIBSEQUENCE_H
