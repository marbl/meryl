
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

  uint32   nCols = 170;
  uint32   nRows = 80;

  char **histPlot = new char * [nRows];

  for (uint32 rr=0; rr<nRows; rr++) {
    histPlot[rr] = new char [nCols + 1];

    memset(histPlot[rr], '*', nCols);

    histPlot[rr][nCols] = 0;
  }

  uint64   maxLength = lengths[0];
  uint32   maxCount  = 0;

  uint32  *nSeqPerLen = new uint32 [nCols + 1];


  for (uint32 cc=0; cc<nCols+1; cc++)                   //  Clear the histogram.
    nSeqPerLen[cc] = 0;

  for (uint32 ii=0; ii<nSeqs; ii++)                     //  Count number of sequences per column.
    nSeqPerLen[nCols * lengths[ii] / maxLength]++;

  for (uint32 cc=0; cc<nCols+1; cc++)                   //  Find the maximum count.
    if (maxCount < nSeqPerLen[cc])
      maxCount = nSeqPerLen[cc];

  for (uint32 cc=0; cc<nCols+1; cc++)                   //  Scale column count by number of sequences per row.
    nSeqPerLen[cc] = (uint32)ceil(nSeqPerLen[cc] * nRows / (double)maxCount);

  for (uint32 cc=0; cc<nCols; cc++)                     //  Generate histogram in text.  Err, generate the negative
    for (uint32 rr=0; rr<nRows-nSeqPerLen[cc]; rr++)    //  image of the histogram.
      histPlot[rr][cc] = ' ';

  for (uint32 rr=0; rr<nRows; rr++)
    fprintf(stdout, "%s\n", histPlot[rr]);

  for (uint32 rr=0; rr<nRows; rr++)
    delete [] histPlot[rr];

  delete [] histPlot;
  delete [] nSeqPerLen;

  fprintf(stdout, "maxLength " F_U64 "\n", maxLength);
