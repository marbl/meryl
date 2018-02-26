
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

#include "AS_global.H"

#include "mt19937ar.H"

#include <vector>
#include <algorithm>

using namespace std;

//  Split words into two pieces, a (32-N) bit prefix and a (N) bit suffix.
//  The prefix is encoded as offsets to the next.
//  The suffix is encoded directly as N bits.
//
//  Each word is then encoded as the difference between then last prefix
//  and this one (unary encoded), and the suffix (binary encoded).
//
//  ppss
//  00aa -> 1aa
//  00bb -> 1bb
//  01cc -> 01cc
//  01dd -> 1dd
//  05ee -> 00001ee
//
//  The total size of the encoding will then be:
//    2^p         bits (size of prefix)
//    N           bits (end of unary encoding)
//    N * suffix  bits
//
//  There's a little space optimization if the last prefix used
//  is less than the maximum prefix.  This is ignored.

int32
main(int32 argc, char **argv) {
  uint64  width  = strtouint64(argv[1]);
  uint64  values = strtouint64(argv[2]);

  for (uint64 N=1; N<width; N++) {
    uint64  eBits = 0;                //  encoded bits
    uint64  rBits = width * values;   //  raw bits

    //  If N = width-1, there are 2^1 possible prefix, need 2^1-1 bits to store this
    //     N = width-2, there are 2^2 possible prefix, need 2^2-1 bits to store this

    eBits += ((uint64)1 << (width - N)) - 1;
    eBits += 1 * values;
    eBits += N * values;

    fprintf(stderr, "%2lu/%2lu encoded %13lu raw %13lu compression %6.4f\n", width-N, N, eBits, rBits, 1.0 - (double)eBits / rBits);
  }

  exit(0);
}
