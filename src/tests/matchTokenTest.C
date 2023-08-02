
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

#include "matchToken.H"

using namespace merylutil;
using namespace merylutil::kmers::v1;


int
main(int argc, char **argv) {
  char const *optout = nullptr;

  assert(matchToken("nope:",   optout, "output:") == false);
  assert(matchToken("no:",     optout, "output:") == false);
  assert(matchToken("n:",      optout, "output:") == false);
  assert(matchToken("outpuy:", optout, "output:") == false);

  assert((matchToken("output:", optout, "output:") == true) && (*optout == 0));
  assert((matchToken("out:",    optout, "output:") == true) && (*optout == 0));
  assert((matchToken("o:",      optout, "output:") == true) && (*optout == 0));

  assert((matchToken("output:database", optout, "output:") == true) && (strcmp(optout, "database") == 0));
  assert((matchToken("out:database",    optout, "output:") == true) && (strcmp(optout, "database") == 0));
  assert((matchToken("o:database",      optout, "output:") == true) && (strcmp(optout, "database") == 0));

  assert((matchToken("output:database", optout, "output:", true) == true) && (strcmp(optout, "database") == 0));
  assert((matchToken("out:database",    optout, "output:", true) == false));
  assert((matchToken("o:database",      optout, "output:", true) == false));

  assert((matchToken("database",        optout, "database") == true));
  assert((matchToken("database",        optout, "database") == true) && (*optout == 0));

  assert(matchToken("database:", optout, "database:") == true);
  assert(matchToken("database:", optout, "database")  == false);
  assert(matchToken("database",  optout, "database:") == false);
  assert(matchToken("database",  optout, "database")  == true);

  assert(matchToken("database=", optout, "database=") == true);
  assert(matchToken("database=", optout, "database")  == false);
  assert(matchToken("database",  optout, "database=") == false);
  assert(matchToken("database",  optout, "database")  == true);

  assert(matchToken("database:", optout, "database=") == false);
  assert(matchToken("database=", optout, "database:") == false);


  assert(matchToken("database:", optout, "database:", true) == true);
  assert(matchToken("database:", optout, "database",  true) == false);
  assert(matchToken("database",  optout, "database:", true) == false);
  assert(matchToken("database",  optout, "database",  true) == true);

  assert(matchToken("database=", optout, "database=", true) == true);
  assert(matchToken("database=", optout, "database",  true) == false);
  assert(matchToken("database",  optout, "database=", true) == false);
  assert(matchToken("database",  optout, "database",  true) == true);

  assert(matchToken("database:", optout, "database=", true) == false);
  assert(matchToken("database=", optout, "database:", true) == false);


  assert(matchToken("data:", optout, "database:", true) == false);
  assert(matchToken("data:", optout, "database",  true) == false);
  assert(matchToken("data",  optout, "database:", true) == false);
  assert(matchToken("data",  optout, "database",  true) == false);


  if (argc > 1) {
    if (matchToken(argv[1], optout, argv[2]) == true)
      fprintf(stderr, "POSITIVE! '%s' matches '%s' with optout '%s'.\n", argv[1], argv[3], optout);
    else
      fprintf(stderr, "NEGATIVE. '%s' doesn't match '%s'.\n", argv[1], argv[3]);
  }

  fprintf(stdout, "matchToken() tests PASS!\n");

  return 0;
}
