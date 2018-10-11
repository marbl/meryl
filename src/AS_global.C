
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
 *    Brian P. Walenz from 2007-AUG-03 to 2013-AUG-01
 *      are Copyright 2007-2009,2011-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2009-MAR-06
 *      are Copyright 2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-MAR-03
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-08
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "meryl_version.H"

#include "files.H"
#include "system.H"

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#ifdef _GLIBCXX_PARALLEL
#include <parallel/algorithm>
#include <parallel/settings.h>
#endif


//  We take argc and argv, so, maybe, eventually, we'll want to parse
//  something out of there.  We return argc in case what we parse we
//  want to remove.
//
int
AS_configure(int argc, char **argv) {


#ifdef X86_GCC_LINUX
  //  Set the x86 FPU control word to force double precision rounding
  //  rather than `extended' precision rounding. This causes base
  //  calls and quality values on x86 GCC-Linux (tested on RedHat
  //  Linux) machines to be identical to those on IEEE conforming UNIX
  //  machines.
  //
  fpu_control_t fpu_cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( fpu_cw );
#endif


#ifdef _GLIBCXX_PARALLEL_SETTINGS_H
  __gnu_parallel::_Settings s = __gnu_parallel::_Settings::get();

  //  Force all algorithms to be parallel.
  //  Force some algs to be sequential by using a tag, eg:
  //    sort(a, a+end, __gnu_parallel::sequential_tag());
  //
  //s.algorithm_strategy = __gnu_parallel::force_parallel;

  //  The default seems to be 1000, way too small for us.
  s.sort_minimal_n = 128 * 1024;

  //  The default is MWMS, which, at least on FreeBSD 8.2 w/gcc46, is NOT inplace.
  //  Then again, the others also appear to be NOT inplace as well.
  //s.sort_algorithm = __gnu_parallel::MWMS;
  //s.sort_algorithm = __gnu_parallel::QS_BALANCED;
  //s.sort_algorithm = __gnu_parallel::QS;

  __gnu_parallel::_Settings::set(s);
#endif


  //  Default to one thread.  This is mostly to disable the parallel sort,
  //  which seems to have a few bugs left in it.  e.g., a crash when using 48
  //  threads, but not when using 47, 49 or 64 threads.

  omp_set_num_threads(1);


  //  Install a signal handler to catch seg faults and errors.

  AS_UTL_installCrashCatcher(argv[0]);


  //  Set the start time.

  getProcessTime();


  //
  //  Et cetera.
  //

  for (int32 i=0; i<argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fputs(MERYL_VERSION, stderr);
      exit(0);
    }
  }

  return(argc);
}
