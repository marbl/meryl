.. _usage:

====================
COMMAND LINE OPTIONS
====================

Global Options
--------------

Global options apply to every action.

+---------------+---------------------------------------------------------+
| Option        | Description                                             |
+===============+=========================================================+
|``-h``         | show command line usage                                 |
+---------------+                                                         +
|``-help``      |                                                         |
+---------------+                                                         +
|``--help``     |                                                         |
+---------------+                                                         +
|``help``       |                                                         |
+---------------+---------------------------------------------------------+
|``-V``         | increase verbosity of algorithms                        |
+---------------+---------------------------------------------------------+
|``-Q``         | run quietly, don't report any algorithm chatter         |
+---------------+---------------------------------------------------------+
|``-P``         | show progress of algorithms (mostly not implemented)    |
+---------------+---------------------------------------------------------+
|``-C``         | only configure the computation, don't do any processing |
+---------------+---------------------------------------------------------+
|``-k <bases>`` | set the kmer size, necessary for counting actions       |
+---------------+---------------------------------------------------------+
|``-l <bits>``  | set the label size                                      |
+---------------+---------------------------------------------------------+
|``-m``         | set maximum memory usage in gigabytes                   |
+---------------+                                                         +
|``-memory``    |                                                         |
+---------------+                                                         +
|``--memory``   |                                                         |
+---------------+---------------------------------------------------------+
|``-t``         | set the number of CPUs to use                           |
+---------------+                                                         +
|``-threads``   |                                                         |
+---------------+                                                         +
|``--threads``  |                                                         |
+---------------+---------------------------------------------------------+

Obsolete forms of some of these are still allowed.  The kmer size could be
set with `k=<kmer-size>`, memory limits with `memory=<memory-in-gigabytes>`,
thread usage with `threads=<thread-count>`.

Obsolete option `-E` has been removed.  This used to estimate the size of an
imput that could be counted in some given memory size.





Deprecated Options
------------------

  n=

  d=
  distinct=
  f=
  word-frequency=
  t=
  threshold=

  segment=



DEBUG ACTIONS
-------------

dumpIndex <meryl-database>
  Report the parameters used for storing data in this meryl database.

dumpFile <meryl-database-file>
  Dump the raw data from a merylData file in a meryl database.

