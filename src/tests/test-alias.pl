#!/usr/bin/env perl

use strict;

# usage: test-alias.pl in1.fasta  in2.fasta
#        test-alias.pl in1.import in2.import
#
#  If fasta is supplied, will count kmers using meryl.
#  If import is supplied, will import directly to databases.
#
#  Then will run meryl on each of the 10 aliases and compare results against
#  a naive perl implementation of the operation.

my $meryl  = "/work/meryl-redo/build/bin/meryl";
my $import = "/work/meryl-redo/build/bin/meryl-import";

#my $meryl  = "/work/meryl/build/bin/meryl";
#my $import = "/work/meryl/build/bin/meryl-import";


my $in1 = shift @ARGV;
my $db1 = shift @ARGV;
my @db1;

my $in2 = shift @ARGV;
my $db2 = shift @ARGV;
my @db2;

my @dbM;
my @dbT;

if (!defined($db2)) {
    print STDERR "usage: $0 in1 db1 in2 db2\n";
    print STDERR "  in1 and in2 are either FASTA or meryl-import inputs.\n";
    print STDERR "  They are used to create db1 and db2.\n";
    print STDERR "  All 10 meryl alias operations are tested against the\n";
    print STDERR "  naive implementation in this program.\n";
    exit(1);
}

#  Create and then load databases used for testing.

sub createDatabase ($$$) {
    my $idx = shift @_;
    my $inp = shift @_;
    my $db  = shift @_;

    #  Create a database if one doesn't exist.

    if (! -e $db) {
        if ($inp =~ m/fast/i) {
            #system("$meryl -k 30 count output $db $inp");
            system("$meryl k=30 count output $db $inp");
        } else {
            #system("$import -k 30 -kmers $inp -valuewidth 6 -labelwidth 10 -output $db");
            system("$import -k 30 -kmers $inp -output $db");
        }
    }

    #  Load the database.

    open(F, "$meryl threads=1 print $db 2> /dev/null |");

    if ($idx == 1) {
        undef @db1;
        @db1 = <F>;   chomp @db1;

        #foreach my $i (@db1) {
        #    print "DB1: '$i'\n";
        #}
    }
    else {
        undef @db2;
        @db2 = <F>;   chomp @db2;

        #foreach my $i (@db1) {
        #    print "DB1: '$i'\n";
        #}
    }

    close(F);
}

#  Test!

sub computeValue ($@) {
    my $test = shift @_;
    my $mer  = shift @_;
    my $res;

    my $count = 0;
    my $min   = 9999999;
    my $max   = 0;
    my $sum   = 0;
    my $dif   = 2 * @_[1];   #  Computes 2 * $v1 - $v1 - $v2 - ...
    my $index = undef;
    my $first = undef;
    my $val   = 0;

    while (defined($_[0])) {
        my $i = shift @_;
        my $v = shift @_;
        my $l = shift @_;

        $count += 1;
        $min    = ($min < $v) ? ($min) : ($v);
        $max    = ($max > $v) ? ($max) : ($v);
        $sum   += $v;
        $dif    = ($dif > $v) ? ($dif - $v) : (0);
        $index  = (defined($index)) ? ($index) : ($i);
        $first  = (defined($first)) ? ($first) : ($v);
    }

    if    ($test eq "union") {
        if ($count > 0) {
            $val = $count;
        }
    }
    elsif ($test eq "union-min") {
        if ($count > 0) {
            $val = $min;
        }
    }
    elsif ($test eq "union-max") {
        if ($count > 0) {
            $val = $max;
        }
    }
    elsif ($test eq "union-sum") {
        if ($count > 0) {
            $val = $sum;
        }
    }

    elsif ($test eq "intersect") {
        if ($count == 2) {
            $val = $first;
        }
    }
    elsif ($test eq "intersect-min") {
        if ($count == 2) {
            $val = $min;
        }
    }
    elsif ($test eq "intersect-max") {
        if ($count == 2) {
            $val = $max;
        }
    }
    elsif ($test eq "intersect-sum") {
        if ($count == 2) {
            $val = $sum;
        }
    }

    elsif ($test eq "subtract") {
        if (($index == 1)) {
            $val = $dif;
        }
    }

    elsif ($test eq "difference") {
        if (($index == 1) && ($count == 1)) {
            $val = $sum;
        }
    }

    else {
        die "Unknown test '$test'.\n";
    }

    if ($val > 0) {
        return("$mer\t$val");
    } else {
        return(undef);
    }
}


sub compareKmers ($) {
    my $testName = shift @_;
    my @diff;
    my $dbMlen   = scalar(@dbM);
    my $dbTlen   = scalar(@dbT);

    if ($dbMlen != $dbTlen) {
        push @diff, "number of kmers differs; meryl $dbMlen, perl $dbTlen\n";
    }
    else {
        my $d = 0;

        for (my $ii=0; $ii<$dbMlen; $ii++) {
            my ($km, $vm, $lm) = split '\t', $dbM[$ii];
            my ($kt, $vt, $lt) = split '\t', $dbT[$ii];

            if (($km ne $kt) || ($vm ne $vt) || ($lm ne $lt)) {
                $d++;
            }
        }

        if ($d > 0) {
            push @diff, "$d kmers differ\n";
        }
    }

    return(@diff);
}


sub testMethod ($) {
    my $testName = shift @_;

    print STDERR "Testing $testName\n";

    #  Run meryl.

    undef @dbM;

    #open(F, "$meryl print output test-alias.meryl $testName $db1 $db2 2> /dev/null |");
    open(F, "$meryl threads=1 print $testName $db1 $db2 2> /dev/null |");
    while (<F>) {
        chomp;
        push @dbM, $_;
    }
    close(F);

    #  Load input.

    createDatabase(1, $in1, $db1);
    createDatabase(2, $in2, $db2);

    #  Compute.

    undef @dbT;

    while ((scalar(@db1) > 0) ||
           (scalar(@db2) > 0)) {
        my ($k1, $v1, $l1) = split '\t', $db1[0];
        my ($k2, $v2, $l2) = split '\t', $db2[0];
        my  $res;

        if      (!defined($k1)) {
            $res = computeValue($testName, $k2, 2, $v2, $l2);
            shift @db2;
        }
        elsif (!defined($k2)) {
            $res = computeValue($testName, $k1, 1, $v1, $l1);
            shift @db1;
        }
        elsif ($k1 eq $k2) {
            $res = computeValue($testName, $k1, 1, $v1, $l1, 2, $v2, $l2);
            shift @db1;
            shift @db2;
        }
        elsif ($k1 lt $k2) {
            $res = computeValue($testName, $k1, 1, $v1, $l1);
            shift @db1;
        }
        else {
            $res = computeValue($testName, $k2, 2, $v2, $l2);
            shift @db2;
        }

        if (defined($res)) {
            #print STDERR "RESULT: $res\n";
            push @dbT, $res;
        }
    }

    #  Compare.

    my @errors = compareKmers($testName);

    if (scalar(@errors) > 0) {
        foreach my $e (@errors) {
            print STDERR "  $e\n";
        }

        foreach my $t (@dbM) {print "dbM '$t'\n";}
        foreach my $t (@dbT) {print "dbT '$t'\n";}

        print STDERR "\n";
    }
    else {
        print STDERR "  Success!\n";
        print STDERR "\n";
    }
}





createDatabase(1, $in1, $db1);
createDatabase(2, $in2, $db2);

testMethod("union");
testMethod("union-min");
testMethod("union-max");
testMethod("union-sum");

testMethod("intersect");
testMethod("intersect-min");
testMethod("intersect-max");
testMethod("intersect-sum");

testMethod("subtract");
testMethod("difference");
