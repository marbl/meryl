#!/usr/bin/env perl

use strict;

my $nmeryl = "/work/meryl/FreeBSD-amd64/bin/meryl";
my $omeryl = "/work/canu/FreeBSD-amd64/bin/meryl";


sub makeSequence ($$$) {
    my $name   = shift @_;
    my $seq    = shift @_;
    my $seqLen = length($seq);
    my $lines  = shift @_;
    my $bases  = 0;

    if (-e $name) {
        for (my $ii=1; $ii <= $lines; $ii++) {
            $bases += $seqLen * $ii;
        }
    }

    else {
        open(F, "> $name");
        print F ">sequence\n";
        for (my $ii=1; $ii <= $lines; $ii++) {
            print F $seq x $ii;
            print F "\n";

            $bases += $seqLen * $ii;
        }
        close(F);
    }

    print STDERR "Generated '$name' with $bases bp of '$seq'.\n";

    return($bases);
}


{
    my $bases = makeSequence("A.fasta", "A", 5000);
    system("meryl -n $bases -k 22 print count A.fasta output A > A.count.dump");
    system("meryl print A > A.dump");
}

{
    my $bases = makeSequence("AC.fasta", "AC", 5000);
    system("meryl -n $bases -k 22 print count AC.fasta output AC > AC.count.dump");
    system("meryl print AC > AC.dump");
}

{
    my $bases = makeSequence("ACG.fasta", "ACG", 5000);
    system("meryl -n $bases -k 22 print count ACG.fasta output ACG > ACG.count.dump");
    system("meryl print ACG > ACG.dump");
}

{
    my $bases = makeSequence("ACGT.fasta", "ACGT", 5000);
    system("meryl -n $bases -k 22 print count ACGT.fasta output ACGT > ACGT.count.dump");
    system("meryl print ACGT > ACGT.dump");
}
