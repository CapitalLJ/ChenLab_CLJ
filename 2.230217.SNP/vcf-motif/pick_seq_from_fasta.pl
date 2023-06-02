#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub SEQ_REV_COMP {
    my $SEQ = reverse shift;
    $SEQ =~ tr/Uu/Tt/;
    return ( $SEQ =~ tr/AGCTagct/TCGAtcga/r );
}

sub SEQ_TR_TU {
    my $SEQ = shift;
    return ( $SEQ =~ tr/Uu/Tt/r );
}

open( my $FASTA, "<", $ARGV[0] );
open( my $SEG,   "<", $ARGV[1] );

my %fasta;
my $title_name;
while (<$FASTA>) {
    if (/^>(\S+)/) {
        $title_name = $1;
    }
    else {
        $_ =~ s/\r?\n//;
        $fasta{$title_name} .= $_;
    }
}
close($FASTA);

while (<$SEG>) {
    s/\r?\n//;
    my ( $chr, $start, $end, $dir, $name ) = split( /\s+/, $_ );

    if ( exists( $fasta{$chr} ) ) {
        my $length = abs( $end - $start ) + 1;
        my $seq    = substr( $fasta{$chr}, $start - 1, $length );

        if ( $dir eq "-" ) {
            $seq = SEQ_REV_COMP($seq);
        }
        else {
            $seq = SEQ_TR_TU($seq);
        }

        if ( defined($name) ) {
            print ">$chr:$start-$end($dir)$name\n$seq\n";
        }
        else {
            print ">$chr:$start-$end($dir)\n$seq\n";
        }

    }
    else {
        warn("Sorry, there is no such a segment: $_\n");
    }
}
close($SEG);

__END__
