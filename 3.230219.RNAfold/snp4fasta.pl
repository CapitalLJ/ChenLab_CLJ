#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h' => sub { Getopt::Long::HelpMessage(0) },
    'in|i=s' => \my $in_list,
    'fa|f=s' => \my $in_fa,
    'stdin'  => \my $stdin,
) or Getopt::Long::HelpMessage(1);

my %fasta;
my $title_name;
open my $FA, "<", $in_fa;
while (<$FA>) {
    if (/^>(\S+)/) {
        $title_name = $1;
    }
    else {
        $_ =~ s/\r?\n//;
        $fasta{$title_name} .= $_;
    }
}
close($FA);

if ( defined($in_list) ) {
    open( my $SEG, "<", $in_list );
    while (<$SEG>) {
        s/\r?\n//;
        my ( $title, $pos, $origin, $snp ) = split( /\s+/, $_ );
        if ( exists( $fasta{$title} ) ) {
            my @seq = split( "", $fasta{$title} );
            if ( $origin eq $seq[ $pos - 1 ] ) {
                splice( @seq, $pos - 1, 1, $snp );
                print(">$title\_$pos\_$origin\_$snp\n");
                print( join( "", @seq ) );
                print("\n");
            }
            else {
                warn("Sorry, $title-$pos-$origin doesn't match $seq[$pos - 1]\n"
                );
            }
        }
        else {
            warn("Sorry, there is no such a segment: $_\n");
        }
    }
    close($SEG);
}
elsif ( defined($stdin) ) {
    while (<>) {
        s/\r?\n//;
        my ( $title, $pos, $origin, $snp ) = split( /\s+/, $_ );
        if ( exists( $fasta{$title} ) ) {
            my @seq = split( "", $fasta{$title} );
            if ( $origin eq $seq[ $pos - 1 ] ) {
                splice( @seq, $pos - 1, 1, $snp );
                print(">$title\_$pos\_$origin\_$snp\n");
                print( join( "", @seq ) );
                print("\n");
            }
            else {
                warn("Sorry, $title-$pos-$origin doesn't match $seq[$pos - 1]\n"
                );
            }
        }
        else {
            warn("Sorry, there is no such a segment: $_\n");
        }
    }
}
else {
    die("==> You should provide the LIST!!\n");
}

__END__

