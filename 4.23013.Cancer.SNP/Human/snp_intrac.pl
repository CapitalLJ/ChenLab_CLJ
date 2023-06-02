#!/usr/bin/perl

use strict;
use warnings;

my $filename = shift @ARGV;
open my $fh, '<', $filename or die "Could not open file: $!\n";

while (my $line = <$fh>) {
    chomp $line;
    my @columns = split /\t/, $line;
    my $ones_count = 0;
    for (my $i = 9; $i < scalar @columns; $i++) {
        if ($columns[$i] eq '0|0') {
            $columns[$i] = '0';
        }
        elsif ($columns[$i] eq '1|0' or $columns[$i] eq '0|1') {
            $columns[$i] = '1';
        }
        elsif ($columns[$i] eq '1|1') {
            $columns[$i] = '2';
        }
        $ones_count += () = $columns[$i] =~ /1/g;
    }
    unshift @columns, $ones_count;
    my $new_line = join "\t", @columns;
    print "$new_line\n";
}
close $fh;