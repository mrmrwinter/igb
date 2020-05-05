#! /usr/bin/perl
use warnings;
use strict;
my $headerfile = 'snakemake.input[1]';
my $input = 'snakemake.input[0]';
open( HEADERFILE, '&lt;', $headerfile ) or die $!;
chomp ( my @headers = map { split } &lt;$headerfile&gt; );    #splitting lines on whitespaces.
close HEADERFILE;
my %seqs;
open( INPUTFILE, '&lt;', $input ) or die $!;
{
local $/ = '';         #Reading until blank line
while ( &lt;$input&gt; ) {
     my ( $header, $sequence ) = m/&gt;\s*(\S+)\n(.*)/ms;
     $sequences{$header} = $sequence;
}
open( my $seqsfile, "&gt;", "input.fa" );
foreach my $header (@headers) {
             if ( $sequences{$header} ) {
                       print $header, "\n";
                       print $sequences{$header}, "\n";
             }
}

close( $seqsfile );
}

close INPUTFILE;
exit;
