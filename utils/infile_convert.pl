#!/usr/bin/perl

#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @update Dec 07 2020
# @license artistic license 2.0
#

use warnings;
use strict;
use Getopt::Std;


sub HELP_MESSAGE { die "
Usage:
   $0 [options] genes.txt blast.txt ... > blast_mytaxa.txt

   genes.gff2		File containing the genes in any supported format
   			(see option -f).
   blast.txt ...	One or more tabular BLAST files.
   blast_metaxa.txt	Input file for MyTaxa.

   Options:
   -l <float>		Minimum fraction of the gene aligned to consider a
   			hit.  By default: 0.75.  Ignored if -f 'no'.
   -f <str>		Format of the genes prediction.  Any of:
   			 o gff2: GFF v2 as produced by MetaGeneMark.hmm.
			 o gff3: GFF v3 with id field in the last column.
			 o tab: Tabular file with columns gene, gene length,
			   and contig.
			 o no: Ignores genes file.
			By default: gff2.
   -q			Run quietly.
   -h			Display this message and exit.

";}

my %o;
getopts('l:f:qh',\%o);
my($gff, @blasts) = @ARGV;
($gff and $#blasts>=0) or &HELP_MESSAGE;
$o{h} and &HELP_MESSAGE;
$o{f} ||= "gff2";
$o{f} = lc $o{f};
$o{l} ||= 0.75;

my %gene;
if($o{f} ne 'no'){
   print STDERR "Reading genes collection.\n" unless $o{q};
   open GFF, "<", $gff or die "Cannot read file: $gff: $!\n";
   while(<GFF>){
      next if /^#/;
      next if /^\s*$/;
      chomp;
      my @ln = split /\t/;
      if($o{f} eq 'gff2'){
	 exists $ln[8] or die "Cannot parse line $.: $_\n";
	 my $id = $ln[8];
	 $id =~ s/(.+[ ,])?gene_id[ =]/gene_id_/;
	 $id =~ s/[, ].*//;
	 $ln[0] =~ s/ .*//;
	 $gene{$id} = [$ln[0], (1+$ln[4]-$ln[3])/3];
      }elsif($o{f} eq 'gff3'){
	 exists $ln[8] or die "Cannot parse line $.: $_\n";
	 $ln[8] =~ m/id=([^;]+)/i or die "Cannot parse line $.: $_\n";
	 my $id = $1;
	 $ln[0] =~ s/ .*//;
	 $gene{$id} = [$ln[0], (1+$ln[4]-$ln[3])/3];
      }elsif($o{f} eq 'tab'){
         exists $ln[2] or die "Cannot parse line $.: $_\n";
	 $ln[1]+0 or die "$ln[0]: Length zero.\n";
         $gene{$ln[0]} = [$ln[2], $ln[1]/3];
      }else{
         die "Unsupported format: ".$o{f}.".\n";
      }
   }
   close GFF;
}

my $i=0;
my $p=0;
print STDERR "Generating MyTaxa input.\n" unless $o{q};
for my $blast (@blasts){
   print STDERR "  o $blast\n" unless $o{q};
   open BLAST, "<", $blast or die "Cannot read file: $blast: $!\n";
   while(<BLAST>){
      chomp;
      my @l = split /\t/;
      $i++;
      my $ctg;
      if($o{f} eq 'no'){
         $ctg = $l[0];
      }else{
	 my $gene_id = $l[0];
	 $gene_id =~ s/\|.*// if not exists $gene{$gene_id};
	 if(not exists $gene{$gene_id} and $gene_id =~ m/^gene_(\d+)$/){
	   $gene_id = "gene_id_$1";
	 }
	 exists $gene{$gene_id} or die "Cannot find contig for gene $l[0] ($gene_id).\n";
	 next unless $l[3] >= $o{l}*$gene{$gene_id}->[1];
	 $ctg = $gene{$gene_id}->[0];
      }
      $l[1] =~ m/gi\|(\d+)\|/ or die "Cannot parse GI in $l[1].\n";
      print "".join("\t", @l, $ctg, $l[0], $1)."\n";
      $p++;
   }
   close BLAST;
}
print STDERR " Found $i results, reported $p.\n" unless $o{q};

