# File information
# cd ~/UsV/R172_Mutation
# perl ReadGenename.pl

#use strict;
#use warnings;
#use List::Util qw( min max sum);# statement for using the min & max function
#use File::Basename;
#use Data::Dumper;
#use HTML::TableExtract;
#use File::Slurp;
#use Statistics::Lite qw(mean);

#############################################
#Files
$reffile="NC_000913.gff3"; #from genbank AP012306.gff3
$outfile="MG1655BiotypeList.txt";
#$outfile_temp="MDS42CDSList_temp.txt";
#############################################

my %Synname=();
open (REF,"<$reffile") or die "$!";
open (OUT,">$outfile") or die "$!";
print OUT "Start\tEnd\tBiotype\tStrand\tName\tLocus_tag\tProduct\tECKnb\n";

$line=<REF>;
$line=<REF>;
$line=<REF>;
my $nb=0;
while($line=<REF>)
{
	chomp $line;
	@b=split(/\t/,$line);
	if ($b[0]=~/\#\#/) # ##FASTA
	{
		last;
	}
	if ($b[2]eq"gene")
	{
		my $locus_tag=$b[8];
		$locus_tag=~/locus_tag=(.*)/;
		$locus_tag=$1;
		my $name=$b[8];
		$name=~/gene=(.*?);/;
		$name=$1;
		my $gene_biotype=$b[8];
		$gene_biotype=~/gene_biotype=(.*?);/;
		$gene_biotype=$1;
		my $strand="C";
		if ($b[6]eq"+")
		{
			$strand="D";
		}
		my $tag=$b[3].$b[4].$gene_biotype.$strand.$name.$locus_tag;
		#print "$tag\n";
		my $gene_synonym=$b[8];
		$gene_synonym=~/gene_synonym=(.*?);/;
		$Synname{$tag}=$1;
		#print "$Synname{$tag}\n";
		@Syns=split(/\,/,$Synname{$tag});
		print "@Syns\n";
		my @ECKNumb=grep(/^ECK/,@Syns);
		print OUT "$b[3]\t$b[4]\t$gene_biotype\t$strand\t$name\t$locus_tag\t$product\t$ECKNumb[0]\n";
		$nb++;
	}
}
close REF;
close OUT;
