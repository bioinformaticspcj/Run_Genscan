#!/usr/bin/perl -w
#author pcj 2019.4.20 
#2019.8.7 modified
use strict;
use Getopt::Long;
my($file,$thresh0,$h,$dir,$spe);
GetOptions( 
        "fa|f:s"=>\$file, 
        "thresh|t:i"=>\$thresh0, 
        "help|?"=>\&USAGE,
        "dir:s"=>\$dir,
	"species|s:s"=>\$spe,
);
&USAGE unless ($file && $dir);
mkdir "work.seqs";
mkdir "genscan.out";
$spe="HumanIso.smat",unless($spe);
open(FA,"$file") or die "please input <*.fa file> <treshold>( num defualt:7000000bp)";
$/=">";
my $thresh=$thresh0||=7000000;
my($head,$seq,@body,$i,$j);
while (<FA>) {
    chomp;
    next,unless($_);
    @body=split/\n/,$_;
    $body[0]=~s/\s+/_/g;
    $head=$body[0];
    $seq=join"",@body[1..$#body];
    my$len=length($seq);
    if ($len > $thresh) {
        my $st=$len%$thresh;
        for($i=0;$i<=$len;$i=$i+$thresh){
            $j++;
            my $seq1=substr($seq,$i,$thresh);#if($i+$thresh<=$len);
            my $seq_c=&cut($seq1);
            my $t=length($seq1);
            open(OUT,">work.seqs/${head}-$j.fa") or die "can not open ";
            print OUT">${head}-$j\n$seq_c\n";
        }
    
    }
    else{open(OUT1, ">work.seqs/$head.fa") or die "can not open ";print OUT1">$head\n$seq\n;"}
    $j=0;
}
sub cut{
my $seqq=shift;
my $te=length($seqq)%80;
my @out=$seqq=~/.{80}/g;
   my$last=substr($seqq,length($seqq)-$te,$te);
   push @out,$last;
my $ou=join"\n",@out;
    return $ou;
}
 
mkdir "works";
my @files=glob("work.seqs/*.fa");
for(@files){
	system("${dir}genscan ${dir}$spe $_ >> genscan.out/total.fa.genscan 2>> genscan.out/total.fa.genscan.log");
}
`perl Deal_outs.pl ./genscan.out/total.fa.genscan $thresh`;
######################################print usage
sub USAGE{
    my $usage=<<"USAGE";
Name:
$0 ----Genscanpip.pl   
Description:
You can use this script to run genscan and cut the too long sequences to shorter segments and format outputs to gff format.  
Note: Predicts span two or more sequences in your input file can not be formatted !All needed files should be 
     putted in the same directory(*.fa,Deal_outs.pl,Genscanpip.pl )! When rerun the script, you should rm works.seqs,works,
     and genscan.out directorys! All results are putted into genscan.out directory.
Usage:
  options:
  -fa|f     <only one>          infile(*.fasta)
  -thresh|t <int number>(bp)        ( example : 8000000|default:7000000) 
  -d        <directory where genscan and its .samt file can be found> 
  -s        <select specie model you need>(default: HumanIso.smat)
  -h        Help                (print usage)

Example:
perl $0 -f *.fa -d *.Genscan/ -l -s HumanIso.smat

USAGE
  print $usage;
  exit;    
}
#system("qsub-seg.pl --resurce vf=1G --lines 2 --maxjobs 200 works/work.sh");
######################################################################################################
######################################################################################################

