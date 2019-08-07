#!/usr/bin/perl -w
#author pcj 2019.4.21-5.7
use strict;
my ($seq,%lines,$head,$j,$k,$st,$p,$n);my($s,$ttt,$ttt1,$r,%lines1,$cds,$pcj,$id1,$id2,@f,@s,$for,@end,@sort1,@sort2);
open(IN,"$ARGV[0]")or die "please input genscan.out file! and split_size(default:7000000) example: perl $0 *.genscan.out 8000000\n";
open(OUT,">$ARGV[0].gff")or die "please input genscan.out file! and size(default:7000000) example: perl $0 *.genscan.out 8000000\n";
 my $flag=0;
 my @all=<IN>;
 my $o=0;
for($s=0;$s<=$#all;$s++){                                                             #read genscan.out file
    chomp($all[$s]);
    if($all[$s]=~/^Sequence\s+(\S+)/){                                                #find seq name
        my @na=split(/-/,$1);
        $ttt=scalar(@na);
        if ($ttt==2) {
            $r=0;
            for( $o=$s;$o<=$#all;$o++){
                chomp($all[$o]);
                if($all[$o]=~/^Sequence\s+(\S+)/){$r++;$seq=$1;}
                
                my @na1=split(/-/,$seq);
                 $ttt1=scalar(@na1);                                                 #judge divided or not
                if($ttt1!=2){
                    $r--;
                    last;
                    }
                last,if($r>$na1[1]);
                 next,unless($all[$o]=~/^\d+/||$all[$o]=~/^\s+\d+/);
                 
            my @line=split(" ",$all[$o]);
            $head=$na1[0]."\t"."Genscan" ,if($na1[0]);
            push @{$lines{$r}},$all[$o];    
            }
            $s=$o-1;
           my @out2=&sort1(%lines);                                                 #output
           undef %lines;
           &put_out($head,@out2);          
        }
        if($ttt==1){
            undef %lines1;                                                          #judge divided or not 
            for( $o=$s;$o<=$#all;$o++){
                chomp($all[$o]);
                if($all[$o]=~/^Sequence\s+(\S+)/){$pcj++;$seq=$1;}
                
                my @na1=split(/-/,$seq);
                 $ttt1=scalar(@na1);
                if($ttt1!=1){
                    last;
                    }
                if($pcj>1){
                    $pcj=0;
                    last;
                    }
                 next,unless($all[$o]=~/^\d+/||$all[$o]=~/^\s+\d+/);
             my @line=split(" ",$all[$o]);
             
            $head=$seq."\t"."Genscan" ,if($seq);
            @sort2=split/\./,$line[0];
            $lines1{$sort2[0]}{$sort2[1]}=$all[$o];         
            }
            $s=$o-1;
            my @result1=&sort2(%lines1),if(%lines1);
            &put_out($head,@result1),if(%lines1);
        }
         
    }  
}
################################################### sort the lines by the first col in divided seq
sub sort1{
    my %li=@_;
    my @out=();
    my $treshold=$ARGV[1]||7000000;
    foreach  my $rr (sort{$a<=>$b}keys %li){
        my %arry1=();
          foreach my $g(@{$li{$rr}}){
            my @arry=split" ",$g;
            $arry[3]=$arry[3]+($rr-1)*$treshold;
            $arry[4]=$arry[4]+($rr-1)*$treshold;
            push(@arry,$rr);
            my $c=join"\t",@arry;
            @sort1=split/\./,$arry[0];
           $arry1{$sort1[0].$rr}{$sort1[1]}=$c;
          }
          foreach  my $rr1 (sort{$a<=>$b} keys %arry1){
            foreach my $rr2 (sort{$a<=>$b} keys %{$arry1{$rr1}}){
                
                 push @out,${arry1{$rr1}}{$rr2};
            }
        }
    }
     return @out;
}
################################################### sort the lines by the first col in undivided seq
sub sort2{
    my %li=@_;
    my @result=();
    foreach my $l (sort{$a<=>$b} keys %li){
        foreach my $ll(sort{$a<=>$b}keys %{$li{$l}}){
            push @result, "${$li{$l}}{$ll}\n";
        }
    }
       return @result;
}
################################################### formatt the out put
sub put_out{
    my ($head1,@out1)=@_;#read data from main function
    my($na,$sorft)=split(" ",$head1);
        for($j=0;$j<=$#out1;$j++){
            if($out1[$j]=~/Prom/||$out1[$j]=~/Init/){ # get the lines between "Prom" and "PolyA"
                $j++,if($out1[$j]=~/Prom/);   
                $n=0;
                $p=0;                                # flags for positive seq and negative seq
                $cds=0;
                    if($out1[$j]){
                            my @out2=();
                            $k++; 
                            for($for=$j;$for<=$#out1;$for++){
                            last,if($out1[$for]=~/PlyA/);
                            @f=split" ",$out1[$for];
                            $cds++;
                           $id1=$1,if($f[0]=~/(\d+)\./);
                            my @re=split(" ",$out1[$for]);
                            my @re1=split(" ",$out1[$for-1]),if($out1[$for-1]);
                                if ($re[1]eq "Init"||$re[1]eq"Sngl") {
                                    if ($re[2]eq"+") {
                                        $p=1;
                                        $st=$re[3];
                                        push(@out2,"$head\tCDS\t$re[3]\t$re[4]\t.\t$re[2]\t0\tID=CDS$cds.$na.t$k;Parent=$na.t$k");  
                                    }
                                    if ($re[2]eq"-") {
                                        $n=1;
                                        $st=$re[3];
                                        push(@out2,"$head\tCDS\t$re[4]\t$re[3]\t.\t$re[2]\t0\tID=CDS$cds.$na.t$k;Parent=$na.t$k") ;
                                    }
                                }
                                else{
                                    push(@out2,"$head\tCDS\t$re[3]\t$re[4]\t.\t$re[2]\t$re1[7]\tID=CDS$cds.$na.t$k;Parent=$na.t$k"),if($p==1); 
                                    push(@out2,"$head\tCDS\t$re[4]\t$re[3]\t.\t$re[2]\t$re1[7]\tID=CDS$cds.$na.t$k;Parent=$na.t$k"),if($n==1);
                                }
                            if(!$out1[$for]){$k--;
                                    last;
                            }
                            @s=split" ",$out1[$for+1],if($out1[$for+1]);
                            $id2=$1,if($s[0]=~/(\d+)\./);
                            if ($id1!=$id2) {
                                $k--;
                               undef @out2;
                               last;
                            }                          
                        }   
                        if($for>$#out1){
                           $j=--$for;
                        }
                        else{ $j=$for;}
                        if($out1[$j]=~/PlyA/){
                             @end=split(" ",$out1[$j-1]);
                        }
                        else{@end=split(" ",$out1[$j]);}
                        unshift(@out2,"$head\tmRNA\t$st\t$end[4]\t.\t$end[2]\t.\tID=$na.t$k;Parent=$na.gene$k"),if($end[2]eq"+"&&@out2);# add title for each trans
                        unshift(@out2,"$head\tgene\t$st\t$end[4]\t.\t$end[2]\t.\tID=$na.gene$k"),if($end[2]eq"+"&&@out2);# add title for each trans
                        unshift(@out2,"$head\tmRNA\t$end[4]\t$st\t.\t$end[2]\t.\tID=$na.t$k;Parent=$na.gene$k"),if($end[2]eq"-"&&@out2); # add title for each trans
                        unshift(@out2,"$head\tgene\t$end[4]\t$st\t.\t$end[2]\t.\tID=$na.gene$k"),if($end[2]eq"-"&&@out2); # add title for each trans
                        print OUT"$_\n",for(@out2);
                        undef @out2;
                    }
                }   
        }    
}

    

