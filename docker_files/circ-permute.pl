#!/usr/bin/env perl

 #use lib '/home/ec2-user/miniconda3/lib/perl5/site_perl/5.22.0';
 use strict;
 use Bio::SeqIO;
 use Getopt::Long;
 my $splitPoint=600;
 my $inFile=undef;
 my $reverse=0;
 my $outFile="permute_scaffold.fasta";
 my $prefix="scaffold_";
 my $maxPos=undef;
 my $regexp=".";
 my $targetId=undef;

 &GetOptions('split=i'=>\$splitPoint,'in=s'=>\$inFile, 'out=s'=>\$outFile,'rev'=>\$reverse, 'reg=s'=>\$regexp,
       'pre=s'=>\$prefix,'max=i'=>\$maxPos,'id=s'=>\$targetId);
 $inFile=shift if (!defined $inFile && scalar(@ARGV)==1);

 unless ( defined $inFile )
 {
   die "circ-perm.pl -in inFile -out outFile -split spsplitPoint (-r)\n";
 }

 my $rdr=new Bio::SeqIO(-file=>$inFile);
 my $writer=new Bio::SeqIO(-format=>'fasta',-file=>">$outFile");
 while (my $rec=$rdr->next_seq)
 {
   if ($rec->id=~/$regexp/ || (defined $targetId && $rec->id eq $targetId))
   {
     print $rec->id, "  ", $rec->length,  "\n";
     my $headSeq=$rec->subseq(1,$splitPoint);
     $maxPos=$rec->length; #unless (defined $maxPos);
     #print $maxPos, "\n";
     my $tailSeq=$rec->subseq($splitPoint+1,$maxPos);
     my $revFlag=""; $revFlag=".rc" if ($reverse);
     my $newSeq=new Bio::Seq(-id=>$prefix.$rec->id.".cp.$splitPoint$revFlag",-seq=>$tailSeq."N".$headSeq);
     $newSeq=$newSeq->revcom if ($reverse);
     $writer->write_seq($newSeq);
   }
   else
   {
     $writer->write_seq($rec);
   }
 }
