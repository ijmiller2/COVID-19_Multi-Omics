#!/usr/bin/env perl
use strict;
use FileHandle;
use IO::Uncompress::Gunzip qw(gunzip) ;
use IO::Compress::Gzip qw(gzip) ;
use IO::Handle;
use List::Util qw[min max]; 

print STDERR join(">\t<", time, "\n");
  my @nRE=(0,0,0,0,0,0);
  my %g = ();
  my $g_jStr=("J" x 500);
  my %g_octamerRE=();
  my %g_octamerSingleMismatches=();
  my @g_kSELen;

  my $nFailOverlap = 0;
  my $nFailQ = 0 ;
  my $nFailNs = 0;
  my $nFailAdapter = 0;
  my $nFailQ1 = 0;
  my $nFailNs1 = 0;
  my $nFailAdapter1 = 0;
  my $nFailQ2 = 0;
  my $nFailNs2 = 0;
  my $nFailAdapter2 = 0;
  my $nReadsMT = 0;
  my $nReadsMTSE = 0;
  my $nReadsSE = 0;
  my $nAdapter = 0;
  my $nAdapter1 = 0;
  my $nAdapter2 = 0;
  my $fhOutR1 = openOut("r1.pe.fq.gz");
  my $fhOutR2 = openOut("r2.pe.fq.gz");
  my $fhOutUnmergedR1 = openOut("r1.pe.um.fq.gz");
  my $fhOutUnmergedR2 = openOut("r2.pe.um.fq.gz");
  my $fhOutR1R2 = openOut("r1.r2.seq.gz");
  my $fhOutR1SE = openOut("r1.se.fq.gz");
  my $fhOutR1SESeq = openOut("r1.se.seq.gz");
  my $fhOutMT = openOut("r1.r2.mt.seq.gz");
  my $fhOutMTSE = openOut("r1.r2.mt.se.seq.gz");
  my $fhOutOverlapR1 = openOut("r1.adapter.fq.gz");
  my $fhOutOverlapR2 = openOut("r2.adapter.fq.gz");
  my $fhOutFailQR1 = openOut("failedQuality.r1.fq.gz");
  my $fhOutFailQR2 = openOut("failedQuality.r2.fq.gz");
  my $fhOutFailNsR1 = openOut("failedN.r1.fq.gz");
  my $fhOutFailNsR2 = openOut("failedN.r2.fq.gz");
  my $fhOutFailAdapterR1 = openOut("failedAdapter.r1.fq.gz");
  my $fhOutFailAdapterR2 = openOut("failedAdapter.r2.fq.gz");
  my $fhOutFailOverlapR1 = openOut("failedOverlap.r1.fq.gz");
  my $fhOutFailOverlapR2 = openOut("failedOverlap.r2.fq.gz");


print STDERR join(">\t<", time, "\n");
  ($g{re1}, $g{re2}, $g{re1Long}, $g{re2Long}, $g{re1_pre}, $g{re2_pre}, $g{re1_pre_pos_h}, $g{re2_pre_pos_h}, $g{ad1_short},$g{ad2_short}) = makeREs();
  my %g_baseRC = ("A" => "T", "C" => "G", "G" => "C", "T" => "A", "N" => "N");
  my $g_re1 = $g{re1};
  my $g_re2 = $g{re2};
  my $g_re1Long = $g{re1Long};
  my $g_re2Long = $g{re2Long};
  my $g_re1_pre = $g{re1_pre};
  my $g_re2_pre = $g{re2_pre};
  my $g_re1_pre_pos_h = $g{re1_pre_pos_h};
  my $g_re2_pre_pos_h = $g{re2_pre_pos_h};
  my $g_ad1_short = $g{ad1_short};
  my $g_ad2_short = $g{ad2_short};
  print STDERR join(">\t<", "ad1_short", $g_ad1_short, "\n");
  print STDERR join(">\t<", "ad2_short", $g_ad2_short, "\n");

  my $g_qMin = $g{qMin} = 29;
  my $g_qMinChar = $g{qMinChar} = chr(29 + 33);
  my $g_qMinLen=$g{qMinLen}=31;
  my $g_qRE = $g{qRE}= "(?:[".$g_qMinChar."-P]{".$g_qMinLen.",})";
  my $g_qRE_LOW = $g{qRE_LOW}= "(?:[".chr(12+33)."-P]{".$g_qMinLen.",})";

  my $g_writeUnmerged = 1;

print STDERR join(">\t<", "qRE", $g_qRE, "\n");
print STDERR join(">\t<", "qRE_LOW", $g_qRE_LOW, "\n");
print STDERR $g{re1}."\n";
print STDERR $g{re2}."\n";

sub writeInsertHisto {
  my $fhOutInserts = openOut("insertHistogram.tab");
  my $i = 0;
  for my $k (@g_kSELen) {
    $fhOutInserts->print( "$i\t".($k||0)."\n");
    $i++;
  }
  $fhOutInserts->close();
}

sub writeLog {
  my $nReads = shift @_;
  print STDERR time."\n";
  print STDERR join("\t", "N reads", $nReads)."\n";
  print STDERR join("\t", "Failed Adapter", $nFailAdapter, $nFailAdapter1, $nFailAdapter2)."\n";
  print STDERR join("\t", "Failed Ns", $nFailNs, $nFailNs1, $nFailNs2)."\n";
  print STDERR join("\t", "Failed Q", $nFailQ, $nFailQ1, $nFailQ2)."\n";
  print STDERR join("\t", "Failed Overlap", $nFailOverlap)."\n";
  print STDERR join("\t", "Failed", ($nFailAdapter + $nFailNs + $nFailQ + $nFailOverlap))."\n";
  print STDERR join("\t", "Mito", $nReadsMT)."\n";
  print STDERR join("\t", "Mito SE", $nReadsMTSE)."\n";
  print STDERR join("\t", "N reads PE", $nReads-$nReadsSE)."\n";
  print STDERR join("\t", "N reads SE", $nReadsSE)."\n";
  print STDERR join("\t", "Passed Adapter", $nAdapter-$nFailAdapter1). "\n";
  print STDERR join("\t", "Passed Non-Mito", $nReads - ($nFailAdapter + $nFailNs + $nFailQ + $nFailOverlap + $nReadsMT))."\n";

  print STDERR join("\t", "nRE", @nRE)."\n";
  writeInsertHisto();
}

mainSub();
$fhOutFailQR1->close();
$fhOutFailQR2->close();

$fhOutR1->close();
$fhOutR1SE->close();
$fhOutR1SESeq->close();
$fhOutOverlapR1->close();
$fhOutOverlapR2->close();
$fhOutFailNsR1->close();
$fhOutFailAdapterR1->close();
$fhOutFailOverlapR1->close();

$fhOutR2->close();
$fhOutFailNsR2->close();
$fhOutFailAdapterR2->close();
$fhOutFailOverlapR2->close();

$fhOutR1R2->close();
$fhOutMT->close();
$fhOutMTSE->close();
$fhOutUnmergedR1->close();
$fhOutUnmergedR2->close();
exit;

sub openReadFiles {
  my $fnR1 = shift @_;
  my $fnR2 = shift @_;

  my $fhR1 = openIn($fnR1);
  my $fhR2 = openIn($fnR2);

  return ($fhR1, $fhR2);
}

sub readFQRecord {
  my $fhIn = shift @_;

  my @result = (undef,undef,undef,undef);
  my $endOfFile = undef;

  for (my $i = 0; $i < 4; $i++) {
    $result[$i]=$fhIn->getline();
  }

  if (! $result[3]){$endOfFile = 1}

  return ($endOfFile, @result);
}

sub loadMTSeq {
  my $fnMT = shift @_;

  my (@c) = qw (A C G T);
  my %mtSeq = ();
  my $fhMT = openIn($fnMT);
  my $mtLine  = $fhMT->getline();

  while ($mtLine) {
    chomp $mtLine;
    $mtSeq{$mtLine}++;
    $mtLine  = $fhMT->getline();
  }

  $fhMT->close();

  return \%mtSeq;
}

sub getR1FileNames {
  my @fnR1 = ();
  foreach  (@_) {
    if (m/_R1_001.fastq.gz$/) {
      push @fnR1, $_;
    }
  }
  return (@fnR1);
}

sub getR2FileNames {
  my @fnR2 = ();
  foreach  (@_) {
    if (m/_R2_001.fastq.gz$/) {
      push @fnR2, $_;
    }
  }
  return (@fnR2);
}

sub loadOctamerSingleMismatches {
  my $fhIn=openIn("../filter_refdata/octamerSingleMismatch.tab");
  my $nextLine = $fhIn->getline();
  while ($nextLine) {
    chop $nextLine;
    my ($octamer, @f)=split /\t/, $nextLine;
    $g_octamerSingleMismatches{$octamer}->{$octamer}++;
    for my $f (@f){$g_octamerSingleMismatches{$octamer}->{$f}++}
    $nextLine = $fhIn->getline();
  }
  $fhIn->close();
}

sub mainSub {
  my (@args)=(@ARGV);
  my (@fnR1) = getR1FileNames(@args);
  my (@fnR2) = getR2FileNames(@args);
  my $nReadsFiles = (scalar @fnR1)+(scalar @fnR2);
  my $fnMT = $args[$nReadsFiles];
  my $mtSeqLength = $args[$nReadsFiles + 1] || 46;

  my $qMin = $args[$nReadsFiles + 2] || 0;
  my $qMinLen = $args[$nReadsFiles + 3] || 0;
  my $nToSkip = $args[$nReadsFiles + 4] || 0;
  my $nToRead = $args[$nReadsFiles + 5] || 10000000000;

  my ($nReads1, $nReads2) = (0,0);
  my $nReadsMT = 0;
  my ($mtSeq_h) = loadMTSeq($fnMT);
  loadOctamerSingleMismatches();
  my $nToLog = 1000000;
  my $log10NtoRead = log($nToRead)/log(10);
  if ($log10NtoRead <= 1.0){$nToLog = 1}
  elsif ($log10NtoRead <= 6.0) {
    $nToLog = 10**(int($log10NtoRead-0.21));
  }

  if ($qMin || $qMinLen) {
    if ($qMin) {
      $g_qMin = $g{qMin} = $qMin;
      $g_qMinChar = $g{qMinChar} = chr($qMin + 33);
    }
    if ($qMinLen) {
      $g_qMinLen=$g{qMinLen}=$qMinLen;
    }
    $g_qRE=$g{qRE}= "(?:[".$g_qMinChar."-P]{".$g_qMinLen.",})";
    $g_qRE_LOW = $g{qRE_LOW}= "(?:[".chr(33+12)."-P]{".$g_qMinLen.",})";
print STDERR join(">\t<", "qRE", $g_qRE, "\n");
print STDERR join(">\t<", "qRE_LOW", $g_qRE_LOW, "\n");
  }

  while (my $fnR1 = shift @fnR1) {
    my $nSkipped = 0;
    my $skipRecord = ($nSkipped < $nToSkip);
    my $fnR2 = shift @fnR2;
    my $nReads1ThisFile = 0;
    my ($fhR1, $fhR2)=openReadFiles($fnR1, $fnR2);
    if ($fhR1 && $fhR2) {
      while (($nToSkip > $nSkipped++) && processPair($fhR1, $fhR2, $mtSeqLength, $mtSeq_h, $skipRecord)) {;}
      $skipRecord = 0; 
      writeLog($nSkipped);
      while (($nToRead > $nReads1ThisFile++) && processPair($fhR1, $fhR2, $mtSeqLength, $mtSeq_h, $skipRecord)) {
        $nReads1++;
        $nReads2++;
        if (0 == ($nReads1 % $nToLog)){ writeLog($nReads1) };
      }
    }
    writeLog($nReads1);
  }
}

sub makeRE {
  my $ad = shift @_;
  my $fragLen = shift @_ || 15;
  my %subsSeen = ();
  my $re = "(?:";
  for my $i (1) {
    my $iA = substr($ad, $i, $fragLen);
    for my $j (0..($fragLen - 2)) {
      my $jA = $iA;
      substr($jA, $j, 1) = ".";
      for my $k (($j+1)..($fragLen - 1)) {
        my $kA = $jA;
        substr($kA, $k, 1) = ".";
        if (! $subsSeen{$kA}++) {
          $re .= "$kA|";
        }
      }
    }
    for my $ij (1..($fragLen - 2)) {
      my $ijA = $iA;
      chop $ijA;
      substr($ijA, $ij, 1) = substr($ijA, $ij, 1) . substr($ijA, $ij, 1);
      for my $j (0..($fragLen - 1)) {
        my $jA = $ijA;
        substr($jA, $j, 1) = ".";
          if (! $subsSeen{$jA}++) {
            $re .= "$jA|";
          }
      }
    }
    for my $ij (1..($fragLen - 1)) {
      my $ijA = substr($ad, $i, ($fragLen + 1));
      substr($ijA, $ij, 1) = "";
      for my $j (0..($fragLen - 2)) {
        my $jA = $ijA;
        substr($jA, $j, 1) = ".";
          if (! $subsSeen{$jA}++) {
            $re .= "$jA|";
          }
      }
    }
  }
  substr($re,-1,1)=")";

  return $re;
}

sub makeREs {

## Tru-seq Adapter sequences
  my $a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
  my $a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
  my $re1_pre = '(?:ATCGG|.{5}AAGAG|.{10}CACAC|.{15}GTCTG)';
  my $re2_pre = '(?:ATCGG|.{5}AAGAG|.{10}CGTCG|.{15}TGTAG)';
  my %re1_pre_pos = (AGATC => 0, GGAAG => 5, AGCAC => 10, ACGTC => 15);
  my %re2_pre_pos = (AGATC => 0, GGAAG => 5, AGCGT => 10, CGTGT => 15);

  my $re1 = makeRE($a1);
  my $re2 = makeRE($a2);
  my $re1Long = makeRE($a1, 24);
  my $re2Long = makeRE($a2, 24);

  return ($re1, $re2, $re1Long, $re2Long, $re1_pre, $re2_pre, \%re1_pre_pos, \%re2_pre_pos,
          "(?:".substr($a1,0,15)."|".substr($a1,1,15).")", 
          "(?:".substr($a2,0,15)."|".substr($a2,1,15).")");
}

sub trimQs {
  my $qual_r = shift @_;
  my $qLen = length $$qual_r;

  my $qs = 0;
  my $qe = 0;
  my $ql = 0;

  if ($$qual_r =~ m/$g_qRE/omgs) {
    $qs = $-[0];
    $qe = $+[0];
    $ql = $qe - $qs;
    if (($qLen - $qe) > $ql) {
      if ($$qual_r =~ m/$g_qRE/omgs) {
        my $qs2 = $-[0];
        my $qe2 = $+[0];
        my $ql2 = $qe2-$qs2;
        if ($ql2 > $ql) {
          ($qs, $qe, $ql) = ($qs2, $qe2, $ql2);
        }
      }
    }
    --$qe;
  }
  pos($$qual_r) = 0;
  return ($qs, $qe, $ql);
}

sub checkQuality {
  my ($s1_r, $s2_r, $q1_r, $q2_r) = (@_);

  my $seqLen = length $$s1_r;
  my ($junkIt, $failQual1, $failQual2) = (0,0,0);
  my ($iqBest1, $jqBest1, $lqBest1) = trimQs($q1_r, $seqLen);
  my ($iqBest2, $jqBest2, $lqBest2) = trimQs($q2_r, $seqLen);
  if (($lqBest1 < $g_qMinLen) || ($lqBest2 < $g_qMinLen)) {
    $junkIt = 1;
    if ($lqBest1 <= $g_qMinLen) { $failQual1 = 1}
    if ($lqBest2 <= $g_qMinLen) { $failQual2 = 1}
  }
  else {
  }

  return ($junkIt, $failQual1, $failQual2,$iqBest1, $jqBest1, $lqBest1,$iqBest2, $jqBest2, $lqBest2);
}

sub checkAdapterREs {
  my $s1 = shift @_;
  my $s2 = shift @_;

  my $junkIt = 0;
  my ($re1Match, $re1MatchPosition, $re1LongMatch, $re1LongMatchPosition) = (0,-1,0,-1);
  my ($re2Match, $re2MatchPosition, $re2LongMatch, $re2LongMatchPosition) = (0,-1,0,-1);
  my $re1PreMatch = ($s1 =~ m/$g_re1_pre/o);
  my $re1PreMatchPosition = $re1PreMatch ? $-[0] : -1;
  my $re2PreMatch = ($s2 =~ m/$g_re2_pre/o);
  my $re2PreMatchPosition = $re2PreMatch ? $-[0] : -1;
  my $preMatchLength = 6;

  if (! ($re1PreMatch || $re2PreMatch)) {
    # no further matches possible, so nothing more to do !
  }
  else {
    my $ssREMatchStart = 0;
    if ($re1PreMatch && $re2PreMatch) {
      $re1PreMatchPosition = $re1PreMatchPosition - $g_re1_pre_pos_h->{substr($s1, $re1PreMatchPosition, $preMatchLength)};
      $re2PreMatchPosition = $re2PreMatchPosition - $g_re2_pre_pos_h->{substr($s2, $re2PreMatchPosition, $preMatchLength)};
      $ssREMatchStart = ($re1PreMatchPosition < $re2PreMatchPosition)
                      ? $re2PreMatchPosition
                      : $re1PreMatchPosition ;
    }
    elsif ($re1PreMatch) {$ssREMatchStart = $re1PreMatchPosition - $g_re1_pre_pos_h->{substr($s1, $re1PreMatchPosition, $preMatchLength)};}
    else {$ssREMatchStart = $re2PreMatchPosition - $g_re2_pre_pos_h->{substr($s2, $re2PreMatchPosition, $preMatchLength)};}

    $ssREMatchStart = ($ssREMatchStart > 2)
                    ? $ssREMatchStart - 2
                    : 0;

    my $ss1 = ($ssREMatchStart == 0) ? $s1 : substr($s1, $ssREMatchStart);
    my $ss2 = ($ssREMatchStart == 0) ? $s2 : substr($s2, $ssREMatchStart);

    $re1MatchPosition = index($ss1,$g_ad1_short);
    if ($re1MatchPosition > 0) {
      $re1Match = 1;
      $re1MatchPosition--;
    }
    elsif ($re1MatchPosition == 0) {
      $re1Match = 1;
    }
    else {
      $re1Match = ($ss1 =~ m/$g_re1/o);
      $re1MatchPosition = $re1Match ? $-[0] : -1;
      if ($re1MatchPosition > 0) {$re1MatchPosition--}
    }
    $re2MatchPosition = index($ss2,$g_ad2_short);
    if ($re2MatchPosition > 0) {
      $re2Match = 1;
      $re2MatchPosition--;
    }
    elsif ($re2MatchPosition == 0) {
      $re2Match = 1;
    }
    else {
      $re2Match = ($ss2 =~ m/$g_re2/o);
      $re2MatchPosition = $re2Match ? $-[0] : -1;
      if ($re2MatchPosition > 0) {$re2MatchPosition--}
    }

    if ($re1Match && $re2Match ) {
# we have found matching adapter hits, YAY
      if (2 >= abs($re1MatchPosition - $re2MatchPosition)) {
        $re1MatchPosition += $ssREMatchStart; #$re1PreMatchPosition;
        $re2MatchPosition += $ssREMatchStart; #$re2PreMatchPosition;
      }
      else {
# we have found  non-matching adapter hits
# whatever this fragment is, it is a mess, and we don't want to deal with it
        $junkIt = 1;
      }
    }
    elsif ($re1Match) {
      my $sss1 = ($re1MatchPosition == 0) ? $ss1 : substr($ss1, $re1MatchPosition);
      $re1LongMatch = ($sss1 =~ m/^$g_re1Long/o);
      $re1LongMatch = ($ss1 =~ m/^$g_re1Long/o);
      $re1LongMatchPosition = $re1LongMatch ? $-[0] : -1;
      if ($re1LongMatchPosition >= 0) {
        $re1LongMatchPosition += $re1MatchPosition + $re1PreMatchPosition;
      }
      $re1MatchPosition += $ssREMatchStart;
    }
    elsif ($re2Match) {
      my $sss2 = ($re2MatchPosition == 0) ? $ss2 : substr($ss2, $re2MatchPosition);
      $re2LongMatch = ($sss2 =~ m/^$g_re2Long/o);
      $re2LongMatchPosition = $re2LongMatch ? $-[0] : -1;
      if ($re2LongMatchPosition >= 0) {
        $re2LongMatchPosition += $re2MatchPosition + $re2PreMatchPosition;
      }
      $re2MatchPosition += $ssREMatchStart;
    }
  }

  return ($junkIt, $re1Match, $re1MatchPosition, $re2Match, $re2MatchPosition,
          $re1LongMatch, $re1LongMatchPosition, $re2LongMatch, $re2LongMatchPosition);
}

# DETECT READS WHERE THE FRAGMENT LENGTH < READLENGTH!
# ALLOWS US TO TRIM ADAPTER SEQUENCE EVEN WHEN THE ADAPTER BASECALLS
# ARE CORRUPTED, OR TOO SHORT TO DISTINGUISH!
# Look for reads where the fragment length > 32 but less than the read length
# In this case, both reads will have adapter sequence at their tails.
# Thus, R2RC will begin with some chunk of the 3' end of the 5' adapter (A1),
#   followed by a sequence that will should map 1-1 with the 5' end of R1.
# Detecting this situation is complicated by the presence of sequencing
#   artifacts (basecall errors, plus indels at the start of reads).
#
# We hunt for overlapping sequence between R1 and R2RC
# If we find overlap, we need to discard the trailing bases
#   that are actually adapter sequence.
# Things are often hinky near both ends of any read, so we're going to start
#   3 bases into R1, and search for overlaps in R2 that start in bases 1+.
sub checkInsertShorterThan2xRead {
  my $s1_r = shift @_;
  my $rcSeq2_r = shift @_;
  my $minOverlap = shift @_;

  my $seqLen = length $$s1_r;
  my $insertLength = 0;
  my $totalMatches = 0;
  my $nMismatchesInStreak = 0;
  my $iLastMatch = 0;

#somewhat experimental
  my $tempR = $s1_r;
  $s1_r=$rcSeq2_r;
  $rcSeq2_r=$tempR;

  my $ss1A=substr($$s1_r,1,8);
  my $ss2A=substr($$rcSeq2_r,-9,8);
  my $ixR1inR2=index(substr($$rcSeq2_r,1,$seqLen-$minOverlap+8), $ss1A);
  my $ixR2inR1 = index(substr($$s1_r,$minOverlap-9),$ss2A);

  if (($ixR1inR2 >= 0) || ($ixR2inR1 >= 0)){
    if ($ixR2inR1 >= 0){
      $ixR2inR1 += $minOverlap - 9 + 1;
    }

    if ($ixR1inR2 == ($seqLen - $ixR2inR1 - 8)) {
      # it's good (check midInsert!)
      $insertLength = $seqLen - $ixR1inR2;
    }
    elsif (($ixR1inR2>=0) && ($g_octamerSingleMismatches{$ss2A}->{substr($$s1_r, ($seqLen - $ixR1inR2 - 8 - 1), 8)})) {
      # it's good!
      $insertLength = $seqLen - $ixR1inR2;
    } 
    elsif (($ixR2inR1>=0) && ($g_octamerSingleMismatches{$ss1A}->{substr($$rcSeq2_r,($seqLen - $ixR2inR1 - 8 + 1), 8)})) {
      # it's  good!
      $insertLength = $ixR2inR1 + 8;
    }
  }
  if ($insertLength > 16) {
    my $midInsert = int ($insertLength / 2);
    if ($g_octamerSingleMismatches{substr($$s1_r,$midInsert-4,8)}->{substr($$rcSeq2_r, $seqLen-$insertLength+$midInsert-4,8)}) {
       # it's good
    }
    else {
      $insertLength = 0;
    }
  }

  return ($insertLength);
}

sub checkInsertShorterThanRead {
  my $s1_r = shift @_;
  my $rcSeq2_r = shift @_;
  my $minOverlap = shift @_;

  my $seqLen = length $$s1_r;
  my $insertLength = 0;
  my $totalMatches = 0;
  my $nMismatchesInStreak = 0;
  my $iLastMatch = 0;

  my $ss1A=substr($$s1_r,1,8);
  my $ss2A=substr($$rcSeq2_r,-9,8);
  my $ixR1inR2=index(substr($$rcSeq2_r,1,$seqLen-$minOverlap+8), $ss1A);
  my $ixR2inR1 = index(substr($$s1_r,$minOverlap-9),$ss2A);

  if (($ixR1inR2 >= 0) || ($ixR2inR1 >= 0)){
    if ($ixR2inR1 >= 0){
      $ixR2inR1 += $minOverlap - 9 + 1;
    }

    if ($ixR1inR2 == ($seqLen - $ixR2inR1 - 8)) {
      # it's good (check midInsert!)
      $insertLength = $seqLen - $ixR1inR2;
    }
    elsif (($ixR1inR2>=0) && ($g_octamerSingleMismatches{$ss2A}->{substr($$s1_r, ($seqLen - $ixR1inR2 - 8 - 1), 8)})) {
      # it's good!
      $insertLength = $seqLen - $ixR1inR2;
    } 
    elsif (($ixR2inR1>=0) && ($g_octamerSingleMismatches{$ss1A}->{substr($$rcSeq2_r,($seqLen - $ixR2inR1 - 8 + 1), 8)})) {
      # it's  good!
      $insertLength = $ixR2inR1 + 8;
    }
  }
  if ($insertLength > 11) {
    my $midInsert = int ($insertLength / 2);
    if ($g_octamerSingleMismatches{substr($$s1_r,$midInsert-4,8)}->{substr($$rcSeq2_r, $seqLen-$insertLength+$midInsert-4,8)}) {
       # it's good
    }
    else {
      $insertLength = 0;
    }
  }

  return ($insertLength);
}

sub mergeReadsWithInsertShorterThan2xReadLength {
  my ($seqLen, $s1_r, $s2_r, $q1_r, $q2_r, $s2RC_r, $overlapOffset, $overlapLength) = (@_);
  
  if ($overlapOffset >= 0) {
    if (substr($$s2RC_r, 0, $overlapLength) eq substr($$s1_r, -$overlapLength)) {
      substr($$q1_r, -$overlapLength) = substr($g_jStr, 0, $overlapLength);
      substr($$q2_r, 0, $overlapLength) = substr($g_jStr, 0, $overlapLength);
    }
    else {
      for (my ($i,$jR,$j) = ($seqLen-$overlapLength, 0, $seqLen-1); $i < $seqLen; $i++, $jR++, $j--) {
        if (substr($$s1_r,$i,1) eq substr($$s2RC_r,$jR,1)) {
          substr($$q1_r,$i,1) = "J";
        }
        elsif (substr($$q1_r,$i,1) lt substr($$q2_r, $j, 1)) {
          substr($$s1_r,$i,1) = substr($$s2RC_r, $jR, 1);
          substr($$q1_r,$i,1) = substr($$q2_r, $j, 1);
        }
        else {
          substr($$s2RC_r,$jR,1) = substr($$s1_r, $i, 1);
          substr($$s2_r,$j,1) = substr($$s1_r, $i, 1);
          substr($$q2_r,$j,1) = substr($$q1_r, $i, 1);
        }
      }
    }
  }
}

sub mergeReadsWithInsertShorterThanReadLength {
  my ($seqLen, $s1_r, $s2_r, $q1_r, $q2_r, $s2RC_r, $re2MatchPosition, $overlapOffset, $overlapLength) = (@_);
  
  if ($overlapOffset == 0) {
    if (substr($$s1_r, 0, $overlapLength) eq substr($$s2RC_r, -$overlapLength)) {
      $$q1_r=substr($g_jStr, 0, $seqLen);
    }
    else {
      for (my ($i,$j,$jR)= (0, $seqLen-$overlapLength, $seqLen-1); $i < $overlapLength; $i++, $j++, $jR--) {
        if (substr($$s1_r,$i,1) eq substr($$s2RC_r,$j,1)) {
          substr($$q1_r,$i,1) = "J";
        }
        elsif (substr($$q1_r,$i,1) lt substr($$q2_r, $jR, 1)) {
          substr($$s1_r,$i,1) = substr($$s2RC_r, $j, 1);
          substr($$q1_r,$i,1) = substr($$q2_r, $jR, 1);
        }
      }
    }
  }
}

# processPair is the primary workhorse of this script
# It picks up an R1/R2 pair from the input streams and attempts
#   to determine whether:
#   A. Both reads are of sufficient quality of the desired length
#   B. Whether the insert was
#     1. shorter than the read length (implying a need to trim adapter sequence)
#     2. or less than twice the read length, minus a minimum overlap (implying the reads overlap)
#     3. or greater than twice the read length
# If the pair appear to overlap, the code will attempt to
#   A. Repair low quality basecalls in the overlapping region of either read
#   B. Bump the quality score for basecalls in the overlapping region that agree
#   C. Merge the pair into a single read that is output to a specific file.
sub processPair {
  my ($fhR1, $fhR2, $mtLength, $mtSeq_h, $skipRecord) = (@_);

  my $junkIt = 0;
  
  my $iBest1 = 0; my $jBest1 = 0; my $lBest1 = 0;
  my $iBest2 = 0; my $jBest2 = 0; my $lBest2 = 0;
  my $iqBest1 = 0; my $jqBest1 = 0; my $lqBest1 = 0;
  my $iqBest2 = 0; my $jqBest2 = 0; my $lqBest2 = 0;
  my $seqLen = 0;

  my $failQual1 = 0; my $failAdapter1 = 0; my $failNs1 = 0;
  my $failQual2 = 0; my $failAdapter2 = 0; my $failNs2 = 0;
  my $failOverlap = 0;
  my $overlap = 0;
  my $overlapOffset = 0;
  my $overlapLength = 0;
 
  my ($endOfFileR1, $h1, $s1, $d1, $q1) = readFQRecord($fhR1);
  my ($endOfFileR2, $h2, $s2, $d2, $q2) = readFQRecord($fhR2);

  if ($endOfFileR1 || $endOfFileR2) {
    print STDERR join("", "failed to getline:\t",
                 $h1||"no h1", $s1||"no s1", $d1||"no d1", $q1||"no q2",
                 $h2||"no h2", $s2||"no s2", $d2||"no d2", $q2||"no q2",
                 "\n");
    return 0;
  }
  elsif ($skipRecord) {
    return 1;
  }
  else {
    chomp $h1; chomp $s1; chomp $d1; chomp $q1;
    chomp $h2; chomp $s2; chomp $d2; chomp $q2;
    $h1 =~ s/\s.*//mgs;
    $h2 =~ s/\s.*//mgs;
    if ($h1 ne $h2) {
      print STDERR join("", "R1 and R2 headers do not match:\t", $h1, $s1, $d1, $q1, $h2, $s2, $d2, $q2, "\n");
      return 0;
    }
    elsif ($mtLength >= 0) {
      my $iLast = (length $s1) - $mtLength;
      my ($mt1, $mt2) = (0,0);
      if ($iLast > 25){$iLast = 25;}
      for (my $i = 0; ($i <= $iLast) && (! ($mt1 && $mt2)); $i+=5) {
        if ((! $mt1) && $mtSeq_h->{substr($s1,$i,$mtLength)}) {$mt1=1}
        if ((! $mt2) && $mtSeq_h->{substr($s2,$i,$mtLength)}) {$mt2=1}
      }
      if ($mt1 && $mt2) {
        $fhOutMT->print(join("\t", $h1, substr($s1, $iBest1, $lBest1), substr($s2, $iBest2, $lBest2))."\n");
        $nReadsMT++;
        return 1;
      }
    }
  }

  $seqLen = length $s1;

  my ($re1Match, $re1MatchPosition, $re2Match, $re2MatchPosition) = (0,-1,0,-1);
  my ($re1LongMatch, $re1LongMatchPosition, $re2LongMatch, $re2LongMatchPosition) = (0,-1,0,-1);

  if (! $junkIt) {
    ($junkIt, $re1Match, $re1MatchPosition, $re2Match, $re2MatchPosition,
          $re1LongMatch, $re1LongMatchPosition, $re2LongMatch, $re2LongMatchPosition) =
      checkAdapterREs($s1, $s2);
  }

  my $s2RC = reverse $s2;
  $s2RC =~ tr/ACGT/TGCA/;

  my $mergedReadsIntoSE = 0;
  if ((! $junkIt) && $re1Match && $re2Match) {
    $mergedReadsIntoSE = 1;
    $overlapOffset = $re2MatchPosition - $re1MatchPosition;
    $overlapLength = min($re1MatchPosition, $re2MatchPosition);
    mergeReadsWithInsertShorterThanReadLength($seqLen, \$s1, \$s2, \$q1, \$q2, \$s2RC, $re2MatchPosition, $overlapOffset, $overlapLength);
  }
  elsif (! $junkIt) {
## Need to check for fragments that contain partial adapter, but not enough for our adapter spotting code
    my ($insertLength) = checkInsertShorterThanRead(\$s1, \$s2RC, ($seqLen - 18));
    if ($insertLength >= $g_qMinLen) {
      $mergedReadsIntoSE = 1;
      $overlapOffset = 0;
      $overlapLength = $insertLength;
      mergeReadsWithInsertShorterThanReadLength($seqLen, \$s1, \$s2, \$q1, \$q2, \$s2RC, $insertLength, $overlapOffset, $overlapLength);
    } 
  }

  if ((! $junkIt) && $mergedReadsIntoSE) {
    if ($overlapLength < $g_qMinLen) {
      $junkIt = 1;
    }
    $nAdapter++;
    $g_kSELen[$overlapLength]++;
  }
  if ((! $junkIt) && (! $mergedReadsIntoSE)) {
    ($overlapLength) = checkInsertShorterThan2xRead(\$s1, \$s2RC, 10);
    if ($overlapLength > 0) {
      $overlapOffset = $seqLen - $overlapLength;
      mergeReadsWithInsertShorterThan2xReadLength ($seqLen, \$s1, \$s2, \$q1, \$q2, \$s2RC, $overlapOffset, $overlapLength);
      ($junkIt, $failQual1, $failQual2,$iqBest1, $jqBest1, $lqBest1,$iqBest2, $jqBest2, $lqBest2) =
        checkQuality(\$s1, \$s2, \$q1, \$q2);
$g_kSELen[$seqLen+$seqLen-$overlapLength]++;
      $mergedReadsIntoSE = 1;
    }
  }

  if ($junkIt && (! ($failQual1 && $failQual2))) { # need to save the fact that it's already known to be junk
    my ($junkIt_orig, $failQual1_orig, $failQual2_orig) = ($junkIt, $failQual1, $failQual2);
    ($junkIt, $failQual1, $failQual2,$iqBest1, $jqBest1, $lqBest1,$iqBest2, $jqBest2, $lqBest2) =
       checkQuality(\$s1, \$s2, \$q1, \$q2);
    $junkIt = 1;
    $failQual1 ||= $failQual1_orig;
    $failQual2 ||= $failQual2_orig;
  }
  elsif ($junkIt && (! ($failQual1 || $failQual2))) {
    my ($junkIt_orig, $failQual1_orig, $failQual2_orig) = ($junkIt, $failQual1, $failQual2);
    ($junkIt, $failQual1, $failQual2,$iqBest1, $jqBest1, $lqBest1,$iqBest2, $jqBest2, $lqBest2) =
       checkQuality(\$s1, \$s2, \$q1, \$q2);
    $junkIt = 1;
  }
  else {
    ($junkIt, $failQual1, $failQual2,$iqBest1, $jqBest1, $lqBest1,$iqBest2, $jqBest2, $lqBest2) =
       checkQuality(\$s1, \$s2, \$q1, \$q2);
  }
  if ($re1Match && ($mergedReadsIntoSE || $junkIt)) {
    $nAdapter1++;
    $failAdapter1 = $junkIt;
  }
  if ($re2Match && ($mergedReadsIntoSE || $junkIt)) {
    $nAdapter2++;
    $failAdapter2 = $junkIt;
  }
  if ((! $mergedReadsIntoSE) && ($failAdapter1 || $failAdapter2)) {
    $nAdapter++;
  }
  

  my $s1Merged = undef;
  my $q1Merged = undef;
  if ((! $junkIt) && $mergedReadsIntoSE) {
    if ($overlapOffset <= 1) {
      $s1Merged = substr($s1,0,$overlapLength);
      $q1Merged = substr($q1,0,$overlapLength);
    }
    else {
      my $q2R = reverse $q2;
      $s1Merged = $s1.substr($s2RC,$overlapLength);
      $q1Merged = $q1.substr($q2R, $overlapLength);
    }
    ($iqBest1, $jqBest1, $lqBest1) = trimQs(\$q1Merged, (length $s1Merged));
    $failQual1 = ($lqBest1 < $g_qMinLen);
    $failQual2 = $failQual1;
    $junkIt = $failQual1;
  }

  if ((! $junkIt) && (!$mergedReadsIntoSE)) {
    $failQual1 = ($lqBest1 < $g_qMinLen);
    $failQual2 = ($lqBest2 < $g_qMinLen);
    $junkIt = ($failQual1 || $failQual2);
  }
  if ((! $junkIt) && (!$mergedReadsIntoSE)) {
    $failQual1 = ($jqBest1 < $g_qMinLen);
    $failQual2 = ($jqBest2 < $g_qMinLen);
  } 

  if ((! $junkIt) && ($mergedReadsIntoSE)) {
    if ( ($mtSeq_h->{substr($s1,2,$mtLength)} && $mtSeq_h->{substr($s2,2,$mtLength)})) {
      $fhOutMTSE->print(join("\t", $h1, $s1Merged, $q1Merged)."\n");
      $nReadsMT++;
      $nReadsMTSE++;
    }
    else {
      $iBest1 = ($iqBest1 < 1) ? 1 : $iqBest1;
      $lBest1 = ($jqBest1 - $iBest1 + 1);
      my $trimmedFrom3Prime = (length $s1Merged)-$lBest1-$iBest1;
      if ($trimmedFrom3Prime < 0){$trimmedFrom3Prime = 0}
      $fhOutR1SESeq->print(join("\t", $h1." t5:$iBest1 t3:$trimmedFrom3Prime", substr($s1Merged, $iBest1, $lBest1), substr($q1Merged, $iBest1, $lBest1))."\n") || die ("Could not write to fhOutR1SESeq :$fhOutR1SESeq:\n");
      $fhOutR1SE->print(join("\n", $h1, substr($s1Merged, $iBest1, $lBest1), "+", substr($q1Merged,  $iBest1, $lBest1),"")) || die ("Could not write to fhOutR1SE :$fhOutR1SE:\n");

      if ($g_writeUnmerged) {
        my $mergedLen = length $s1Merged;
        my $overlapLen = (length $s1Merged)-$seqLen;
        my $r1Len = $seqLen-$iBest1;
        my $r2Len = $seqLen-$trimmedFrom3Prime;
        if (($r1Len >= $g_qMinLen) && ($r2Len >= $g_qMinLen) && ($r2Len <= $lBest1) && ($r1Len <= $lBest1)) {
          my $s1MergedTrimmed = substr($s1Merged, $iBest1, $lBest1);
          my $q1MergedTrimmed = substr($q1Merged, $iBest1, $lBest1);
          my $s2MergedTrimmed = reverse $s1MergedTrimmed;
          my $q2MergedTrimmed = reverse $q1MergedTrimmed;
          $s2MergedTrimmed =~ tr/ACGT/TGCA/;
          $fhOutUnmergedR1->print(join("\n", $h1, substr($s1MergedTrimmed, 0, $r1Len), "+", substr($q1MergedTrimmed, 0, $r1Len),"")) || die ("Could not write to fhOutR1 :$fhOutUnmergedR1:\n");
          $fhOutUnmergedR2->print(join("\n", $h2, substr($s2MergedTrimmed, 0, $r2Len), "+", substr($q2MergedTrimmed, 0, $r2Len),"")) || die ("Could not write to fhOutR1 :$fhOutUnmergedR2:\n");
        }
      }
      $nReadsSE++;
    }
  }
  elsif ((! $junkIt) && (!$mergedReadsIntoSE)) {
    $iBest1 = ($iqBest1 < 2) ? 2 : $iqBest1;
    $lBest1 = ($jqBest1 - $iBest1 + 1);
    $iBest2 = ($iqBest2 < 2) ? 2 : $iqBest2;
    $lBest2 = ($jqBest2 - $iBest2 + 1);
    my $h12=$h1;
    $h12.=" t5:$iBest1 t3:$iBest2";
    if ( ($mtSeq_h->{substr($s1,2,$mtLength)} && $mtSeq_h->{substr($s2,2,$mtLength)})) {
      $fhOutMT->print(join("\t", $h12, substr($s1, $iBest1, $lBest1), substr($s2, $iBest2, $lBest2))."\n");
      $nReadsMT++;
    }
    else {
      $fhOutR1R2->print(join("\t", $h12, substr($s1, $iBest1, $lBest1), substr($s2, $iBest2, $lBest2))."\n");
      $fhOutR1->print(join("\n", $h1, substr($s1, $iBest1, $lBest1), "+", substr($q1,  $iBest1, $lBest1),"")) || die ("Could not write to fhOutR1 :$fhOutR1:\n");
      $fhOutR2->print(join("\n", $h2, substr($s2, $iBest2, $lBest2), "+", substr($q2,  $iBest2, $lBest2),"")) || die ("Could not write to fhOutR1 :$fhOutR2:\n");
      if ($g_writeUnmerged) {
        $fhOutUnmergedR1->print(join("\n", $h1, substr($s1, $iBest1, $lBest1), "+", substr($q1,  $iBest1, $lBest1),"")) || die ("Could not write to fhOutR1 :$fhOutUnmergedR1:\n");
        $fhOutUnmergedR2->print(join("\n", $h2, substr($s2, $iBest2, $lBest2), "+", substr($q2,  $iBest2, $lBest2),"")) || die ("Could not write to fhOutR1 :$fhOutUnmergedR2:\n");
      }
    }
  }
  elsif ($failQual1 || $failQual2) {
    $fhOutFailQR1->print(join("\n", $h1, $s1, "+", $q1,""));
    $fhOutFailQR2->print(join("\n", $h2, $s2, "+", $q2,""));
    $nFailQ1 += $failQual1;
    $nFailQ2 += $failQual2;
    $nFailQ++;
  }
  elsif ($overlap||$failAdapter1 || $failAdapter2) {
    $fhOutFailAdapterR1->print(join("\n", $h1, $s1, "+$failAdapter1", $q1,""));
    $fhOutFailAdapterR2->print(join("\n", $h2, $s2, "+$failAdapter2", $q2,""));
    $nFailAdapter1 += $failAdapter1;
    $nFailAdapter2 += $failAdapter2;
    $nFailAdapter++;
  }
  elsif ($failNs1 || $failNs2) {
    $fhOutFailNsR1->print(join("\n", $h1, $s1, "+", $q1,""));
    $fhOutFailNsR2->print(join("\n", $h2, $s2, "+", $q2,""));
    $nFailNs1 += $failNs1;
    $nFailNs2 += $failNs2;
    $nFailNs++;
  }

  return 1;
}

sub openInZ {
  my $fnIn = shift @_;

  my $fhIn = undef;

  open ($fhIn, "-|", "pigz -d -p 3 --stdout $fnIn");
  if (! $fhIn) {
    print STDERR "pigz Failed to open gzipped file '$fnIn' for reading\n";
    $fhIn = undef;
  }
  else {
    print STDERR "pigz Opened gzipped file '$fnIn' for reading\n";
  }

  return $fhIn;
}

sub openIn {
  my $fnIn = shift @_;

  my $fhIn = undef;

  if ($fnIn =~ m/[.]gz$/) {
    $fhIn = openInZ($fnIn);
  }
  else {
    $fhIn = FileHandle->new();
    if (! $fhIn->open("<$fnIn")) {
      print STDERR "Failed to open file '$fnIn' for reading\n";
      $fhIn = undef;
    }
    else {
      print STDERR "Opened file '$fnIn' for reading\n";
    }
  }

  return $fhIn;
}

sub openOutPigz {
  my $fnOut = shift @_;

  my $fhOut = undef;
  open ($fhOut, "|-", "pigz -p 3 --stdout > $fnOut");

  if (! $fhOut) {
    $fhOut = undef;
    print STDERR "Failed to open file '$fnOut' for pigzed writing\n";
  }
  else {
    print STDERR "Opened file '$fnOut' for pigzed writing\n";
  }
  return $fhOut;
}

sub openOutZ {
  my $fnOut = shift @_;

  return openOutPigz($fnOut);
}

sub openOut {
  my $fnOut = shift @_;

  my $fhOut = undef;

  if ($fnOut =~ m/[.]gz$/) {
    $fhOut = openOutZ($fnOut);
  }
  else {
    $fhOut = FileHandle->new();

    if (! $fhOut->open(">$fnOut")) {
      $fhOut = undef;
      print STDERR "Failed to open file '$fnOut' for writing\n";
    }
    else {
      print STDERR "Opened file '$fnOut' for writing\n";
    }
  }
  return $fhOut;
}
