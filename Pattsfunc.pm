#!/usr/bin/perl
package Pattsfunc;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('update_reference_pat','get_consensus','get_number_accessions', 'load_data', 'create_blocks', 'find_links', 'join_consensuses', 'reconstruct_blocks', 'pattern_comparison2', 'clean', 'get_subconsensus', 'get_vertical_consensus');
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use List::Util qw( reduce );
use List::Util qw( min max sum);
use Storable;
use File::Basename;
use Parallel::ForkManager qw( );
use threads;
use threads::shared;

sub get_number_accessions{
  my $rarg = shift;
  my %arg = %$rarg;
  my $number = `bcftools query -l $arg{'file'} | wc -l`;
  chomp $number;
  return $number;
}

sub load_data{
  my $rarg = shift;
  my %arg = %$rarg;
  my @patterns;
  my $command = system("bcftools view $arg{'file'} | vcf-annotate -f StrandBias=0.0001 -f EndDistBias=0.0001 -f SnpGap=10 --hard-filter | vcftools --vcf - --maf 0.05 --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | bcftools query -f '%CHROM %POS [%GT ]\n' | sed 's,[/|],,g' | sed 's, ,,3g' > $arg{'file'}.patt");
  open(my $patfh, '<:encoding(UTF-8)', $arg{'file'}.".patt")
    or die "Could not open file $!";
  while (my $string=<$patfh>){
    chomp $string;
    my ($chrom, $coord, $patty) = split / /, $string;
    my $pattern = &transform($patty);
    if (&validpat(\%arg, $pattern)){
      my $abyss = $pattern =~ tr/.//;
      my $heter = $pattern =~ tr/v//;
      $chrom =~ s/_//g;
      push @patterns, {pat => $pattern, coord => $coord, chrom => $chrom, heter=>$heter, abyss=>$abyss};
    }
  }
  close $patfh;
  return @patterns;
}

sub transform{
  my $str = $_[0];
  my $trstr;
  my $len = length($str);
  for (my $i=0; $i<$len; $i=$i+2){
    my $call = substr($str, $i, 2);
    #print "$str\n$call -> $i\n";
    $trstr .= 'V' if $call eq '11';
    $trstr .= 'N' if $call eq '00';
    $trstr .= 'v' if $call eq '01';
    $trstr .= '.' if $call eq '..';
 }
 if (length($trstr) ne length($str)/2){
   print "WARNING, error while transforming string\n";
 }
 return $trstr
}

sub validpat{
  my $rarg = shift;
  my $str = shift;
  my %arg = %$rarg;
  my $pattern = $str;
  my $bool = 0;
  my $number_abyss = $pattern =~ tr/.//;
  my $number_refer = $pattern =~ tr/N//;
  my $number_alter = $pattern =~ tr/V//;
  if ($number_abyss/length($pattern)<=$arg{'maxnoinfo_raw'} && $number_refer >= $arg{'minre'} && $number_alter >= $arg{'minva'}){
    $bool=1;
  }
  return $bool;
}

sub create_blocks{
  my $klist    = shift;
  my $parr     = shift;
  my $rarg     = shift;
  my @chromosome_list = @$klist;
  my @patterns        = @$parr;
  my %arg             = %$rarg;
  #print Dumper('create_blocks',$rarg);
  my @blocky_mcblock;
  for my $chromosome (@chromosome_list){
    my @subpatterns = grep { $_->{chrom} eq $chromosome } @patterns;
    my @coordinate_list = grep defined, map {$_->{coord}} @subpatterns;
    my @sorted_coordinate_list = sort { $a <=> $b } @coordinate_list;
    my $starting_coordinate = $sorted_coordinate_list[0];
    my $base = $arg{'bsize'};
    my $rounded_starting_coordinate = int($starting_coordinate/$base)*$base;
    my $count=0;
    my $pm = Parallel::ForkManager->new($arg{'threads'});
    for (my $i = $rounded_starting_coordinate; $i<= $sorted_coordinate_list[-1]; $i=$i+$arg{'ssize'}){
      $count++;
      my @internal_blocky_mcblock;
      my $pid = $pm->start and next; # PARA
      my $start = $i;
      my $end = $i+$arg{'bsize'};
      my @block;
      for my $coord (@sorted_coordinate_list){
        if ($coord > $start && $coord <= $end){
          push @block, $coord;
        }elsif ($coord < $end){
          next;
        }
      }
      my @block_pattern;
      for my $c (@block){
        my @pat = grep { $_->{coord} eq $c } @subpatterns;
        push @block_pattern, @pat;
      }
      #print "New block consensus built: chr:$chromosome --> blo:$count.\n";
      my %consens = &get_consensus(\@block_pattern, \%arg);
      my %compress_consens = %consens; #&summarize_consensus_hash(\%consens, \%arg);
      push @internal_blocky_mcblock, {chrom => $chromosome, start => $start, end => $end, cons => \%compress_consens};
      store \@internal_blocky_mcblock, 'tmp'.$arg{'outputb'}.'_'.$chromosome.'_'.$count;
      print "\nNew block consensus built: chr:$chromosome --> blo:$count has cns:".scalar(keys %compress_consens)." .\n";
      $pm->finish;
    }
    $pm->wait_all_children;
    my $temporary = 'tmp'.$arg{'outputb'}.'_'.$chromosome.'_';
    my @files = <$temporary*>;
    for my $storage (@files){
      my $story = retrieve($storage);
      my @story = @$story;
      push @blocky_mcblock , @story;
      unlink $storage or warn "Could not unlink $storage: $!";
    }
  }
  my @sorted_blocky_mcblock =  sort { $a->{chrom} cmp $b->{chrom} || $a->{start} <=> $b->{start} } @blocky_mcblock;
  return @sorted_blocky_mcblock;
}

sub get_consensus{
  my $arr    = shift; my @array = @$arr;
  my $rarg    = shift; my %arg = %$rarg;
  # a pattern will make it into the consensus library if its abys satisfies the condition
  #my $perc_abyss_allowed_in_consensus = 0.15;
  my $paaic = $arg{'maxnoinfo'};
  #print Dumper('get_consensus',$rarg);
  #print Dumper(@array);
  my %consensus;
  if (scalar(@array) < $arg{'minco'} ){
    %consensus = ();
  }else{
    my @sorted_patterns =  sort { $b->{abyss}+$b->{heter} <=> $a->{abyss}+ $a->{heter} or $b->{abyss} <=> $a->{abyss} } @array; 
    while (scalar(grep defined, @sorted_patterns)>0){
      @sorted_patterns = grep defined, @sorted_patterns;
      my $i=$#sorted_patterns;
      my $ref = $sorted_patterns[$i]->{pat};
      if ($ref){
        #print ":NEXT REF:\n$ref:REF\n";
        my @matched=();
        for (my $j = $#sorted_patterns; $j>=0; $j--){
          my $qwr = $sorted_patterns[$j]->{pat};
          if ($qwr){
            my ($wlsref, $pmamref) = &pattern_comparison2(\$ref, \$qwr, \%arg);
            my $wls = $$wlsref; my $pmam = $$pmamref;
            if ($wls <= $pmam){
              push @matched, $qwr;
              my $upref = &update_reference_pat(\$ref, \$qwr, \@matched, \%arg);
              #print "$upref:UPD\n";
              $ref = $upref;
              undef $sorted_patterns[$j];
            }
          }
        }
        if (scalar(@matched)>= $arg{'minco'}){
          #print Dumper('matched '.scalar(@matched), @matched);
          my $consensusref = &get_subconsensus(\@matched, \%arg);
          my $cns = $$consensusref; #$$consensusref;
          #print "consensus $cns\n";
          my $count = scalar(@matched);
          my $abyss = $cns =~ tr/.//;
          my $refer = $cns =~ tr/N//;
          my $alter = $cns =~ tr/V//;
          if ($arg{'verbose'}){
            print 'matched and retained '.scalar(@matched)."\n".$cns.":cns\n".join("\n",@matched)."\n" if $abyss < $arg{'nacc'}*$paaic;
            print 'matched and ditchedd '.scalar(@matched)."\n".$cns.":cns\n".join("\n",@matched)."\n" if $abyss >= $arg{'nacc'}*$paaic;
          }
          $consensus{$cns} = $count if ($abyss < $arg{'nacc'}*$paaic && $refer >=$arg{'minre'} && $alter >= $arg{'minva'});
        }else{
          undef $sorted_patterns[$i];
        }
      }
    }
  }
  #print scalar(@array)." and ".scalar(keys %consensus)."\n".Dumper(@array);
  #print Dumper('hey',%consensus);
  return %consensus;
}

sub update_reference_pat{
  my $rref = shift; my $ref = $$rref;
  my $rqwr = shift; my $qwr = $$rqwr;
  my $rmatched = shift; my @matched = @$rmatched; 
  my $rarg = shift; my %arg = %$rarg;
  #print "$ref : newref\n";
  #print "#matched#\n";
  #print join("\n",@matched), "\n";
  #print "#end#\n\n";
  my $rupref = &get_subconsensus(\@matched, \%arg);;
  my $upref = $$rupref;
  #for(0 .. length($qwr)){
  #  my $charref = substr($ref, $_, 1);
  #  my $charqwr = substr($qwr, $_, 1);
  #  if ($charref eq '.' and $charqwr ne '.'){
  #    $upref.=$charqwr; # filled in by query
  #  }elsif($charqwr eq 'v' and $charref ne 'v'){
  #    $upref.=$charqwr; # updated hetero by query
  #  }else{
  #    $upref.=$charref;
  #  }
  #}
  return $upref;
}

sub summarize_consensus_hash{
  my $hashref = shift;
  my %hash = %$hashref;
  my $rarg    = shift; my %arg = %$rarg;
  my %summary_hash;
  #print Dumper('summarize_consensus_hash',$rarg);
  my @patterns = keys(%hash);
  while (scalar(grep defined, @patterns)>0){
    @patterns = grep defined, @patterns;
    my $i=$#patterns;
    my $ref = $patterns[$i];
    if ($ref){
      my @matched=();
      for (my $j = $#patterns; $j>=0; $j--){
        my $qwr = $patterns[$j];
        if ($qwr){
          my ($wlsref, $pmamref) = &pattern_comparison_consensus(\$ref, \$qwr, \%arg);
          my $wls = $$wlsref; my $pmam = $$pmamref;
          if ($wls <= $pmam){
            my $upref = &update_reference_pat(\$ref, \$qwr, \@patterns, \%arg);
            $ref = $upref; 
            push @matched, $qwr;
            undef $patterns[$j];
          }
        }
      }
      #my $consensusref = &get_subconsensus_summary(\@matched, \%arg);
      my $cns = $ref;#$$consensusref;
      my $real_occurence = 0;
      for my $el (@matched){
        $real_occurence += $hash{$el};
      }
      my $cnt_abys   = $cns =~ tr/.//;
      my $freq_abys  = $cnt_abys / $arg{'nacc'};
      $summary_hash{$cns} = $real_occurence if ($freq_abys < $arg{'maxnoinfo'}); # c mega cho
    }
  }
  return %summary_hash;
}

sub pattern_comparison2 {
  my $refref = shift;
  my $qwrref = shift;
  my $ref = ${$refref};
  my $qwr = ${$qwrref};
  my $rarg    = shift; my %arg = %$rarg;
  my $wls=0; # weighted levenshtein score
  my $emptyness=0;
  my $pmam=0; # pondered maximum allowed mismatches
  for(0 .. length($qwr)){
    my $charref = substr($ref, $_, 1);
    my $charqwr = substr($qwr, $_, 1);
    if ($charref =~ m/[xo.]/ or $charqwr =~ m/[xo.]/){
      $emptyness = $emptyness + 1;
    }
  }
  my $cnt_N      = $qwr =~ tr/N//;
  my $cnt_V      = $qwr =~ tr/V//;
  my $minLET = min($cnt_N, $cnt_V);
  my $homoz_mismatch_penalty = 4000;
  $pmam = ($arg{'nacc'}-($emptyness))*$arg{'diff'};
  if ($minLET<=$arg{'nacc'}*0.20){ # In the case of 50 accessions, i want strict homo penalty # patterns with less than 10 Ns or 10 Vs
    $homoz_mismatch_penalty = $pmam +1;
  }else{
    $homoz_mismatch_penalty = $pmam / $arg{'homi'} if $arg{'homi'} ne 0;
    $homoz_mismatch_penalty = $pmam +1 if $arg{'homi'} eq 0;
  }
  my $hetez_mismatch_penalty = $pmam/($arg{'nacc'}*0.25);#$homoz_mismatch_penalty = 0 if $homoz_mismatch_penalty<0 ;
  
  for(0 .. length($qwr)){
    my $charref = substr($ref, $_, 1);
    my $charqwr = substr($qwr, $_, 1);
    my $combi = $charref . $charqwr;
    if ($combi eq 'NV' || $combi eq 'VN' ){ # mismatch homo
      $wls = $wls + $homoz_mismatch_penalty;              # deal breaker
      if ($wls > $pmam){
        last;
      }
    }elsif($combi eq 'NW' || $combi eq 'VM' ||$combi eq 'MV' ||$combi eq 'WN'){ #mismatch semihomo
      $wls = $wls + $pmam; # danger
      if ($wls > $pmam){
        last;
      }
    }elsif ($combi eq 'Nv' || $combi eq 'Vv' || $combi eq 'vN' ||$combi eq 'vV' || $combi eq 'MW' || $combi eq 'WM'){ # mismatch hetero
      $wls=$wls + $hetez_mismatch_penalty; #allowed differences, pmam is calculated for hetero differences
      if ($wls > $pmam){
        last;
      }
    }elsif($charref eq '.' or $charqwr eq '.'){
      $wls=$wls + 0 ; # if there is an empty spot, I can't compare it, I take a risk, wls takes a hit
      if ($wls > $pmam){
        last;
      }
    }
  }
  return \$wls, \$pmam;
}

sub pattern_comparison_consensus {
  my $refref = shift;
  my $qwrref = shift;
  my $ref = ${$refref};
  my $qwr = ${$qwrref};
  my $rarg    = shift; my %arg = %$rarg;
  my $wls=0; # weighted levenshtein score
  my $emptyness=0;
  my $pmam=0; # pondered maximum allowed mismatches
  for(0 .. length($qwr)){
    my $charref = substr($ref, $_, 1);
    my $charqwr = substr($qwr, $_, 1);
    if ($charref =~ m/[xo.]/ or $charqwr =~ m/[xo.]/){
      $emptyness = $emptyness + 1;
    }
  }
  $pmam = ($arg{'nacc'}-($emptyness))*$arg{'diff'};
  my $hetez_mismatch_penalty = $pmam/($arg{'nacc'}*0.25);
  my $homoz_mismatch_penalty = $pmam +1;
  for(0 .. length($qwr)){
    my $charref = substr($ref, $_, 1);
    my $charqwr = substr($qwr, $_, 1);
    my $combi = $charref . $charqwr;
    if ($combi eq 'NV' || $combi eq 'VN' ){ # mismatch homo
      $wls = $wls + $homoz_mismatch_penalty;              # deal breaker
      if ($wls > $pmam){
        last;
      }
    }elsif($combi eq 'NW' || $combi eq 'VM' ||$combi eq 'MV' ||$combi eq 'WN'){ #mismatch semihomo
      $wls = $wls + $homoz_mismatch_penalty; 
      if ($wls > $pmam){
        last;
      }
    }elsif ($combi eq 'Nv' || $combi eq 'Vv' || $combi eq 'vN' ||$combi eq 'vV' || $combi eq 'MW' || $combi eq 'WM'){ # mismatch hetero
      $wls=$wls + $hetez_mismatch_penalty; #allowed differences, pmam is calculated for hetero differences
      if ($wls > $pmam){
        last;
      }
    }elsif($charref eq '.' or $charqwr eq '.'){
      $wls=$wls + 0 ; # if there is an empty spot, I can't compare it, I take a risk, wls takes a hit
      if ($wls > $pmam){
        last;
      }
    }
  }
  return \$wls, \$pmam;
}

sub get_subconsensus{
  my $arr = shift; my @array = @$arr;
  my $rarg = shift; my %arg = %$rarg;
  my $consensus;
  my $len = $arg{'nacc'};
  #print "$len\n";
  for (0 .. $len-1){
    my $consstring;
    for my $pattern (@array){
      my $patchar = substr($pattern, $_, 1);
      $consstring .= $patchar;
    }
    my $conscharref = &get_vertical_consensus(\$consstring, \%arg);
    my $conschar = $$conscharref;
    $consensus .= $conschar ;
    #print "$consstring --> $conschar\n";
  }
  return \$consensus;
}

sub get_subconsensus_summary{
  my $arr = shift; my @array = @$arr;
  my $rarg = shift; my %arg = %$rarg;
  my $consensus;
  my $len = $arg{'nacc'};
  for (0 .. $len-1){
    my $consstring;
    for my $pattern (@array){
      my $patchar = substr($pattern, $_, 1);
      $consstring .= $patchar;
    }
    my $conscharref = &get_vertical_consensus_summary(\$consstring, \%arg);
    my $conschar = $$conscharref;
    $consensus .= $conschar ;
  }
  return \$consensus;
}

sub get_vertical_consensus_summary {
  my $rstring = shift; my $string = $$rstring;
  my $rarg    = shift; my %arg = %$rarg;
  my $consensus_char = '?';
  my $cnt_N      = $string =~ tr/N//;
  my $cnt_V      = $string =~ tr/V//;
  my $cnt_M      = $string =~ tr/M//;
  my $cnt_W      = $string =~ tr/W//;
  my $cnt_v      = $string =~ tr/v//;
  my $cnt_abys   = $string =~ tr/.//;
  my $tot        = $string =~ tr/NVvMW.//;
  my $info_tot   = $string =~ tr/NVvMW//;
  my %freqhash;
  if ($info_tot != 0){
    %freqhash = (
         "N"  => $cnt_N/$tot,
         "V"  => $cnt_V/$tot,
         "M"  => $cnt_M/$tot,
         "W"  => $cnt_W/$tot,
         "v"  => $cnt_v/$tot,
         "."  => $cnt_abys/$tot,
    );
  }else{
     %freqhash = (
         "N"  => 0,
         "V"  => 0,
         "M"  => 0,
         "W"  => 0,
         "v"  => 0,
         "."  => 1,
    );
  }
  my $max_char = reduce { $freqhash{$a} > $freqhash{$b} ? $a : $b } keys %freqhash;
  if ($max_char ne '.' and $freqhash{$max_char}>=$arg{'consthres'}){
    $consensus_char = $max_char;
  }else{
    if ( $string !~ /N/ and $string =~ /[Vv]/ and ($cnt_V+$cnt_v)>=1){
      if ($string =~ /[v]/ and $freqhash{'v'}>=0.3){
        $consensus_char="v";
      }else{
        $consensus_char="V";
      }
    }elsif( $string !~ /V/ and $string =~ /[Nv]/ and ($cnt_N+$cnt_v)>=1){
       if ($string =~ /[v]/ and $freqhash{'v'}>=0.3){
         $consensus_char="v";
       }else{
         $consensus_char="N";
       }
    }elsif( $string =~ /V/ and $string =~ /N/ and $string =~ /v/ and ($freqhash{'V'}+$freqhash{'N'}+$freqhash{'v'})>=0.40){
      if (($freqhash{'v'})>=0.30){ # one third of positions
        $consensus_char="v";
      }else{
        $consensus_char="v"; # a decision can not be made because it is an even mix from V, N and v.
        #print "$string -> $freqhash{'V'} $freqhash{'N'} $freqhash{'v'}\n"; 
      }
    }elsif($string =~ /V/ and $string =~ /N/ and ($freqhash{'V'}+$freqhash{'N'})>=0.4){
      $consensus_char="v";
    }else{
      $consensus_char=".";
      #print "h --> $string\n";
    }
  }
 #print "Hey: $string $max_char $consensus_char\n";
 return \$consensus_char;
}

sub get_vertical_consensus {
  my $rstring = shift; my $string = $$rstring;
  my $rarg    = shift; my %arg = %$rarg;
  my $consensus_char = '?';
  my $cnt_N      = $string =~ tr/N//;
  my $cnt_V      = $string =~ tr/V//;
  my $cnt_M      = $string =~ tr/M//;
  my $cnt_W      = $string =~ tr/W//;
  my $cnt_v      = $string =~ tr/v//;
  my $cnt_abys   = $string =~ tr/.//;
  my $tot        = $string =~ tr/NVvMW.//;
  my $info_tot   = $string =~ tr/NVvMW//;
  my %freqhash;
  if ($info_tot != 0){
    %freqhash = (
         "N"  => $cnt_N/$tot,
         "V"  => $cnt_V/$tot,
         "M"  => $cnt_M/$tot,
         "W"  => $cnt_W/$tot,
         "v"  => $cnt_v/$tot,
         "."  => $cnt_abys/$tot,
    );
  }else{
     %freqhash = (
         "N"  => 0,
         "V"  => 0,
         "M"  => 0,
         "W"  => 0,
         "v"  => 0,
         "."  => 1,
    );
  }
  my $max_char = reduce { $freqhash{$a} > $freqhash{$b} ? $a : $b } keys %freqhash;
  if ($max_char ne '.' and $freqhash{$max_char}>=$arg{'consthres'}){
    $consensus_char = $max_char;
  }else{
    if ( $string !~ /N/ and $string =~ /[Vv]/ and (($cnt_V+$cnt_v)>=3 or $freqhash{'V'}+$freqhash{'v'}>=0.4)){
      if ($string =~ /[v]/ and $freqhash{'v'}>=0.3){
        $consensus_char="v";
      }else{
        $consensus_char="V";
      }
    }elsif( $string !~ /V/ and $string =~ /[Nv]/ and (($cnt_N+$cnt_v)>=3 or $freqhash{'N'}+$freqhash{'v'}>=0.4)){
       if ($string =~ /[v]/ and $freqhash{'v'}>=0.3){
         $consensus_char="v";
       }else{
         $consensus_char="N";
       }
    }elsif( $string =~ /V/ and $string =~ /N/ and ($freqhash{'V'}+$freqhash{'N'}+$freqhash{'v'})>=0.40){
      if ($string =~ /v/){ # one third of positions
        $consensus_char="v";
      }else{
        if ($cnt_V > $cnt_N){
          $consensus_char="V";
        }elsif($cnt_V < $cnt_N){
          $consensus_char="N";
        }else{
          $consensus_char=".";
        }
          
      }
    }else{
      $consensus_char=".";
      #print "echec . --> $string $freqhash{'N'} \n";
    }
  }
  
  #my %mix_hashocc;
  #$mix_hashocc{'v'} = ($hashocc{'W'}+$hashocc{'V'}+$hashocc{'v'});
  #$mix_hashocc{'v'} = ($hashocc{'M'}+$hashocc{'N'}+$hashocc{'v'});
  #$mix_hashocc{'.'} = ($hashocc{'N'}+$hashocc{'V'});

  #my $max_char        = reduce { $hashocc{$a} > $hashocc{$b} ? $a : $b } keys %hashocc;
  #my $max_char_mix    = reduce { $mix_hashocc{$a} > $mix_hashocc{$b} ? $a : $b } keys %mix_hashocc;

  #if ($hashocc{$max_char}>=$arg{'consthres'}){
  #  $consensus_char=$max_char;
  #}elsif ($mix_hashocc{$max_char_mix}>=$arg{'consthres'} && $hashocc{$max_char}<$arg{'consthres'}){
  #  $consensus_char=$max_char_mix;
  #}elsif (substr( $string, 0, 1 ) eq 'v'){
  #   $consensus_char="v";
  #}else{
  #  if ( $string !~ /N/ and $string =~ /[Vv]/ ){
  #    $consensus_char="v";
  #  }elsif( $string !~ /V/ and $string =~ /[Nv]/ ){
  #     $consensus_char="v";
  #  }else{
  #    $consensus_char=".";
  #  }
 #}
 #print "$string --> $consensus_char : $arg{'consthres'} $freqhash{$max_char} \n";
 return \$consensus_char;
}

sub get_occurence {
  my $string = $_[0];
  my $cnt_N      = $string =~ tr/N//;
  my $cnt_V      = $string =~ tr/V//;
  my $cnt_M      = $string =~ tr/M//;
  my $cnt_W      = $string =~ tr/W//;
  my $cnt_v      = $string =~ tr/v//;
  my $cnt_abys   = $string =~ tr/.//;
  my $tot        = $string =~ tr/NVvMW.//;
  my $info_tot   = $string =~ tr/NVvMW//;
  my %hsh_frq_of;
  if ($info_tot != 0){
    %hsh_frq_of = (
         "N"  => $cnt_N/$tot,
         "V"  => $cnt_V/$tot,
         "M"  => $cnt_M/$tot,
         "W"  => $cnt_W/$tot,
         "v"  => $cnt_v/$tot,
         "."  => $cnt_abys/$tot,
    );
  }else{
     %hsh_frq_of = (
         "N"  => 0,
         "V"  => 0,
         "M"  => 0,
         "W"  => 0,
         "v"  => 0,
         "."  => 1,
    );
  }
  return \%hsh_frq_of;
}

sub find_links{
  my $block = shift;
  my $rarg  = shift;
  my @blocky_mcblock = @$block;
  my %arg = %$rarg;
  print "Finding links...\n";
  #print Dumper(\@blocky_mcblock);
  my $mishp = $arg{'mishp'}; # minimum shared patterns
  my $mioccshp = $arg{'mioccshp'}; # minimum occurence of shared patterns
  my @blocky_mcblock_mclink;
  my $pm = Parallel::ForkManager->new($arg{'threads'});
  for my $index01 (0..$#blocky_mcblock) {
    print "Looking into block $blocky_mcblock[$index01]->{chrom}_$blocky_mcblock[$index01]->{start}:$blocky_mcblock[$index01]->{end}\n";
    my $pid = $pm->start and next; #PARA
    my $blo1r = $blocky_mcblock[$index01];
    my @block_linkage;
    for my $index02 (0..$#blocky_mcblock) {
      my $blo2r = $blocky_mcblock[$index02];
      my @occurence = ();
      my $jaccardi= -2;
      if ($blo1r ne $blo2r){
        my $blo1_chr   = $blo1r->{chrom};
        my $blo2_chr   = $blo2r->{chrom};
        my $blo1_start = $blo1r->{start};
        my $blo2_start = $blo2r->{start};
        my $blo1_end = $blo1r->{end};
        my $blo2_end = $blo2r->{end};
        if ($blo1_chr ne $blo2_chr){ # different chromosomes
          #print "$blocky_mcblock[$index01]->{chrom}_$blocky_mcblock[$index01]->{start}:$blocky_mcblock[$index01]->{end} compared to $blocky_mcblock[$index02]->{chrom}_$blocky_mcblock[$index02]->{start}:$blocky_mcblock[$index02]->{end}\n";
          my $blo1_consensus_hash = $blo1r->{cons};
          my $blo2_consensus_hash = $blo2r->{cons};
          my ($occyref, $jaccardiref) = &get_link_count($blo1_consensus_hash, $blo2_consensus_hash, \%arg);
          my %occy = %$occyref;
          my $count = keys %occy;
          $jaccardi= $$jaccardiref;
          if ($count >= $mishp){ # more than 1 patterns link the two blocks
            push @occurence, \%occy;
          }
        }else{ # within same chromosome
          if (($blo2_start < $blo1_end and $blo2_start > $blo1_start) or ($blo2_end < $blo1_end and $blo2_end > $blo1_start) or ($blo1_start<$blo2_end and $blo1_start>$blo2_end) or ($blo1_end<$blo2_end and $blo1_end>$blo2_end)){
            #do nothing
          }else{
            my $blo1_consensus_hash = $blo1r->{cons};
            my $blo2_consensus_hash = $blo2r->{cons};
            my ($occyref, $jaccardiref) = &get_link_count($blo1_consensus_hash, $blo2_consensus_hash, \%arg);
            $jaccardi= $$jaccardiref;
            my %occy = %$occyref;
            #print Dumper("occy",%occy);
            #print Dumper("jacc",$jaccardi);
            my $count = keys %occy;
            if ($count >= $mishp){ # more than 1 patterns links the block
              push @occurence, \%occy;
            }
          }
        }
        push @block_linkage, {lchrom => $blo2_chr, lstart => $blo2_start,lend => $blo2_end, lforce=>\@occurence, jaccardi=>$jaccardi} if (scalar(@occurence)>0);
      }# no else
    }
    my @single_blocky_mcblock_mclink;
    push @single_blocky_mcblock_mclink, {chrom => $blo1r->{chrom}, start => $blo1r->{start}, end => $blo1r->{end}, cons => $blo1r->{cons}, link => \@block_linkage, linkocc=>scalar(@block_linkage)};
    store \@single_blocky_mcblock_mclink, 'tmplink'.$arg{'outputb'}.'_'.$index01;
    print "New linkage built: $blocky_mcblock[$index01]->{chrom} --> blo:$index01 of $#blocky_mcblock has ".scalar(@block_linkage)." links.\n". Dumper(@single_blocky_mcblock_mclink) . "\n";
    $pm->finish;  
  }
  
  $pm->wait_all_children;
  my $temporary = 'tmplink'.$arg{'outputb'}.'_';
  my @tmpfiles = <$temporary*>;
  for my $storage (@tmpfiles){
    my $story = retrieve($storage);
    my @story = @$story;
    push @blocky_mcblock_mclink , @story;
    unlink $storage or die "Could not delete the temporary link files\n";
  } 
  return @blocky_mcblock_mclink; 
}

sub get_link_count{
  my $count = 0;
  my %occy ;
  my $refref = shift;
  my $qwrref = shift;
  my $rarg   = shift ; my %arg = %$rarg;
  my %blo1_consensus_hash = %$refref;
  my %blo2_consensus_hash = %$qwrref;
  # I am sorting the consensus hash by occurence so the first patterns to be compared are those with higher occurence
  # In this way I avoid hitting patterns with low count if a very similar one but with higher count is present
  # As the higher count will be iterated over first
  my $ref_patcount =0 ;
  my $qwr_patcount = 0;
  my $linked_ref_patcount = 0;
  my $linked_qwr_patcount = 0;
  my @blo1_consensus = sort { $blo1_consensus_hash{$b} <=> $blo1_consensus_hash{$a} } keys %blo1_consensus_hash;
  my @blo2_consensus = sort { $blo2_consensus_hash{$b} <=> $blo2_consensus_hash{$a} } keys %blo2_consensus_hash;
  for my $index1 (0..$#blo1_consensus){
    my $cns1 = $blo1_consensus[$index1];
    my $occ1 = $blo1_consensus_hash{$cns1};
    $ref_patcount+=$occ1; # counting the total amount of pattern occurences in reference blocks
    for my $index2 (0..$#blo2_consensus){
      my $cns2 = $blo2_consensus[$index2];
      my $occ2 = $blo2_consensus_hash{$cns2};
      $qwr_patcount+=$occ2 if ($index1 eq 0); # same in query block
      if ($occ1 >= $arg{'minco'} && $occ2 >= $arg{'minco'}) {
        my ($wlsref, $pmamref) = &pattern_comparison2(\$cns1, \$cns2, \%arg);
        my $wls = $$wlsref; my $pmam = $$pmamref;
        if ($wls <= $pmam){
          #$blo1_consensus[$index1] = 'itiswednesdaymydudes'; # you can not use the same pattern
          #$blo2_consensus[$index2] = 'itiswednesdaymydudes'; # to link two blocks twice
          $blo1_consensus_hash{$cns1} = -1;
          $blo2_consensus_hash{$cns2} = -1; 
          $linked_ref_patcount+=$occ1; # shared patterns in reference block
          $linked_qwr_patcount+=$occ2; # shared patterns in query block
          $count = $count + 1;
          print Dumper($refref);
          print Dumper($qwrref);
          print "ref:$cns1 :: $occ1\nqwr:$cns2 :: $occ2 $wls $pmam PASSED\n";
          my $key = $cns1."::".$occ1;
          my $value = $cns2."::".$occ2;
          $occy{$key} = $value;
        }else{
          #print "ref:$cns1 :: $occ1\nqwr:$cns2 :: $occ2 $wls $pmam NOT\n";
        }
      }
    }
  }
  my $jaccardi = -1;
  if ($linked_ref_patcount ne 0 && $linked_qwr_patcount ne 0){
    my $normcoef = $linked_ref_patcount / $linked_qwr_patcount;
    my $intersection = $linked_ref_patcount ;
    my $union = $ref_patcount + $qwr_patcount*$normcoef - $intersection;
    $jaccardi = $intersection/$union;
    print "ref:$ref_patcount\nqwr:$qwr_patcount\nlref:$linked_ref_patcount\nlqwr:$linked_qwr_patcount\nc:$normcoef\nunion:$union\ninter:$intersection\njaccardi:$jaccardi\n"; 
  }
  return \%occy, \$jaccardi;
}

sub join_consensuses{
  my $arrref = $_[0];
  my @array  = @$arrref;
  my @allc;
  for my $block (@array){
    my @cons = $block->{cons};
    my $hash = $cons[0];
    my @consensuses = keys %$hash;
    push @allc, @consensuses;
  }
  #print "There are ", scalar(@allc), " inital consensus patterns in ", scalar(@array), " blocks.\n";
  return @allc;
}

sub reconstruct_blocks{
  my $bloref = shift;
  my $curref = shift;
  my @blocks = @$bloref;
  my @newblocks;
  my @currat = @$curref;
  my %currat = map { $_ => 1 } @currat;
  #print Dumper(%currat);
  for my $bloky (@blocks){
    my $chrom  = $bloky->{chrom};
    my $start  = $bloky->{start};
    my $end  = $bloky->{end};
    my @cons  = $bloky->{cons};
    my $hash = $cons[0];
    my %cons = %$hash;
    my %newcons;
    for my $el (keys %cons){
      $newcons{$el} = $cons{$el} if (exists($currat{$el}));
    }
    push @newblocks, {chrom => $chrom, start => $start, end => $end, cons => \%newcons};
  }
  return @newblocks;
}

sub clean{
  my $blocktmp = shift;
  my @block = @$blocktmp;
  my $rarg   = shift ; 
  my %arg = %$rarg;

  my @consensus = &join_consensuses(\@block);
  print "There are ", scalar(@consensus), " inital consensus patterns in ", scalar(@block), " blocks.\n";
  my $pm = Parallel::ForkManager->new($arg{'threads'});

  # What I should rather do is paste all the patterns and compare them without keeping track of where they are
  for (my $j = 0 ;  $j < $#block ; $j ++){
    my $pid = $pm->start and next;
    #print Dumper($block[$j]);
    my @futurehash;
    my $eblock = $block[$j];
    my @cons = $eblock->{cons};
    my $hash = $cons[0];
    my @cons_block = keys %$hash;
    for (my $i = 0 ;  $i < $#cons_block ; $i ++){
      my $ref = $cons_block[$i];
      my $close_scores = 0;
      for (my $k = 0; $k < $#consensus; $k++){
        my $qwr = $consensus[$k];
        my ($wlsref, $pmamref) = &pattern_comparison2(\$ref, \$qwr, \%arg);
        my $wls = $$wlsref; my $pmam = $$pmamref;
        #print "co:\n$ref\nto:\n$qwr $wls $pmam\n";
        if ($wls <= $pmam){
          $close_scores++;
          print "clean match:\n$ref\n$qwr $wls $pmam\n";
          #last to go faster
        }
      }
      push @futurehash, "$eblock->{chrom}_$eblock->{start}_$eblock->{end}_${ref}_$close_scores";
      print "$eblock->{chrom} $eblock->{start} $eblock->{end}\t$ref $close_scores\n";
    }
    print "NEXT BLOCK\n";
    store \@futurehash, 'tmp_curation'.$arg{'outputb'}.'_'.$j;
    $pm->finish;
  }
  $pm->wait_all_children();

  my %occ_hash;
  my @occ_info;
  my $temporary = 'tmp_curation'.$arg{'outputb'}; #print "hey\n";
  my @tmpfiles = <$temporary*>;
  for my $storage (@tmpfiles){
    my $story = retrieve($storage);
    my @story = @$story; #print $story[0]," --> ", $story[1], "\n";
    for my $el (@story){
      my ($chrom, $start, $end, $ref, $count) = split/_/, $el;
      $occ_hash{$ref} = $count;
      push @occ_info, "$chrom $start $end $ref $count";
    }
    unlink $storage or die "Could not delete the temporary link file $storage\n";
  }

  my @curated_consensus = grep{ $occ_hash{$_} <= $arg{'maxwgocc'} && $occ_hash{$_} >= $arg{'minwgocc'} } keys %occ_hash;
  my @dumpster_consensus = grep{ $occ_hash{$_} > $arg{'maxwgocc'} || $occ_hash{$_}  < $arg{'minwgocc'} } keys %occ_hash;

  my @curated_blocks = &reconstruct_blocks(\@block, \@curated_consensus);
  my @removed_blocks = &reconstruct_blocks(\@block, \@dumpster_consensus);
  my @curated_blocky_mcblock = @curated_blocks;

  store \@curated_blocks, "curated_dp".$arg{'outputb'};
  store \@removed_blocks, "removed_dp".$arg{'outputb'};

  my $dumpcurated = Dumper(@curated_blocks);
  my $dumpremoved = Dumper(@removed_blocks);
  my $dumpoccinfo = Dumper(@occ_info);
  open(my $fh, '>', "curated_dp".$arg{'outputh'});
  print $fh $dumpcurated;
  close $fh;
  open(my $fr, '>', "removed_dp".$arg{'outputh'});
  print $fr $dumpremoved;
  close $fr;
  open(my $fg, '>', "removed_dp".$arg{'outputh'}.'occurencehash');
  print $fg $dumpoccinfo;
  close $fg;

  my $curated = scalar(&join_consensuses(\@curated_blocks));
  my $removed = scalar(&join_consensuses(\@removed_blocks));
  print "After cleaning step $curated patterns were kept and $removed patterns were removed.\n";

  return \@curated_blocky_mcblock;

}

1;
