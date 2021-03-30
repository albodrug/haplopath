#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable;
use GraphViz;

my $file;
my $jaccard = 0.1;
my $help;
GetOptions( 'file=s'   => \$file,
            'jaccard:f'   => \$jaccard,
            'help|h' => \$help,
          );

=head1
Use the script as follows:
perl haplopath_validator.pl --file test_output/link*.bps --jaccard 0.1

	file: 		file containing the binary link output of haplopath_connect.pl
	jaccard: 	limit of jaccard index for a relation to be considered valid (default 0.1)
	help: 		display help

=cut

# Print help
die `pod2text $0` if $help;


my $genomeref = retrieve($file);
my @genome = @$genomeref;
my %region_name_consensus;
my %connections_sum;
my %connections_count;
my %connections_ord;

#making the hash region name -> region consensus
for my $el (@genome){
  my %region = %$el;
  my $region_name = $region{chrom}."_".$region{start}."_".$region{end};
  my $region_consensus = $region{cons};
  $region_name_consensus{$region_name} = $region_consensus;
}

for my $el (@genome){
  my %region       = %$el;
  my $chrom        = $region{chrom};
  my ($tig, $len)  = split /len=/, $chrom;
  my $start        = $region{start};
  my $end          = $region{end};
  my $region_links = $region{link};
  my $cons         = $region{cons};
  my $ord1 = "n";
  $ord1 = "h" if ($start<=$len/2);
  $ord1 = "t" if ($start>$len/2);
  for my $lref (@$region_links){
    my $sumforce;
    my $nbforce;
    my %l              = %$lref;
    my $lchrom         = $l{lchrom};
    if ($chrom ne $lchrom){
      my ($ltig, $llen)  = split /len=/, $lchrom;
      my $lstart         = $l{lstart};
      my $lend           = $l{lend};
      my $ord2 = "a";
      $ord2 = "h" if ($lstart<=$llen/2);
      $ord2 = "t" if ($lstart>$llen/2);
      my $linkjaccard     = $l{jaccardi};
      my @ordered = sort { lc($a) cmp lc($b) } ($chrom, $lchrom);
      
      $connections_count{$ordered[0].'-'.$ordered[1]} += 1 if $linkjaccard>=$jaccard;   
      $connections_sum{$ordered[0].'-'.$ordered[1]} += $linkjaccard if $linkjaccard>=$jaccard;
      $connections_ord{$ordered[0].'-'.$ordered[1]} .= $ord1.$ord2."_";     
      #print $ordered[0], " ", $ordered[1], " ", $linkjaccard, "\n";
    }
  }
}

while(my($tigs, $count) = each %connections_count){
  my $sum = $connections_sum{$tigs};
  my $ord = $connections_ord{$tigs};
  my ($t1, $t2) = split /-/, $tigs;
  #print "Strong relation between $t1 and $t2 : $count relations with an average of $sum/$count Jaccard index (>=0.1).\n";
  my @orders = split /_/, $ord;
  my %ordcount;
  $ordcount{$_}++ foreach @orders;  

  my @array;
  my @array2;

  foreach my $o (sort { $ordcount{$a} <=> $ordcount{$b} } keys %ordcount) {
   push @array, $ordcount{$o};
   push @array2, $o; 
  }
  if ($array[-1] eq $array[-2]){
    $ord='na'
  }else{
    $ord=$array2[-1];
  }

  print "$t1 $t2 ", $sum/$count, " ", $count/2, " ", $ord, "\n";
  #print Dumper(\%ordcount);
}






