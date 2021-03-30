#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use List::Util qw( reduce );
use List::Util qw( min max sum);
use GD::Simple;
use Storable;
use GraphViz;
use Math::Trig;
use File::Basename;
use POSIX;

=head1
Use as following:


perl haplopath_viz.pl  --file test_output/validator.out --jaccard 0.1 --exclude exclude.list --start qpac010tig00005444len=60711

Usage:


	file 		file containing output of haplopath_validator.pl
	jaccard		limit Jaccard index value to consider relations
	exclude 	file containing list of contig ID to exclude from graph
	start 		ID of starting contig for graph building
        help            display help
=cut


my $file;
my $start;
my $j_thres;
my $exclude;
my @excont;
my $help;
GetOptions( 'file=s'   => \$file,
            'start=s'  => \$start,
            'jaccard=f' => \$j_thres,
            'exclude=s' => \$exclude,
            'help|h' =>\$help,
);


# Print help
die `pod2text $0` if $help;

my $g = GraphViz->new(layout=>'neato',
                        ratio=>"fill");

my $exclfile = $exclude;
open(my $fe, '<:encoding(UTF-8)', $exclfile)
  or die "Could not open file '$exclfile' $!";
while (my $row = <$fe>) {
  chomp $row;
  my $tig = $row;
  push @excont, $row;
}
close $fe;
my %excluded = map { $_ => 1 } @excont;

my $filename = $file;
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";

my @contigs = ($start);
my @connections;
while (my $row = <$fh>) {
  chomp $row;
  my ($tig1, $tig2, $j, $c, $ori) = split(" ", $row);
  if ($j>=$j_thres){
    push @connections, [$tig1, $tig2, $ori, $j];
  }
}
close $fh;

my $count = 0;
my @contigs_linked_to_start = &ctg_lookup(\@contigs, \@connections);
my @connections_linked_to_start = &conn_lookup(\@contigs_linked_to_start, \@connections);

print "\n#",join("\n#",@contigs_linked_to_start),"\n";
#print Dumper(\@connections);

&build_graph(\@contigs_linked_to_start,\@connections_linked_to_start);
print $g->as_text;

sub ctg_lookup{
  my $arr1 = shift; my @contigs = @$arr1;
  my $arr2 = shift; my @connections = @$arr2;
  my @contigs_updated = @contigs;
  #print Dumper(@contigs); 
  foreach my $row (@connections){
    my %presence = map { $_ => 1 } @contigs_updated;
    my $init_contigs_nb = scalar(@contigs_updated);
    my $tig1 = @$row[0];
    my $tig2 = @$row[1];
    if (exists($excluded{$tig1}) or exists($excluded{$tig2})){
      #do noting
    }else{
    if (exists($presence{$tig1}) or exists($presence{$tig2})){
      push @contigs_updated, $tig1 if exists($presence{$tig2});
      push @contigs_updated, $tig2 if exists($presence{$tig1});
      my @u = uniq @contigs_updated; 
      @contigs_updated = @u;
      my $upda_contigs_nb = scalar(@contigs_updated);
      if ($upda_contigs_nb>$init_contigs_nb && $count < 1000000){
        #print Dumper($upda_contigs_nb, $count);
        $count++;
        @contigs_updated = &ctg_lookup(\@contigs_updated, \@connections);
      }elsif ($count>=1000000){
        print "#Deep_recursion\t";
        #last;
      }
    }
  }
  }
  return @contigs_updated;
}

sub conn_lookup{
  my $arr1 = shift; my @contigs = @$arr1;
  my $arr2 = shift; my @connections = @$arr2;
  my @connections_subset;
  my %presence = map { $_ => 1 } @contigs;
  foreach my $row (@connections){
    my $tig1 = @$row[0];
    my $tig2 = @$row[1];
    #print $tig1, $tig2, @$row[2], @$row[3], "\n";
    if (exists($presence{$tig1}) or exists($presence{$tig2})){
      #print Dumper($row);
      push @connections_subset, $row;
    }
  }
  return @connections_subset;
}

sub build_graph{
  my $arr1 = shift; my @contigs = @$arr1;
  my $arr2 = shift; my @connections = @$arr2;
  for my $el (@contigs){
    my ($name, $len) = split("len=", $el);
    my $klen = $len/1000;
    #print "node: $el\n";
    my $node = $el;
    my $node_color = 'cyan';
    $node_color = 'violet' if $node =~ m/lg/;
    $g->add_node($el, shape=>"record", label=>"{$name|${klen}kbp}", fillcolor=>$node_color, style=>"filled");
  }
  for my $tig1 (@contigs){
    for my $tig2 (@contigs){
      my @order;
      my $force = 0;
      my $jaccard = 0;
      for my $conn (@connections){
        if ($tig1 eq @$conn[0] && $tig2 eq @$conn[1]){
          $force++;
          push @order, @$conn[2];
          $jaccard+=@$conn[3];
        }
      }
      if ($force ne 0){
        my $mean_jaccard = sprintf "%.1f", $jaccard / $force;
        my $mostfrequentori = 'na';
        if (scalar(uniq @order) eq 1){
          $mostfrequentori = $order[0];
        }else{
          my %counter;
          foreach my $element (@order){
            $counter{$element}++;
          }
          my @values = values %counter;
          my $max = max(@values);
          my $occurence_of_max = grep { $_ eq $max  } @values;
          if ($occurence_of_max eq 1){
            $mostfrequentori = (sort {$counter{$b} <=> $counter{$a}} @order)[0];
          }
        }
        #print Dumper(\@order);
        #print Dumper(\%count);
        #print "hey $max_value_1, $max_value_2 \n";
        
        my $color = "black";
        $color = "red" if $mostfrequentori eq "hh";
        $color = "blue" if $mostfrequentori eq "tt";
        $color = "violet" if $mostfrequentori eq "ht";
        $color = "cyan4" if $mostfrequentori eq "th";
        print "#edge: $tig1 $tig2 $mean_jaccard $force $mostfrequentori\n";
        my $actual_force = ${force}/2;
        my $pen = 0;
        if ($actual_force*${mean_jaccard}>15){
          $pen=20;
        }else{
          $pen=$actual_force*${mean_jaccard};
        }
        if ($color ne 'black'){
          #label=>"${mean_jaccard};${actual_force}"
          $g->add_edge($tig1 => $tig2, label=>"", color=>$color, penwidth=>$pen);
        }else{
          $g->add_edge($tig1 => $tig2, label=>"", color=>$color, penwidth=>$pen, dir=>'none');
        }
      }
    }
  }
 
  #for my $row (@connections){
  #    my $tig1 = @$row[0]; 
  #    my $tig2 = @$row[1];
  #    my $o = @$row[2];
  #    my $j = sprintf "%.1f", @$row[3];
  #  
  #    my $color = "black";
  #    $color = "red" if $o eq "hh";
  #    $color = "blue" if $o eq "tt";
  #    $color = "violet" if $o eq "ht";
  #    $color = "cyan4" if $o eq "th";
  #    print "#edge: $tig1 $tig2 $j $o\n";
  #    $g->add_edge($tig1 => $tig2, label=>"$j", color=>$color);
  #}
}






























