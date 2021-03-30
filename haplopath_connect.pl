#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use lib 'code_varpat/';
use Pattsfunc;
use List::MoreUtils qw(uniq);
use Data::Dumper;
use Storable;

=head1 
varcon_wrapper.pl
  Perl script that uses variation information from a vcf file to create 
  connections based on shared variation patterns between genomic regions.
Usage
  buildblock|bb 		: boolean, build blocks (default 0, no)
  cleanblock|cb 		: boolean, clean blocks (default 0, no)
  buildlink|bl  		: boolean, build connections (default 0, no)
  help|h        		: boolean, display help
  verbose|v     		: boolean, be verbose
  file|f        		: string, path to file (vcf if starting from --bb, binary storable if starting from --cb or --bl)
  bsize|bs      		: integer, size of genomic segments to build
  ssize|ss      		: integer, sliding window for building the genomic segments
  minva|mv      		: integer, minimum # of accessions carrying an alternative allele
  minre|mr      		: integer, minimum # of accessions carrying a reference-like allele
  minco|mc      		: integer, minimum # of positions carrying the same variation pattern to create fingerprint
  diff|d        		: float, allowed difference between two similar variation patterns in order to build consensus
  consthres|ct  		: float, faction of positions that have to agree to build a consensus variation position for an accession
  maxnoinfo|mi  		: float, maximum fraction of missing information allowed in variation patterns
  maxnoinfo_raw|mir		: float, maximum fraction of missing information allowed in variant positions
  threads|th    		: integer, # of threads to use
  minwgocc|mio  		: integer, minimum whole genome occurence allowed (curation step)
  maxwgocc|mao  		: integer, maximum whole genome occurence allowed (curation step)
  nbacc|nba     		: integer, # of accessions in the dataset (required when input file is not vcf)
  homozmismatches|homi		: integer, allowed # of homozygous mismatches between variation patterns
=cut

#Default variable
my ($buildblock, $cleanblock, $buildlink, $help, $verbose, $file) = (0, 0, 0, 0, 0, 'n/a');
my ($minva, $minre, $minco, $bsize, $ssize, $diff, $consthres, $maxnoinfo, $maxnoinfo_raw,  $homomis) = (2, 2, 3, 100000, 100000, 0.06, 0.90, 0.10, 0.50, 0);
my ($threads, $minwgocc, $maxwgocc, $mishp, $mioccshp, $outid, $nacc) = (1, 2, 100, 1, 1, 'unset', -1);

#Get options
GetOptions( 'buildblock|bb'   => \$buildblock,
            'cleanblock|cb'   => \$cleanblock,
            'buildlink|bl'    => \$buildlink,
            'help|h'         => \$help,
            'verbose|v'      => \$verbose,
            'file|f=s'       => \$file,
            'bsize|bs:i' => \$bsize,
            'ssize|ss:i' => \$ssize,
            'minva|mv:i' => \$minva,
            'minre|mr:i' => \$minre,
            'minco|mc:i' => \$minco,
            'diff|d:f' => \$diff,
            'consthres|ct:f' => \$consthres,
            'maxnoinfo|mi:f' => \$maxnoinfo,
            'maxnoinfo_raw|mir:f' => \$maxnoinfo_raw,
            'threads|th:i' => \$threads,
            'minwgocc|mio:i' => \$minwgocc,
            'maxwgocc|mao:i' => \$maxwgocc,
            'minsharedpat|msp:i' => \$mishp,
            'minsharedocc|mso:i' => \$mioccshp,
            'outputid|oid:s' => \$outid,
            'nbacc|nba:s' => \$nacc,
            'homozmismatches|homi:s' => \$homomis,
);

# Print help
die `pod2text $0` if $help;

#Output intermediate and definitive IDs
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $run_number   = $outid."_".$mday."_".($mon+1)."_".($year+1900)."_".$hour."_".$min."_".$sec;
my $output       = '_'.$outid.'_'.$run_number.'.bps';
my $output_human = '_'.$outid.'_'.$run_number.'.hps';
my $boutb = "block".$output;
my $bouth = "block".$output_human;
my $loutb = "links".$output;
my $louth = "links".$output_human;

#Argument hash to pass to all functions so they can pick what they need inside
my %arg;
$arg{'file'}      = $file;
$arg{'bsize'}      = $bsize;
$arg{'ssize'}      = $ssize;
$arg{'minva'}      = $minva;
$arg{'minre'}      = $minre - 1; # The reference is a reference like accession! If I allow 1N I in fact allow 2 reference-like accessions
$arg{'minco'}      = $minco;
$arg{'diff'}       = $diff;
$arg{'consthres'}  = $consthres;
$arg{'maxnoinfo'}  = $maxnoinfo;
$arg{'maxnoinfo_raw'}  = $maxnoinfo_raw;
$arg{'threads'}    = $threads;
$arg{'minwgocc'}   = $minwgocc;
$arg{'maxwgocc'}   = $maxwgocc;
$arg{'outputb'}    = $output;
$arg{'outputh'}    = $output_human;
$arg{'mishp'}      = $mishp;
$arg{'mioccshp'}   = $mioccshp;
$arg{'verbose'}   = $verbose;
$arg{'nacc'}   = $nacc;
$arg{'homi'}   = $homomis;

# Print options
print <<DOC;
Your output files have the timestamp $run_number .
You initiated the program with:
buildblock:$buildblock
cleanblock:$cleanblock
buildlink: $buildlink
verbose:   $verbose
file:      $file
bsize:     $bsize
ssize:     $ssize
minva:     $minva
minre:     $minre
minco:     $minco
diff:      $diff
homozygous mismatch penalties: $homomis
consthres: $consthres
maxnoinfo: $maxnoinfo
threads:   $threads
minwgocc:  $minwgocc
maxwgocc:  $maxwgocc
mishp      $mishp
mioccshp   $mioccshp (not used, use varcon_connection_validation.pl instead)

DOC

#GLOBALS
my @blocky_mcblock;
my @curated_blocky_mcblock;
my @blocky_mcblock_mclink;


#Building blocks

if ($buildblock){
  my $nb_accessions = &get_number_accessions(\%arg);
  print "You have $nb_accessions accessions in your dataset.\n";
  $arg{'nacc'}       = $nb_accessions;
  print("Loading data.\n");
  my @patterns = &load_data(\%arg);
  print("You have ",scalar(@patterns), " patterns in your dataset.\n");
  my @chromosome_list = uniq(grep defined, map {$_->{chrom}} @patterns);
  @blocky_mcblock = &create_blocks(\@chromosome_list, \@patterns, \%arg);
  my $blockdump = Dumper(@blocky_mcblock);
  open(my $fh, '>', $bouth);
  print $fh $blockdump;
  close $fh;
  store \@blocky_mcblock, $boutb;
}

#Cleaning blocks
if ($cleanblock){
  my $curated_blocky_mcblock_ref = &clean(\@blocky_mcblock, \%arg);
  @curated_blocky_mcblock = @$curated_blocky_mcblock_ref;
}

#Building links
if ($buildlink){
  @blocky_mcblock_mclink = &find_links(\@blocky_mcblock, \%arg) if ($cleanblock == 0);
  @blocky_mcblock_mclink = &find_links(\@curated_blocky_mcblock, \%arg) if ($cleanblock);
  my $linkdump = Dumper(@blocky_mcblock_mclink);
  open(my $fh, '>', $louth);
  print $fh $linkdump;
  close $fh;
  store \@blocky_mcblock_mclink, $loutb;
}
