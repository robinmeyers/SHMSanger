#!/usr/bin/env perl

##
## This program analyzes somatic hypermutation
## from Sanger sequencing data
##
## run with "--help" for usage information
##
## Robin Meyers <robin.meyers@childrens.harvard.edu>, 05mar2013

use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use List::MoreUtils qw(pairwise);
use Data::Dumper;
use threads;
use threads::shared;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "SHMSangerHelper.pl";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );



# Forward declarations
sub parse_command_line;
sub read_in_meta_file;
sub check_existance_of_files;
sub process_experiment ($$);
sub create_summary;


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $indir;
my $outdir;
my $refdir;

my $userscreenopt = "";
my $userphredopt = "";
my $userphrapopt = "";
my $userswopt = "";
my $max_threads = 2;
my $phred;
my $min_qual = 15;
my $min_sw_score = 100;
my $shm_threshold = 0;
my $minplotsubs = 0;
my $ow;

# Global variabless
my %meta_hash;
my $defaultphredopt = "-trim_alt \"\" -trim_cutoff 0.05 -trim_fasta";
my $defaultscreenopt = "-minmatch 12 -minscore 20";
my $defaultphrapopt = "-minscore 100";
my $defaultswopt = "-gapopen 13 -gapextend 0.1 -supper1 -supper2 -datafile /usr/local/emboss/data/EDNACUST";
my $sw_outfmt = "-aformat markx10";
my $exptfile;
my $shmexptfile;
my $delshmexptfile;
my $clonefile;
my $shmclonefile;
my $delshmclonefile;
my $bxexptfile;
my $mutfile;
#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];


my $phredopt = manage_program_options($defaultphredopt,$userphredopt);
my $screenopt = manage_program_options($defaultscreenopt,$userscreenopt);
my $phrapopt = manage_program_options($defaultphrapopt,$userphrapopt);
my $swopt = manage_program_options($defaultswopt,$userswopt) . " $sw_outfmt";


read_in_meta_file;

check_existance_of_files;



my @threads = ();


# Standard multithreading loop - each experiment is a thread
foreach my $expt_id (sort keys %meta_hash) {

    while (1) {

    # joins any threads if possible
        foreach my $thr (@threads) {
            $thr->join() if $thr->is_joinable();
        }

        my @running = threads->list(threads::running);
        
        # if there are open threads, create a new one, push it onto list, and exit while loop
        if (scalar @running < $max_threads) {
            my $thr = threads->create( sub {
                        my $t0_expt = [gettimeofday];
                        print "\nStarting $expt_id\n";
                        
                        process_experiment($expt_id, $meta_hash{$expt_id} );
                        my $t1 = tv_interval($t0_expt);
                        printf("\nFinished %s in %.2f seconds.\n", $expt_id,$t1);
                    });
            push(@threads,$thr);
            sleep(1);
            last;
        }
        sleep(1);
    } 
}

# waits for all threads to finish
while( scalar threads->list(threads::all) > 0) {
  for my $thr (@threads) {
      $thr->join() if $thr->is_joinable;
  }
  sleep(1);
}


create_summary;

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#


sub create_summary {


  # Define stats and base exchange hashes
  my %stats;
  my %bx;


  # Open file handles and print headers
  # Experiment and Clone file for each of three levels
  # 1) All clones
  # 2) Clones having at least n mutations
  # 3) Clones having at least n mutations and a deletion at least n bases long
  my $exptsfh = IO::File->new(">$exptfile");
  my $shmexptsfh = IO::File->new(">$shmexptfile");
  my $delshmexptsfh = IO::File->new(">$delshmexptfile");
  my $clonesfh = IO::File->new(">$clonefile");
  my $shmclonesfh = IO::File->new(">$shmclonefile");
  my $delshmclonesfh = IO::File->new(">$delshmclonefile");
  $exptsfh->print(join("\t",qw(Expt Allele Clones Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN))."\n");
  $shmexptsfh->print(join("\t",qw(Expt Allele Clones Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN))."\n");
  $delshmexptsfh->print(join("\t",qw(Expt Allele Clones Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN))."\n");
  $clonesfh->print(join("\t",qw(Expt Allele Clone Bp Subs Del DelBp LargeDel Ins InsBp RefA RefC RefG RefT RefN Coords))."\n");
  $shmclonesfh->print(join("\t",qw(Expt Allele Clone Bp Subs Del DelBp LargeDel Ins InsBp RefA RefC RefG RefT RefN Coords))."\n");
  $delshmclonesfh->print(join("\t",qw(Expt Allele Clone Bp Subs Del DelBp LargeDel Ins InsBp RefA RefC RefG RefT RefN Coords))."\n");

  # Base exchange file has all mutations listed in a single row to allow for tabulated experiment data
  my $bxexptsfh = IO::File->new(">$bxexptfile");
  my @bx_combs = ();
  my @bases = qw(A C G T);
  foreach my $b1 (@bases) {
    foreach my $b2 (@bases) {
      next if $b1 eq $b2;
      push(@bx_combs,$b1."->".$b2);
    }
  }
  $bxexptsfh->print(join("\t","Expt","Allele",@bx_combs)."\n");

  # For 
  foreach my $expt (sort keys %meta_hash) {

    # Open experiment clonefile
    my $clonefh = IO::File->new("<".$meta_hash{$expt}->{clonefile});
    my $csv = Text::CSV->new({sep_char => "\t"});
    my $header = $csv->getline($clonefh);
    $csv->column_names(@$header);

    # Open experiment base exchange file
    my $bxfh = IO::File->new("<".$meta_hash{$expt}->{base_ex});
    my $bx_csv = Text::CSV->new({sep_char => "\t"});
    my $bx_header = $bx_csv->getline($bxfh);
    $bx_csv->column_names(@$bx_header);

    @{$bx{$expt}}{@bx_combs} = (0) x @bx_combs;

    # Read in each clones base exchange exchange stats
    while (my $clone = $bx_csv->getline_hr($bxfh)) {
      $bx{$expt}->{$clone->{ID}} = $clone;
    }
    $bxfh->close;


    $stats{$expt}->{Clones} = 0;
    $stats{$expt}->{Bp} = 0;
    $stats{$expt}->{Subs} = 0;
    $stats{$expt}->{Del} = 0;
    $stats{$expt}->{DelBp} = 0;
    $stats{$expt}->{Ins} = 0;
    $stats{$expt}->{InsBp} = 0;
    $stats{$expt}->{RefA} = 0;
    $stats{$expt}->{RefC} = 0;
    $stats{$expt}->{RefG} = 0;
    $stats{$expt}->{RefT} = 0;
    $stats{$expt}->{RefN} = 0;

    $stats{$expt}->{SHMClones} = 0;
    $stats{$expt}->{SHMBp} = 0;
    $stats{$expt}->{SHMSubs} = 0;
    $stats{$expt}->{SHMDel} = 0;
    $stats{$expt}->{SHMDelBp} = 0;
    $stats{$expt}->{SHMIns} = 0;
    $stats{$expt}->{SHMInsBp} = 0;
    $stats{$expt}->{SHMRefA} = 0;
    $stats{$expt}->{SHMRefC} = 0;
    $stats{$expt}->{SHMRefG} = 0;
    $stats{$expt}->{SHMRefT} = 0;
    $stats{$expt}->{SHMRefN} = 0;

    $stats{$expt}->{DelSHMClones} = 0;
    $stats{$expt}->{DelSHMBp} = 0;
    $stats{$expt}->{DelSHMSubs} = 0;
    $stats{$expt}->{DelSHMDel} = 0;
    $stats{$expt}->{DelSHMDelBp} = 0;
    $stats{$expt}->{DelSHMIns} = 0;
    $stats{$expt}->{DelSHMInsBp} = 0;
    $stats{$expt}->{DelSHMRefA} = 0;
    $stats{$expt}->{DelSHMRefC} = 0;
    $stats{$expt}->{DelSHMRefG} = 0;
    $stats{$expt}->{DelSHMRefT} = 0;
    $stats{$expt}->{DelSHMRefN} = 0;


    # Read through clone file for each experiment
    while (my $clone = $csv->getline_hr($clonefh)) {

      # Clones are marked as 0 Bp if they are a repeat (or unaligned)
      next unless $clone->{Bp} > 0;
      # !!!!!!!What is this next line doing?
      # next if $clone->{LargeDel} > 10 && $clone->{Subs} < 2;
      # Add clone stats to overall experiment stats
      $stats{$expt}->{Clones}++;
      $stats{$expt}->{Bp} += $clone->{Bp};
      $stats{$expt}->{Subs} += $clone->{Subs};
      $stats{$expt}->{Del} += $clone->{Dels};
      $stats{$expt}->{DelBp} += $clone->{DelBp};
      $stats{$expt}->{Ins} += $clone->{Ins};
      $stats{$expt}->{InsBp} += $clone->{InsBp};
      $stats{$expt}->{RefA} += $clone->{RefA};
      $stats{$expt}->{RefC} += $clone->{RefC};
      $stats{$expt}->{RefG} += $clone->{RefG};
      $stats{$expt}->{RefT} += $clone->{RefT};
      $stats{$expt}->{RefN} += $clone->{RefN};


      # Experiment bx stats
      my @tmp1 = @{$bx{$expt}}{@bx_combs};
      # Clone bx stats
      my @tmp2 = @{$bx{$expt}->{$clone->{ID}}}{@bx_combs};

      our ($a,$b);

      # Update experiment bx stats by adding in clone bx stats
      @{$bx{$expt}}{@bx_combs} = pairwise { $a + $b } @tmp1, @tmp2;


      # Print clones stats to Clone summary file
      $clonesfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},
                                $clone->{ID},
                                $clone->{Bp},
                                $clone->{Subs},
                                $clone->{Dels},
                                $clone->{DelBp},
                                $clone->{LargeDel},
                                $clone->{Ins},
                                $clone->{InsBp},
                                $clone->{RefA},
                                $clone->{RefC},
                                $clone->{RefG},
                                $clone->{RefT},
                                $clone->{RefN},
                                $clone->{Coords})."\n");

      # Apply Level 2 criterion - minimum of n mutations
      next unless $clone->{Subs} > $shm_threshold;

      # Add mutated clones to experiment Level 2 summary stats
      $stats{$expt}->{SHMClones}++;
      $stats{$expt}->{SHMBp} += $clone->{Bp};
      $stats{$expt}->{SHMSubs} += $clone->{Subs};
      $stats{$expt}->{SHMDel} += $clone->{Dels};
      $stats{$expt}->{SHMDelBp} += $clone->{DelBp};
      $stats{$expt}->{SHMIns} += $clone->{Ins};
      $stats{$expt}->{SHMInsBp} += $clone->{InsBp};
      $stats{$expt}->{SHMRefA} += $clone->{RefA};
      $stats{$expt}->{SHMRefC} += $clone->{RefC};
      $stats{$expt}->{SHMRefG} += $clone->{RefG};
      $stats{$expt}->{SHMRefT} += $clone->{RefT};
      $stats{$expt}->{SHMRefN} += $clone->{RefN};
      # Print clone stats to clone Level 2 summary stats
      $shmclonesfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},
                                $clone->{ID},
                                $clone->{Bp},
                                $clone->{Subs},
                                $clone->{Dels},
                                $clone->{DelBp},
                                $clone->{LargeDel},
                                $clone->{Ins},
                                $clone->{InsBp},
                                $clone->{RefA},
                                $clone->{RefC},
                                $clone->{RefG},
                                $clone->{RefT},
                                $clone->{RefN},
                                $clone->{Coords})."\n");

      # Apply Level 3 criterion - Must contain a deletion at least n bases large
      next unless $clone->{LargeDel} > 0;
      # Add mutated/deleted clones to experiment Level 3 summary stats
      $stats{$expt}->{DelSHMClones}++;
      $stats{$expt}->{DelSHMBp} += $clone->{Bp};
      $stats{$expt}->{DelSHMSubs} += $clone->{Subs};
      $stats{$expt}->{DelSHMDel} += $clone->{Dels};
      $stats{$expt}->{DelSHMDelBp} += $clone->{DelBp};
      $stats{$expt}->{DelSHMIns} += $clone->{Ins};
      $stats{$expt}->{DelSHMInsBp} += $clone->{InsBp};
      $stats{$expt}->{DelSHMRefA} += $clone->{RefA};
      $stats{$expt}->{DelSHMRefC} += $clone->{RefC};
      $stats{$expt}->{DelSHMRefG} += $clone->{RefG};
      $stats{$expt}->{DelSHMRefT} += $clone->{RefT};
      $stats{$expt}->{DelSHMRefN} += $clone->{RefN};
      # Pring clone stats to clone Level 3 summary stats
      $delshmclonesfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},
                                $clone->{ID},
                                $clone->{Bp},
                                $clone->{Subs},
                                $clone->{Dels},
                                $clone->{DelBp},
                                $clone->{LargeDel},
                                $clone->{Ins},
                                $clone->{InsBp},
                                $clone->{RefA},
                                $clone->{RefC},
                                $clone->{RefG},
                                $clone->{RefT},
                                $clone->{RefN},
                                $clone->{Coords})."\n");

    }
    # Done reading through experiment's clone file
    $clonefh->close;


    # Calculate experiment mutation rate for Level 1 and Level 2 clones
    # Substitutions/(Total basepairs aligned to - deleted basepairs)
    my $mutrate = $stats{$expt}->{Bp} - $stats{$expt}->{DelBp} > 0 ? $stats{$expt}->{Subs}/($stats{$expt}->{Bp} - $stats{$expt}->{DelBp}) : "";
    my $shmmutrate = $stats{$expt}->{SHMBp} - $stats{$expt}->{SHMDelBp} > 0 ? $stats{$expt}->{SHMSubs}/($stats{$expt}->{SHMBp} - $stats{$expt}->{SHMDelBp}) : "";

    # Print experiment summary stats for Levels 1, 2, and 3
    $exptsfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},$stats{$expt}->{Clones},
                                  $stats{$expt}->{Bp},
                                  $stats{$expt}->{Subs},
                                  $stats{$expt}->{Del},                                  
                                  $stats{$expt}->{DelBp},
                                  $stats{$expt}->{Ins},
                                  $stats{$expt}->{InsBp},                                  
                                  $stats{$expt}->{RefA},
                                  $stats{$expt}->{RefC},
                                  $stats{$expt}->{RefG},
                                  $stats{$expt}->{RefT},
                                  $stats{$expt}->{RefN})."\n");

    $shmexptsfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},$stats{$expt}->{SHMClones},
                                  $stats{$expt}->{SHMBp},
                                  $stats{$expt}->{SHMSubs},
                                  $stats{$expt}->{SHMDel},                                  
                                  $stats{$expt}->{SHMDelBp},
                                  $stats{$expt}->{SHMIns},
                                  $stats{$expt}->{SHMInsBp},                                  
                                  $stats{$expt}->{SHMRefA},
                                  $stats{$expt}->{SHMRefC},
                                  $stats{$expt}->{SHMRefG},
                                  $stats{$expt}->{SHMRefT},
                                  $stats{$expt}->{SHMRefN})."\n");

    $delshmexptsfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},$stats{$expt}->{DelSHMClones},
                                  $stats{$expt}->{DelSHMBp},
                                  $stats{$expt}->{DelSHMSubs},
                                  $stats{$expt}->{DelSHMDel},                                  
                                  $stats{$expt}->{DelSHMDelBp},
                                  $stats{$expt}->{DelSHMIns},
                                  $stats{$expt}->{DelSHMInsBp},                                  
                                  $stats{$expt}->{DelSHMRefA},
                                  $stats{$expt}->{DelSHMRefC},
                                  $stats{$expt}->{DelSHMRefG},
                                  $stats{$expt}->{DelSHMRefT},
                                  $stats{$expt}->{DelSHMRefN})."\n");
    # Print base exchange stats
    $bxexptsfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},@{$bx{$expt}}{@bx_combs})."\n");

    
  }

  $exptsfh->close;
  $shmexptsfh->close;
  $delshmexptsfh->close;
  $clonesfh->close;
  $shmclonesfh->close;
  $delshmclonesfh->close;
  $bxexptsfh->close;

  # All mutations included in mutations file
  System("cat $outdir/experiments/*/*_muts.txt > $mutfile");

  mkdir "$outdir/group_viz";
  # Create mutation visualizations by groups
  my $viz_cmd = "Rscript $FindBin::Bin/../R/mutationVizGrouped.R $meta_file $mutfile $clonefile $refdir $outdir/group_viz/ minsubs=$minplotsubs";
  System("$viz_cmd > $outdir/group_viz/R.out 2>&1");


  # Create grouped stats for experiments for Levels 1, 2, and 3
  my $group_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $exptfile $clonefile $meta_file $outdir/Group.txt grouping=\"genotype,allele,tissue,pna\"";
  System("$group_cmd >> $outdir/group_viz/R.out 2>&1");
  my $shm_group_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $shmexptfile $shmclonefile $meta_file $outdir/GroupSHM.txt grouping=\"genotype,allele,tissue,pna\"";
  System("$shm_group_cmd >> $outdir/group_viz/R.out 2>&1");
  my $del_shm_group_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $delshmexptfile $delshmclonefile $meta_file $outdir/GroupSHMDel.txt grouping=\"genotype,allele,tissue,pna\"";
  System("$del_shm_group_cmd >> $outdir/group_viz/R.out 2>&1");

  # Group base exchange stats
  my $bx_group_cmd = "Rscript $FindBin::Bin/../R/baseExGroup.R $bxexptfile $meta_file $outdir/GroupBXTabular.txt";
  System("$bx_group_cmd");

  # Un-tabulate GroupBXTabular file
  my $infh = IO::File->new("<$outdir/GroupBXTabular.txt");
  my $outfh = IO::File->new(">$outdir/GroupBX.txt");

  my $csv = Text::CSV->new({sep_char => "\t"});
  my $header = $csv->getline($infh);
  $csv->column_names(@$header);


  while (my $group = $csv->getline_hr($infh)) {
    $outfh->print(join(" ",$group->{genotype},$group->{allele},$group->{tissue},$group->{pna})."\n");
    $outfh->print(join("\t","Base",@bases)."\n");
    foreach my $b1 (@bases) {
      my @row = ($b1);
      foreach my $b2 (@bases) {
        if ($b1 eq $b2) {
          push(@row, "-");
        } else {
          push(@row,$group->{$b1."..".$b2})
        }
      }
      $outfh->print(join("\t",@row)."\n");
    }
    $outfh->print("\n");

  }

  # Same as grouping stats but aggregate on mouse number as well
  my $mouse_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $exptfile $clonefile $meta_file $outdir/Mouse.txt grouping=\"genotype,allele,mouse,tissue,pna\"";
  System("$mouse_cmd >> $outdir/group_viz/R.out 2>&1");
  my $shm_mouse_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $shmexptfile $shmclonefile $meta_file $outdir/MouseSHM.txt grouping=\"genotype,allele,mouse,tissue,pna\"";
  System("$shm_mouse_cmd >> $outdir/group_viz/R.out 2>&1");
  my $del_shm_mouse_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $delshmexptfile $delshmclonefile $meta_file $outdir/MouseSHMDel.txt grouping=\"genotype,allele,mouse,tissue,pna\"";
  System("$del_shm_mouse_cmd >> $outdir/group_viz/R.out 2>&1");
}



sub process_experiment ($$) {

	my $expt_id = shift;
	my $expt_hash = shift;

  # 'raw' key will have been set if fasta and qual files are not found
  if ($phred || defined $expt_hash->{raw}) {
    # Apply phred on the experiment trace files for each clone to create fasta and qual files
    # Consult phred documentation for further info
    my $phred_cmd = join(" ","phred -id",$expt_hash->{raw},"-sa",$expt_hash->{phred},"-qa",$expt_hash->{phred}.".qual",$phredopt,">>",$expt_hash->{exptdir}."/phred.out 2>&1");
    System($phred_cmd);
  }

  # Remove any vector sequence from clone fasta files 
  my $screen_cmd = join(" ","cross_match",$expt_hash->{phred},$expt_hash->{vector},$screenopt,"-screen",">>",$expt_hash->{exptdir}."/screen.out 2>&1");
  System($screen_cmd);

  $expt_hash->{seqdir} = $expt_hash->{exptdir}."/sequences";

  mkdir $expt_hash->{seqdir} unless -d $expt_hash->{seqdir};

  # Move phredded and screened (fasta/qual) files to results directory
  $expt_hash->{phred} =~ s/\/\//\//g;
  my ($phr_path,$phr_name,$phr_ext ) = parseFilename($expt_hash->{phred});
  System(join(" ","mv",$expt_hash->{phred}.".screen",$expt_hash->{seqdir}."/".$phr_name.$phr_ext),1);
  $expt_hash->{phred} = $expt_hash->{seqdir}."/".$phr_name.$phr_ext;

  $expt_hash->{qual} =~ s/\/\//\//g;
  my ($qual_path,$qual_name,$qual_ext) = parseFilename($expt_hash->{qual});
  System(join(" ","cp",$expt_hash->{qual},$expt_hash->{seqdir}."/".$qual_name.$qual_ext),1);
  $expt_hash->{qual} = $expt_hash->{seqdir}."/".$qual_name.$qual_ext;


  $expt_hash->{phrapdir} = $expt_hash->{exptdir}."/phrap";
  mkdir $expt_hash->{phrapdir} unless -d $expt_hash->{phrapdir};

  # Use phrap to assembly the sequence pairs for each clone
  phrap_sequence_pairs($expt_id,$expt_hash,$phrapopt);

  $expt_hash->{water} = $expt_hash->{exptdir} . "/${expt_id}.water";
  $expt_hash->{waternice} = $expt_hash->{exptdir} . "/${expt_id}_nice.water";

  # Run the smith-waterman program from EMBOSS to align the contig to the reference
  (my $niceswopt = $swopt) =~ s/-aformat markx10//;
  System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{contig},"-outfile",$expt_hash->{water},$swopt));
  # Run again with a pretty output format
  System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{contig},"-outfile",$expt_hash->{waternice},$niceswopt));

  # Tabulate mutations from the alignment
  parse_smith_water_file($expt_id,$expt_hash,$min_sw_score,$min_qual);

  # Remove duplicate clones
  System(join(" ","Rscript $FindBin::Bin/../R/removeDupClones.R",$expt_hash->{mutfile},$expt_hash->{clonefile}));

  my $vizdir = $outdir."/viz";
  mkdir $vizdir;
  # Create mutation vizualisations
  System(join(" ","Rscript $FindBin::Bin/../R/mutationVizSuite.R",$expt_hash->{mutfile},$expt_hash->{clonefile},$expt_hash->{reference},"$vizdir/$expt_id","tstart=".$expt_hash->{start},"tend=".$expt_hash->{end},"minsubs=$minplotsubs",">",$expt_hash->{exptdir}."/R.out 2>&1"));


}


sub check_existance_of_files {
	print "\nSearching for sequence files...\n";


	foreach my $expt_id (sort keys %meta_hash) {
		my $input = $indir."/".$expt_id;
		my @exts = qw(.fa .fasta);
		foreach my $ext (@exts) {
      # Set fasta and qual files if they can be found
			if (-r $input.$ext && -r $input.$ext.".qual") {
				$meta_hash{$expt_id}->{phred} = $input.$ext;
        $meta_hash{$expt_id}->{qual} = $input.$ext.".qual";
				last;
			}
		}
    # If phred-ed fasta file not found (or phred is explicitly set in options)
    if (! defined $meta_hash{$expt_id}->{phred} || $phred) {
      if (-d $input && scalar listFilesInDir($input) > 0 ) {
        # Set input trace file directory in experiment hash
        $meta_hash{$expt_id}->{raw} = $input;
        $meta_hash{$expt_id}->{phred} = $input.".fa";
        $meta_hash{$expt_id}->{qual} = $meta_hash{$expt_id}->{phred}.".qual";

      } else {
        croak "Error: No valid experiment $expt_id in $indir" ;
      }
    }
	}
	

	print "\nSearching for reference files...\n";
  # Make sure the reference files exist in the indicated directory
	foreach my $expt_id (sort keys %meta_hash) {
		my $reffile = $refdir."/".$meta_hash{$expt_id}->{reference};
		croak "Error: Could not locate reference file $reffile in $indir" unless (-r $reffile);
		$meta_hash{$expt_id}->{reference} = $reffile;
    my $vecfile = $refdir."/".$meta_hash{$expt_id}->{vector};
    croak "Error: Could not locate reference file $vecfile in $indir" unless (-r $vecfile);
    $meta_hash{$expt_id}->{vector} = $vecfile;
	}
	
	print "Done.\n";
  mkdir "$outdir/experiments";
  mkdir "$outdir/mergenice";
  # Define and check on experiment output files
	print "\nChecking on output files...\n";
	foreach my $expt_id (sort keys %meta_hash) {
    my $exptdir = "$outdir/experiments/$expt_id/";
    $meta_hash{$expt_id}->{exptdir} = $exptdir;
    mkdir $meta_hash{$expt_id}->{exptdir};
		my $file = $exptdir."/".$expt_id."_muts.txt";
    my $stats = $exptdir."/".$expt_id."_clones.txt";
    my $base_ex = $exptdir."/".$expt_id."_baseEx.txt";
    my $mergenice = $outdir."/mergenice/".$expt_id."_mergenice.txt";
    # Only actually check mutation file
		croak "Error: Output file $file already exists and overwrite switch --ow not set" if (-r $file && ! defined $ow);
		
    # Try to create file and croak unless it's writable
    System("touch $file",1);
		croak "Error: Cannot write to $file" unless (-w $file);

    # Set output files in the experiment hash
		$meta_hash{$expt_id}->{mutfile} = $file;
    $meta_hash{$expt_id}->{clonefile} = $stats;
    $meta_hash{$expt_id}->{base_ex} = $base_ex;
    $meta_hash{$expt_id}->{mergenice} = $mergenice;

	}

  # Define summary files
	$exptfile = "$outdir/Expts.txt";
  $shmexptfile = "$outdir/ExptsSHM.txt";
  $delshmexptfile = "$outdir/ExptsSHMDel.txt";
  $clonefile = "$outdir/Clones.txt";
  $shmclonefile = "$outdir/ClonesSHM.txt";
  $delshmclonefile = "$outdir/ClonesSHMDel.txt";
  $bxexptfile = "$outdir/ExptsBX.txt";

  # Check if Mutations summary file can be created and is writable
  $mutfile = "$outdir/Mutations.txt";
	croak "Error: Output file $exptfile already exists and overwrite switch --ow not set" if (-r $exptfile && ! defined $ow);
	System("touch $exptfile",1);
	croak "Error: Cannot write to $exptfile" unless (-w $exptfile);
	print "Done.\n";

}


sub read_in_meta_file {
	System("perl -pi -e 's/\\r/\\n/g' $meta_file");

	print "\nReading in meta file...\n";

	my $meta = IO::File->new("<$meta_file");
	my $csv = Text::CSV->new({sep_char => "\t"});
	my $header = $csv->getline($meta);
	$csv->column_names(@$header);

	while (my $row = $csv->getline_hr($meta)) {
    next unless $row->{experiment} =~ /\S/;
		my $expt_id = $row->{experiment};
		$meta_hash{$expt_id} = $row;
	}

}


sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( 
														"meta=s" => \$meta_file ,
														"in=s" => \$indir ,
														"out=s" => \$outdir ,
                            "ref=s" => \$refdir ,
                            "screenopt=s" => \$userscreenopt ,
														"phredopt=s" => \$userphredopt ,
                            "phrapopt=s" => \$userphrapopt ,
                            "swopt=s" => \$userswopt ,
                            "minscore=i" => \$min_sw_score ,                            
                            "minqual=i" => \$min_qual ,
														"threads=i" => \$max_threads ,
                            "minplotsubs=i" => \$minplotsubs,
                            "phred" => \$phred ,
														"ow" => \$ow ,
														"help" => \$help

				            			);
	
	usage() if ($help);


  #Check options

  croak "Error: must specifiy --meta" unless (defined $meta_file);
  croak "Error: cannot find $meta_file" unless (-r $meta_file);
  croak "Error: must specify --in" unless (defined $indir);
  croak "Error: cannot find $indir" unless (-d $indir);
  croak "Error: must specify --out" unless (defined $outdir);
	unless (-d $outdir) {
  	mkdir $outdir or croak "Error: output directory $outdir does not exist and cannot be created";
  }
  croak "Error: changing -aformat is not allowed" if ($userswopt =~ /aformat/);
	exit unless $result;
}


sub usage()
{
print<<EOF;
SHMSanger.pl, by Robin Meyers, 05mar2013

This program analyzes somatic hypermutation from Sanger sequencing data.

Usage: $0 --metafile FILE (--in FILE | --indir DIR) --outdir DIR
		[--bcmismatch N] [--blastopt "-opt val"] [--threads N] [--ow]

Arguments (defaults in parentheses):

$arg{"--meta","Tab-delimited file containing experiment information"}
$arg{"--in","Input directory"}
$arg{"--out","Output directory"}
$arg{"--ref","Reference directory"}
$arg{"--screenopt","Specify cross_match options for vector screening",$defaultscreenopt}
$arg{"--threads","Number of processing core threads to run on",$max_threads}
$arg{"--phred","Force phred to run"}
$arg{"--minscore","Minimum score to accept alignment",$min_sw_score}
$arg{"--minqual","Minimum quality to accept mutation",$min_qual}
$arg{"--minplotsubs","Minimum number of substitutions a clone must have to be included in the R plots",$minplotsubs}
$arg{"--ow","Overwrite output files if they already exist"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
