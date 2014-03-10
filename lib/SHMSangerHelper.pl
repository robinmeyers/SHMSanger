#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use File::Basename;
use IO::File;
use Text::CSV;
use List::Util qw(min max sum);
use List::MoreUtils qw(pairwise);


sub argument {
	my $var = shift;
	my $description = shift;
	my $default = shift;

	return sprintf("  \%\-16s - %s\n",
		(defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
		$description);
}

sub listFilesInDir ($)
{
	my $path = shift;
	croak "Error: path $path not found" unless -d $path;
	my ($record, @data, @filenames);
	my @list = `ls -l $path`;

	foreach $record (@list) {
		next if $record =~ /^d/;
		@data = split(/\s+/,$record);
		next unless defined($data[8]);
		push(@filenames, join("/",$path,$data[8]));
	}
	
	carp "Warning: nothing found in path $path" unless @filenames > 0;

	return @filenames;
}

sub parseFilename ($) {
	my $fullname = shift;
	my ($name, $path, $ext) = fileparse($fullname, qr/\.\w{2,5}$/);
	return ($path, $name, $ext);
}

#invert exit status of system call
sub System ($;$) {
	my $cmd = shift;
	my $quiet = shift;
	$quiet = 0 unless defined $quiet;
	print "$cmd\n" unless $quiet;
	my $status = system($cmd);
	return !$status;
}

sub mean {
    return sum(@_)/@_;
}

sub read_fasta ($) {

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = $fh->getline();
	return () unless defined $seq_name;
	chomp($seq_name);
	croak "Error: unexpected format - missing '\>' symbol in sequence header $seq_name" unless $seq_name =~ s/^\>\s*(\S.*)$/$1/;
	my $seq_bases = "";
	while (my $line = $fh->getline()) {
		if ($line =~ /^>/) {
			seek($fh,-length($line),1);
			last;
		}
		chomp($line);
		croak "Error: unexpected base in sequence $line" unless $seq_bases =~ /^[AGCTagctNnXx\-\.]*$/;
		$seq_bases .= $line;
	}
	
	return ($seq_name,$seq_bases);
}

sub write_fasta ($$$) {
	croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 3;

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = shift;
	my $seq_bases = shift;
 

	$seq_name = ">" . $seq_name unless $seq_name =~ /^\>.*$/;
	croak "Error: sequence contains unexpected base" unless $seq_bases =~ /^[AGCTagctNnXx\-\.]+$/;

	print $fh $seq_name,"\n";
	print $fh $seq_bases,"\n";


}

sub read_qual ($) {

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = $fh->getline();
	return () unless defined $seq_name;
	chomp($seq_name);
	croak "Error: unexpected format - missing '\>' symbol in sequence header $seq_name" unless $seq_name =~ s/^\>\s*(\S.*)$/$1/;
	my @seq_quals = ();
	while (my $line = $fh->getline()) {
		if ($line =~ /^>/) {
			seek($fh,-length($line),1);
			last;
		}
		chomp($line);
		croak "Error: unexpected character in quality $line" unless $line =~ /^[\d\s]*$/;
		push(@seq_quals,split(/\s+/,$line));
	}
	
	return ($seq_name,\@seq_quals);
}

sub write_qual ($$$;$) {
	croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 3;

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = shift;
	my $quals = shift;
	my @seq_quals = @$quals;
	my $line_size = shift;
	$line_size = 20 unless defined $line_size;

	$seq_name = ">" . $seq_name unless $seq_name =~ /^\>.*$/;

	print $fh $seq_name,"\n";

	while (scalar @seq_quals > 0) {
		my $line = join(" ",splice(@seq_quals,0,$line_size));
		croak "Error: qualities contains unexpected character" unless $line =~ /^[\d\s]*$/;
		print $fh $line,"\n";
	}

}

sub reverseComplement ($) {
	my $seq = shift;
	(my $rc = reverse($seq)) =~ tr/ACGTacgtNn/TGCAtgcaNn/;
	return $rc;
}




sub manage_program_options ($$) {

	my $defaultoptstr = shift;
	my $useroptstr = shift;

	return $defaultoptstr unless $useroptstr =~ /\S/;
	return $useroptstr unless $defaultoptstr =~ /\S/;

	my @defaultopt = split(/\s+-/,$defaultoptstr);
	my @useropt = split(/\s+-/,$useroptstr);

	my %opt_hash;
	my $optkey;
	my $optval;

	my @return_opt = ();

	my $i = 1;
	while ($i < scalar @defaultopt) {
		if ($defaultopt[$i] =~ /^\d/) {
			$defaultopt[$i-1] .= " -".$defaultopt[$i];
			splice(@defaultopt,$i,1);
			next;
		}
		$i++;
	}
	$i = 1;
	while ($i < scalar @useropt) {
		if (defined $useropt[$i] && $useropt[$i] =~ /^\d/) {
			$useropt[$i-1] .= " -".$useropt[$i];
			splice(@useropt,$i,1);
			next;
		}
		$i++;
	}

	foreach (@defaultopt) {
		if ( /^-?(\S+)\s+(\S+)$/ ) {
			$opt_hash{$1} = $2;
		} elsif ( /^-?(\S+)$/ ) {
			$opt_hash{$1} = "";
		} else {
			croak "Error: incorrect default program option string $defaultoptstr";
		}
	}

	foreach (@useropt) {
		if ( /^-?(\S+)\s+(\S+)$/ ) {
			$opt_hash{$1} = $2;
			delete $opt_hash{$1} if $2 eq "off";
		} elsif ( /^-?(\S+)$/ ) {
			$opt_hash{$1} = "";
		} else {
			croak "Error: incorrect program option string $useroptstr";
		}
	}

	foreach $optkey (sort keys %opt_hash) {
		$optval = $opt_hash{$optkey};
		if ( $optval eq "" ) {
			push(@return_opt,"-$optkey");
		} else {
			push(@return_opt,"-$optkey $optval");
		}
	}

	return(join(" ",@return_opt));


}


sub smith_water_to_reference ($$$) {
	my $expt_id = shift;
	my $expt_hash = shift;
	my $swopt = shift;

	$expt_hash->{water} = $expt_hash->{exptdir} . "/${expt_id}.water";
	$expt_hash->{waternice} = $expt_hash->{exptdir} . "/${expt_id}_nice.water";

	(my $niceswopt = $swopt) =~ s/-aformat markx10//;
	System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{fa},"-outfile",$expt_hash->{water},$swopt));
	System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{fa},"-outfile",$expt_hash->{waternice},$niceswopt));

}

sub parse_smith_water_file ($$$$) {
	my $expt_id = shift;
	my $expt_hash = shift;
	my $min_score = shift;
	my $min_qual = shift;

	print "\n$expt_id: Parsing smith-waterman alignments...\n";

	my $swfh = IO::File->new("<".$expt_hash->{water}) or croak "Error: could not open water file for parsing";
	my $outfh = IO::File->new(">".$expt_hash->{mutfile}) or croak "Error: could not write to output";
	my $statfh = IO::File->new(">".$expt_hash->{clonefile}) or croak "Error: could not write to stats file";
	my $bxfh = IO::File->new(">".$expt_hash->{base_ex}) or croak "Error: could not write to base exchange file";
	my $mnfh = IO::File->new(">".$expt_hash->{mergenice}) or croak "Error: could not write to merge nice file";


 	# Define hash
 	my %q = %{$expt_hash->{q}};

	# Read in smith-water file
	while ( my $line = $swfh->getline ) {
		chomp($line);
		last if $line =~ /^#-+$/;

		my $score;
		my $ident;
		my $overlap;
		my $refid;
		my $tstart;
		my $tend;
		my $tseq = "";
		my $clone;
		my $qstart;
		my $qend;
		my $qseq = "";

		if ($line =~ /^\>\>#\d+/) {


			until ($line =~ /^\>(\S+)\s\.\./) {
				last unless defined $line;
				$score = $1 if $line =~ /^; sw_score: (\d+)/ ;
				$ident = $1 if $line =~ /^; sw_ident: ([\d\.]+)/;
				$overlap = $1 if $line =~ /^; sw_overlap: (\d+)/;


				$line = $swfh->getline;
				chomp($line);
			}

			$refid = $1 if $line =~ /^\>(\S+)\s\.\./;

			$line = $swfh->getline;
			chomp($line);

			until ($line =~ /^\>(\S+)\s\.\./) {
				last unless defined $line;
				$tstart = $1 if $line =~ /^; al_start: (\d+)/;
				$tend = $1 if $line =~ /^; al_stop: (\d+)/;
				$tseq .= $1 if $line =~ /^([ACGTacgtNnXx-]+)/;

				$line = $swfh->getline;
				chomp($line);
			}

			$clone = $1 if $line =~ /^\>(\S+)\s\.\./;

			$line = $swfh->getline;
			chomp($line);

			until ($line =~ /^#/) {
				last unless defined $line;
				$qstart = $1 if $line =~ /^; al_start: (\d+)/;
				$qend = $1 if $line =~ /^; al_stop: (\d+)/;
				$qseq .= $1 if $line =~ /^([ACGTacgtNnXx-]+)/;

				$line = $swfh->getline;
				chomp($line);
			}

			#print join(" ",$clone,$score,$qstart,$qend,$tstart,$tend)."\n";
			if ($score >= $min_score) {
				unless ($clone =~ /\.[fr]$/) {

					$q{$clone}->{score} = $score;
					$q{$clone}->{ident} = $ident;
					$q{$clone}->{overlap} = $overlap;

					$q{$clone}->{refid} = $refid;
					$q{$clone}->{tstart} = $tstart;
					$q{$clone}->{tend} = $tend;
					$q{$clone}->{tseq} = $tseq;

					$q{$clone}->{qstart} = $qstart;
					$q{$clone}->{qend} = $qend;
					$q{$clone}->{qseq} = $qseq;

				} else {
					($clone, my $read) = $clone =~ /^(.*)\.([fr])$/;
					my $ori = $read eq 'f' ? 'for' : 'rev' ;
					$q{$clone}->{$ori}->{score} = $score;
					$q{$clone}->{$ori}->{ident} = $ident;
					$q{$clone}->{$ori}->{overlap} = $overlap;

					$q{$clone}->{$ori}->{refid} = $refid;
					$q{$clone}->{$ori}->{tstart} = $tstart;
					$q{$clone}->{$ori}->{tend} = $tend;
					$q{$clone}->{$ori}->{tseq} = $tseq;

					$q{$clone}->{$ori}->{qstart} = $qstart;
					$q{$clone}->{$ori}->{qend} = $qend;
					$q{$clone}->{$ori}->{qseq} = $qseq;
				}
			}
		}

	}
	$swfh->close;

	my @bases = qw(A C G T N);


	my $reffh = IO::File->new("<".$expt_hash->{reference}) or croak "Error: could not open reference file";
	my ($refid,$refseq) = read_fasta($reffh);
	$reffh->close;

	$statfh->print(join("\t",qw(ID Bp Subs Dels LargeDel DelBp Ins InsBp RefA RefC RefG RefT RefN Coords Notes))."\n");

	my @bx_header = ("ID");
	foreach my $b1 (@bases) {
		next if $b1 eq "N";
		foreach my $b2 (@bases) {
			next if $b2 eq "N";
			next if $b1 eq $b2;
			push(@bx_header,$b1."->".$b2);
		}
	}
	$bxfh->print(join("\t",@bx_header)."\n");

	foreach my $clone (sort keys %q) {

		# Parse alignment for each clone
		parse_sw_alignment($outfh,$expt_id,$q{$clone},$min_qual,$expt_hash->{start},$expt_hash->{end},$mnfh);
		

		if (defined $q{$clone}->{bps}) {
			$statfh->print(join("\t",$clone,$q{$clone}->{bps},
																		$q{$clone}->{sub},
																		$q{$clone}->{del},
																		$q{$clone}->{largedel},
																		$q{$clone}->{delbp},
																		$q{$clone}->{ins},
																		$q{$clone}->{insbp},
																		$q{$clone}->{$bases[0]},
																		$q{$clone}->{$bases[1]},
																		$q{$clone}->{$bases[2]},
																		$q{$clone}->{$bases[3]},
																		$q{$clone}->{$bases[4]},
																		$q{$clone}->{coords})."\n");

			# $bxfh->print(join("\t",$clone,@bases)."\n");
			my @bx_row = ($clone);
			foreach my $b1 (@bases) {
				next if $b1 eq "N";
				foreach my $b2 (@bases) {
					next if $b2 eq "N";
					next if $b1 eq $b2;
					push(@bx_row,$q{$clone}->{bx}->{$b1}->{$b2});
				}
			}
			$bxfh->print(join("\t",@bx_row)."\n");
			# 	$bxfh->print(join("\t",$base,$q{$clone}->{bx}->{$base}->{$bases[0]},
			# 															 $q{$clone}->{bx}->{$base}->{$bases[1]},
			# 															 $q{$clone}->{bx}->{$base}->{$bases[2]},
			# 															 $q{$clone}->{bx}->{$base}->{$bases[3]})."\n");
			
			# $bxfh->print(join("\n"));

		} else {
			$statfh->print(join("\t", $clone, 0,"","","","","","","","","","","","","bad alignment")."\n");
			$bxfh->print("$clone\n")
		}
	}

	$statfh->close;
	$bxfh->close;

}


sub merge_alignments ($$) {

	my $q = shift;
	my $min_qual = shift;

	#print "\nmerging alignments for " . $q->{clone} . "\n";


	return if defined $q->{score};
	my $ori;
	if (!defined $q->{for}->{score}) {
		$ori = "rev";
	}
	if (!defined $q->{rev}->{score}) {
		$ori = "for";
	}
	if (defined $ori) {
		$q->{score} = $q->{$ori}->{score};
		$q->{tseq} = $q->{$ori}->{tseq};
		$q->{qseq} = $q->{$ori}->{qseq};
		$q->{contigqual} = $q->{$ori}->{quals};
		$q->{tstart} = $q->{$ori}->{tstart};
		$q->{tend} = $q->{$ori}->{tend};
		$q->{qstart} = $q->{$ori}->{qstart};
		$q->{qend} = $q->{$ori}->{qend};

		# print "only using $ori orientation\n";
		# print "tseq: ".$q->{tstart}."-".$q->{tend}."\n";
		# print "qseq: ".$q->{qstart}."-".$q->{qend}."\n";
	} else {
		
		# print "using both alignments\n";
		# print "forward: tseq ".$q->{for}->{tstart}."-".$q->{for}->{tend}." qseq ".$q->{for}->{qstart}."-".$q->{for}->{qend}."\n";
		# print "tseq: ".$q->{for}->{tseq}."\n";
		# print "qseq: ".$q->{for}->{qseq}."\n";
		# print "reverse: tseq ".$q->{rev}->{tstart}."-".$q->{rev}->{tend}." qseq ".$q->{rev}->{qstart}."-".$q->{rev}->{qend}."\n";
		# print "tseq: ".$q->{rev}->{tseq}."\n";
		# print "qseq: ".$q->{rev}->{qseq}."\n";

		$q->{score} = $q->{for}->{score} + $q->{rev}->{score};
		$q->{tstart} = min($q->{for}->{tstart},$q->{rev}->{tstart});
		$q->{tend} = max($q->{for}->{tend},$q->{rev}->{tend});
		

		my $t = $q->{tstart};

		my $f_i = -1;
		my $r_i = -1;
		my $f_qual_i = $q->{for}->{qstart};
		my $r_qual_i = $q->{rev}->{qstart};

		my $qseq = "";
		my $tseq = "";
		my @quals = ();

		my $qpos_f;
		my $qpos_r;

		# Loop until tpos exceeds tend
		until ($t > $q->{tend}) {

			$f_i = 0 if $t == $q->{for}->{tstart};
			$r_i = 0 if $t == $q->{rev}->{tstart};



			my $tbase_f;
			my $tbase_r;
			my $qbase_f;
			my $qbase_r;
			my $qual_f;
			my $qual_r;

			unless ($f_i < 0 || $f_i >= length($q->{for}->{tseq})) {
				$tbase_f = substr($q->{for}->{tseq},$f_i,1);
				$qbase_f = substr($q->{for}->{qseq},$f_i,1);
				$qual_f = $q->{for}->{quals}->[$f_qual_i-1];
			}
			unless ($r_i < 0 || $r_i >= length($q->{rev}->{tseq})) {
				$tbase_r = substr($q->{rev}->{tseq},$r_i,1);
				$qbase_r = substr($q->{rev}->{qseq},$r_i,1);
				$qual_r = $q->{rev}->{quals}->[$r_qual_i-1];
			}

			my @t_bases = ();
			push(@t_bases, $tbase_f) if defined $tbase_f;
			push(@t_bases, $tbase_r) if defined $tbase_r;

			my @q_bases = ();
			push(@q_bases, $qbase_f) if defined $qbase_f;
			push(@q_bases, $qbase_r) if defined $qbase_r;

			# Add Ns and 0 to quality if no bases defined
			if (@t_bases == 0) {
				$tseq .= "N";
				$qseq .= "N";
				push(@quals, 0);
				# and push ahead
				$t++;
				next;
			}

			# If there is an insertion
			if (@t_bases ~~ /-/) {
				# Use forward read if reverse read not defined
				if (! defined $tbase_r) {
					$tseq .= "-";
					$qseq .= $qbase_f;
					push(@quals,$qual_f);
					$f_i++;
					$f_qual_i++;
					next;
				}
				# Use reverse read if forward read not defined
				if (! defined $tbase_f) {
					$tseq .= "-";
					$qseq .= $qbase_r;
					push(@quals,$qual_r);
					$r_i++;
					$r_qual_i++;
					next;
				}

				# Otherwise only include insertion if they both reads have it
				if ($tbase_f eq "-" && $tbase_r eq "-") {
					$tseq .= "-";
					# Include the inserted base if they match on both strands, otherwise use N
					$qseq .= $qbase_f eq $qbase_r ? $qbase_f : "N";
					push(@quals, $qbase_f eq $qbase_r ? $qual_f + $qual_r : 0);
					$f_i++;
					$r_i++;
					$f_qual_i++;
					$r_qual_i++;
					next;

				} else {
					# If only one read has it, do not include
					# Push ahead on the read that had the insertion
					$f_i++ if $tbase_f eq "-";
					$f_qual_i++ if $tbase_f eq "-";
					$r_i++ if $tbase_r eq "-";
					$r_qual_i++ if $tbase_r eq "-";
					next;
				}

			}

			# If there is a deletion
			if (@q_bases ~~ /-/) {
				# Use forward read if reverse read not defined
				if (! defined $qbase_r) {
					$tseq .= $tbase_f;
					$qseq .= "-";
					$t++;
					$f_i++;
					next;
				}
				# Use reverse read if forward read not defined
				if (! defined $qbase_f) {
					$tseq .= $tbase_r;
					$qseq .= "-";
					$t++;
					$r_i++;
					next
				}

				#Otherwise only include deletion if both reads have it
				if ($qbase_f eq "-" && $qbase_r eq "-") {
					$tseq .= $tbase_f;
					$qseq .= "-";
					$t++;
					$f_i++;
					$r_i++;
					next;
				} else {
					$tseq .= $tbase_f;
					$qseq .= "N";
					push(@quals,0);
					$t++;
					$f_i++;
					$r_i++;
					$f_qual_i++ unless $qbase_f eq "-";
					$r_qual_i++ unless $qbase_r eq "-";
					next;
				}
			}

			if (! defined $qbase_r) {
				$tseq .= $tbase_f;
				$qseq .= $qbase_f;
				push(@quals,$qual_f);
				$t++;
				$f_i++;
				$f_qual_i++;
				next;
			}

			if (! defined $qbase_f){
				$tseq .= $tbase_r;
				$qseq .= $qbase_r;
				push(@quals,$qual_r);
				$t++;
				$r_i++;
				$r_qual_i++;
				next;
			}

			if ($qbase_f eq $qbase_r) {
				$tseq .= $tbase_f;
				$qseq .= $qbase_f;
				push(@quals,$qual_f + $qual_r);				
			} else {
				$tseq .= $tbase_f;
				$qseq .= "N";
				push(@quals,0);
			}
			$t++;
			$f_i++;
			$r_i++;
			$f_qual_i++;
			$r_qual_i++;

		}


		$q->{qstart} = 1;
		$q->{qend} = scalar @quals;

		$q->{tseq} = $tseq;
		$q->{qseq} = $qseq;

		@quals = map {$_ > 80 ? 80 : $_ } @quals;

		$q->{contigqual} = \@quals;

		# print "finished with tseq ".$q->{tstart}."-".$q->{tend}." qseq ".$q->{qstart}."-".$q->{qend}."\n";
		# print "tseq: ".$q->{tseq}."\n";
		# print "qseq: ".$q->{qseq}."\n";
		# print "quals: ".join(" ",@quals)."\n";
	}

		
}


sub pretty_print_merge ($$) {

	my $mnfh = shift;
	my $q = shift;

	return unless defined $q->{for}->{qseq} && defined $q->{rev}->{qseq};

	$mnfh->print("\n\n".$q->{clone}."\n");
	$mnfh->print($q->{tseq}."\n");
	$mnfh->print($q->{qseq}."\n");
	$mnfh->print(join("",map {chr($_ + 33)} @{$q->{contigqual}})."\n");
	if (defined $q->{for}->{qseq}) {
		$mnfh->print($q->{for}->{qseq}."\n");
		$mnfh->print(join("",map {chr($_+33)} @{$q->{for}->{quals}})."\n");
	} else {
		$mnfh->print("no for\n");
	}
	if (defined $q->{rev}->{qseq}) {
		$mnfh->print($q->{rev}->{qseq}."\n");
		$mnfh->print(join("",map {chr($_+33)} @{$q->{rev}->{quals}})."\n");
	} else {
		$mnfh->print("no rev\n");
	}



	my $t = $q->{tstart};

	my $i = 0;
	my $line_width = 60;
	my $t_seq = "";
	my $q_seq = "";
	my $f_seq = "";
	my $r_seq = "";
	my $q_qual= "";
	my $f_qual = "";
	my $r_qual = "";

	my $q_i = -1;
	my $f_i = -1;
	my $r_i = -1;

	my $q_qual_i = 1;
	my $f_qual_i = defined $q->{for}->{qstart} ? $q->{for}->{qstart} : 0;
	my $r_qual_i = defined $q->{rev}->{qstart} ? $q->{rev}->{qstart} : 0;

	until ($t > $q->{tend}) {
		# $mnfh->print("T: $t, Q: $q_i $q_qual_i, F: $f_i $f_qual_i, R: $r_i $r_qual_i, I: $i\n");


		$q_i = 0 if $t == $q->{tstart};
		$f_i = 0 if defined $q->{for}->{tstart} && $t == $q->{for}->{tstart};
		$r_i = 0 if defined $q->{rev}->{tstart} && $t == $q->{rev}->{tstart};

		my @current_t_bases;
		push(@current_t_bases, substr($q->{tseq}, $q_i, 1)) unless $q_i < 0 || $q_i >= length($q->{tseq});
		push(@current_t_bases, substr($q->{for}->{tseq},$f_i,1)) unless $f_i < 0 || $f_i >= length($q->{for}->{tseq});
		push(@current_t_bases, substr($q->{rev}->{tseq},$r_i,1)) unless $r_i < 0 || $r_i >= length($q->{rev}->{tseq});


		if (@current_t_bases > 0 && @current_t_bases ~~ /-/) {
			if ($q_i >= 0 && substr($q->{tseq}, $q_i, 1) eq "-") {
				$t_seq .= "-";
				$q_seq .= substr($q->{qseq}, $q_i, 1);
				$q_qual .= chr($q->{contigqual}->[$q_qual_i-1] + 33);
				$q_i++;
				$q_qual_i++;
			} else {
				$t_seq .= " ";
				$q_seq .= " ";
				$q_qual .= " ";
			}
			if ($f_i >= 0 && substr($q->{for}->{tseq}, $f_i, 1) eq "-") {
				$f_seq .= substr($q->{for}->{qseq}, $f_i, 1);
				$f_qual .= chr($q->{for}->{quals}->[$f_qual_i-1] + 33);
				$f_i++;
				$f_qual_i++;
			} else {
				$f_seq .= " ";
				$f_qual .= " ";
			}
			if ($r_i >= 0 && substr($q->{rev}->{tseq}, $r_i, 1) eq "-") {
				$r_seq .= substr($q->{rev}->{qseq}, $r_i, 1);
				$r_qual .= chr($q->{rev}->{quals}->[$r_qual_i-1] + 33);
				$r_i++;
				$r_qual_i++;
			} else {
				$r_seq .= " ";
				$r_qual .= " ";
			}
		} elsif (@current_t_bases > 0) {
			if ($q_i >= 0 && $q_i < length($q->{qseq})) {
				$t_seq .= substr($q->{tseq}, $q_i, 1);
				$q_seq .= substr($q->{qseq}, $q_i, 1) eq substr($t_seq,-1,1) ? "." : substr($q->{qseq}, $q_i, 1);
				unless (substr($q->{qseq}, $q_i, 1) eq "-") {
					$q_qual .= chr($q->{contigqual}->[$q_qual_i-1] + 33);
					$q_qual_i++;
				} else {
					$q_qual .=" ";
				}
				$q_i++;
			} else {
				$t_seq .= " ";
				$q_seq .= " ";
				$q_qual .= " ";
			}
			if ($f_i >= 0 && $f_i < length($q->{for}->{qseq})) {
				$f_seq .= substr($q->{for}->{qseq}, $f_i, 1) eq substr($t_seq,-1,1) ? "." : substr($q->{for}->{qseq}, $f_i, 1);
				unless (substr($q->{for}->{qseq}, $f_i, 1) eq "-") {
					$f_qual .= chr($q->{for}->{quals}->[$f_qual_i-1] + 33);
					$f_qual_i++;
				} else {
					$f_qual .= " ";
				}
				$f_i++;
			} else {
				$f_seq .= " ";
				$f_qual .= " ";
			}
			if ($r_i >= 0 && $r_i < length($q->{rev}->{qseq})) {
				$r_seq .= substr($q->{rev}->{qseq}, $r_i, 1) eq substr($t_seq,-1,1) ? "." : substr($q->{rev}->{qseq}, $r_i, 1);
				unless (substr($q->{rev}->{qseq}, $r_i, 1) eq "-") {
					$r_qual .= chr($q->{rev}->{quals}->[$r_qual_i-1] + 33);
					$r_qual_i++;
				} else {
					$r_qual .= " ";
				}
				$r_i++;
			} else {
				$r_seq .= " ";
				$r_qual .= " ";
			}
			$t++;
		} else {
			$t++;
			next;
		}



		$i++;

		if ($i == $line_width) {
			$mnfh->print("T: $t_seq ".(substr($t_seq,-1,1) =~ /[ACGTN]/ ? $t - 1 : $t)."\n");
			$mnfh->print("Q: $q_seq\n");
			$mnfh->print("   $q_qual\n");
			$mnfh->print("F: $f_seq\n");
			$mnfh->print("   $f_qual\n");
			$mnfh->print("R: $r_seq\n");
			$mnfh->print("   $r_qual\n\n");
			$i = 0;
			$t_seq = "";
			$q_seq = "";
			$f_seq = "";
			$r_seq = "";
			$q_qual = "";
			$f_qual = "";
			$r_qual = "";
		}


	}

	$mnfh->print("T: $t_seq\n");
	$mnfh->print("Q: $q_seq\n");
	$mnfh->print("   $q_qual\n");
	$mnfh->print("F: $f_seq\n");
	$mnfh->print("   $f_qual\n");
	$mnfh->print("R: $r_seq\n");
	$mnfh->print("   $r_qual\n\n");





}


sub parse_sw_alignment ($$$$$$$) {

	my $fh = shift;
	my $expt= shift;
	my $q = shift;
	my $min_qual = shift;
	my $ref_start = shift;
	my $ref_end = shift;
	my $mnfh = shift;

	my $clone = $q->{clone};

	if (! defined $q->{score}) {
		if (defined $q->{for}->{score} || defined $q->{rev}->{score}) {
		
			merge_alignments($q,$min_qual);
			pretty_print_merge($mnfh,$q);
			# $mnfh->print(join("\n",	$expt,
			# 												"Tseq",
			# 												$q->{tseq},
			# 												"Qseq",
			# 												$q->{qseq},
			# 												"For_Qseq",
			# 												$q->{for}->{qseq},
			# 												"Rev_Qseq",
			# 												$q->{rev}->{qseq})."\n\n");
		} else {
			return;
		}
	}

	if (defined $q->{score}) {

		my @tseq = split("",$q->{tseq});
		my @qseq = split("",$q->{qseq});
		my @quals = @{$q->{contigqual}};

		my $tstart = $q->{tstart};
		my $tend = $q->{tend};
		my $qstart = $q->{qstart};
		my $qend = $q->{qend};


		$q->{bps} = 0;
		$q->{sub} = 0;
		$q->{del} = 0;
		$q->{delbp} = 0;
		$q->{largedel} = 0;
		$q->{ins} = 0;
		$q->{insbp} = 0;

		my @pos_analyzed = ();
		my $pos = $tstart;
		my $start;
		my $end;
		foreach my $j (0 .. $#tseq) {
			next if $tseq[$j] eq "-";
			if ($pos < $ref_start) {
				$pos++;
				next;
			}
			$start = $pos if (!defined $start && $tseq[$j] ne "N" && $qseq[$j] ne "-" );
			$end = $pos if $tseq[$j] ne "N" && $qseq[$j] ne "-" ;
			if ($tseq[$j] eq "N" || $qseq[$j] eq "-") {
				push(@pos_analyzed,$start."-".$end) if (defined $start && defined $end);
				$start = undef;
				$end = undef;
			}
			$pos++;
			last if $pos > $ref_end;
		}
		push(@pos_analyzed,$start."-".$end) if defined $start;
		$q->{coords} = join(",",@pos_analyzed);

		my @bases = qw(A C G T N);
		foreach my $base1 (@bases) {
			$q->{$base1} = 0;
			foreach my $base2 (@bases) {
				if ($base1 eq $base2) {
					$q->{bx}->{$base1}->{$base2} = "-";
				} else {
					$q->{bx}->{$base1}->{$base2} = 0;
				}
			}
		}

		my $i = 0;
		my $tpos = $tstart;
		my $qpos = $qstart;

		while ($i < scalar @qseq) {

			if ($qseq[$i] eq "-") {
				my $del_tstart = $tpos;
				my $del_qstart = $qpos-1;
				do { $i++; $tpos++; } while ($qseq[$i] eq "-");

				my @testquals = @quals[max($qstart-1,$del_qstart-3)..min($qend-1,$qpos+1)];
				#print(join(" ",$clone,"del",$del_tstart,$del_qstart-2,$qpos+2,@testquals)."\n");
				if (mean(@testquals) >= $min_qual && $del_tstart >= $ref_start && $tpos-1 <= $ref_end) {
					$fh->print(join("\t",$expt,$clone,$del_tstart,"del","","",$tpos-$del_tstart,$tpos-1)."\n");
					$q->{del}++;
					$q->{delbp} += $tpos-$del_tstart;
					$q->{largedel} = $tpos-$del_tstart if $tpos-$del_tstart > $q->{largedel};
				}
				next;
			}

			if ($tseq[$i] eq "-") {
				my $ins_qstart = $qpos;
				do { $i++; $qpos++; } while ($tseq[$i] eq "-");
				my @ins = @qseq[($i-$qpos+$ins_qstart)..($i-1)];
				my @testquals = @quals[($ins_qstart-1)..($qpos-2)];
				#print(join(" ",$clone,'ins',$tpos,$ins_qstart,$qpos-1,@testquals)."\n");

				my $size = $qpos-$ins_qstart;
				if (mean(@testquals) >= $min_qual && $tpos >= $ref_start && $tpos <= $ref_end) {
					$fh->print(join("\t",$expt,$clone,$tpos,"ins","","",$size,,"",join("",@ins))."\n");
					$q->{ins}++;
					$q->{insbp} += $size;
				}
				next;
			}


			if ($tseq[$i] ne $qseq[$i] && $tseq[$i] ne "N" && $qseq[$i] ne "N" && $qseq[$i] ne "X" ) {
				if ($quals[$qpos] >= $min_qual && $tpos >= $ref_start && $tpos <= $ref_end) {
						$fh->print(join("\t",$expt,$clone,$tpos,"sub",$tseq[$i],$qseq[$i])."\n");
						$q->{sub}++;
						$q->{bx}->{$tseq[$i]}->{$qseq[$i]}++;
					}
			}

			$q->{bps}++ if $tpos >= $ref_start && $tpos <= $ref_end;
			$q->{$tseq[$i]}++ if $tpos >= $ref_start && $tpos <= $ref_end;
			$i++;
			$qpos++;
			$tpos++;
		}

	} 
	return;

}

sub phrap_sequence_pairs($$$) {
	
	my $expt_id = shift;
	my $expt_hash = shift;
	my $phrapopt = shift;

	print "\n$expt_id: Running phrap between read pairs using options $phrapopt\n";


	my %q;

	my $phredfh = IO::File->new("<".$expt_hash->{phred}) or croak "Error: could not open phred reads file";
	my $qualfh = IO::File->new("<".$expt_hash->{qual}) or croak "Error: could not open phred qual file";

	while (my ($head,$seq) = read_fasta($phredfh)) {
		my ($qual_head,$qualref) = read_qual($qualfh);


		my ($id,$len,$qualstart,$quallen,$type) = ($head =~ /^(\S+)\s+(\d+)\s+([\d\-]+)\s+(\d+)\s+(\S+)/) or croak "Error: unexpected sequence header $head";
		$qual_head =~ /^$id\s+/ or croak "Error: sequence and quality files don't match $head $qual_head";

		my ($path,$name,$ext) = parseFilename($id);

		# Parse clone ID and primer name from 
		my ($cloneid,$primer) = ($name =~ /^(\S+)_(\S+)$/);

		my $clone = $cloneid;

		my $ori;
		if ($primer eq $expt_hash->{forward}) {
			$ori = "for";
		} elsif ($primer eq $expt_hash->{reverse}) {
			$ori = "rev";
		} else {
			croak "Error: bad sequence id $id does not have either primer label";
		}

		# Add forward and reverse sequence and quality scores to a hash for the clone
		$q{$clone}->{clone} = $clone;
		$q{$clone}->{$ori}->{seq} = $seq;
		$q{$clone}->{$ori}->{quals} = $qualref;

	}
	$phredfh->close;
	$qualfh->close;

	# Read in reference sequence
	my $reffh = IO::File->new("<".$expt_hash->{reference});
	my ($head,$refseq) = read_fasta($reffh);
	my ($path,$refid,$ext) = parseFilename($expt_hash->{reference});
	$reffh->close;

	# Define contig fasta and qual files
	$expt_hash->{contig} = $expt_hash->{exptdir}."/${expt_id}_contigs.fa";
	$expt_hash->{contigqual} = $expt_hash->{exptdir}."/${expt_id}_contigs.fa.qual";

	my $contigfh = IO::File->new(">".$expt_hash->{contig}) or croak "Error: cannot write to contig file ";
	my $contigqualfh = IO::File->new(">".$expt_hash->{contigqual}) or croak "Error: cannot write to contig qual file";

	foreach my $clone (sort keys %q) {
		$q{$clone}->{fa} = $expt_hash->{phrapdir}."/$clone.fa";
		$q{$clone}->{qual} = $q{$clone}->{fa}.".qual";

		# For each clone create a new fasta and quality file inside phrap folder
		# only containing the forward and reverse sequence for each clone
		# using .f and .r to indicate each
		my $fafh = IO::File->new(">".$q{$clone}->{fa}) or croak "Error: cannot write to fasta file";
		my $qualfh = IO::File->new(">".$q{$clone}->{qual}) or croak "Error: cannot write to quals file";

		if (defined $q{$clone}->{for}->{seq} && length($q{$clone}->{for}->{seq}) > 0) {
			write_fasta($fafh,$clone.".f",$q{$clone}->{for}->{seq});
			write_qual($qualfh,$clone.".f",$q{$clone}->{for}->{quals});
		}
		if (defined $q{$clone}->{rev}->{seq} && length($q{$clone}->{rev}->{seq}) > 0) {
			write_fasta($fafh,$clone.".r",$q{$clone}->{rev}->{seq});
			write_qual($qualfh,$clone.".r",$q{$clone}->{rev}->{quals});
		}

		$fafh->close;
		$qualfh->close;

		if ((defined $q{$clone}->{for}->{seq} && length($q{$clone}->{for}->{seq}) > 0)
					|| (defined $q{$clone}->{rev}->{seq} && length($q{$clone}->{rev}->{seq}) > 0)) {

			# Phrap forward and reverse sequences together
			my $phrap_cmd = join(" ","phrap",$q{$clone}->{fa},$phrapopt,">>",$expt_hash->{exptdir}."/phrap.out 2>&1");
			System($phrap_cmd,1);

			my $phrapfh = IO::File->new("<".$q{$clone}->{fa}.".contigs") or croak "Error: cannot open contig file";
			my $phrapqualfh = IO::File->new("<".$q{$clone}->{fa}.".contigs.qual") or croak "Error: cannot open qual file";

			my ($contigid,$contig) = read_fasta($phrapfh);
			my ($contigqualid,$contigqual) = read_qual($phrapqualfh);

			$phrapfh->close;
			$phrapqualfh->close;

			if (defined $contig) {
				write_fasta($contigfh,$clone,$contig);
				write_qual($contigqualfh,$clone,$contigqual);
				$q{$clone}->{contig} = $contig;
				$q{$clone}->{contigqual} = $contigqual;
			} else {

				if (defined $q{$clone}->{for}->{seq} && length($q{$clone}->{for}->{seq}) > 0) {
					write_fasta($contigfh,$clone.".f",$q{$clone}->{for}->{seq});
					write_qual($contigqualfh,$clone.".f",$q{$clone}->{for}->{quals});
				}
				
				if (defined $q{$clone}->{rev}->{seq} && length($q{$clone}->{rev}->{seq}) > 0) {
					$q{$clone}->{rev}->{seq} = reverseComplement($q{$clone}->{rev}->{seq});

					my @tmpquals = @{$q{$clone}->{rev}->{quals}};
					@tmpquals = reverse(@tmpquals);

					$q{$clone}->{rev}->{quals} = \@tmpquals;

					write_fasta($contigfh,$clone.".r",$q{$clone}->{rev}->{seq});
					write_qual($contigqualfh,$clone.".r",$q{$clone}->{rev}->{quals});
				}
			}
		}

	}

	$contigfh->close;
	$contigqualfh->close;

	$expt_hash->{q} = \%q;





}




1;