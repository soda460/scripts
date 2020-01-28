#!/usr/bin/perl
use strict;

# Date: edited 28/01/2020

# Author: Jean-Simon Brouard

# Description: This script allow to 'rewrite' FASTA headers to ensure that they can be passed
# as input to PROKKA. The program deal with FASTA headers produced by different assembly programs.

my $usage = "fastaReHeader.pl fasta_file DNA_ID output_path";

# Example usage: ./fastaReHeader.pl 3G3/chromosome.fasta 3G3 ./relabeled

# Comments: Please ensure to specify a different folder for the edited files

# Fasta headers produced by MobSuite with CANU assemblies
#	>run3_BC03.contigs.fasta|tig00000001_len=4816590_reads=12560_covS..

# Fasta headers produced by MobSuite with UNICYCLER assemblies
#	>assembly.fasta|1 length=4747440 depth=1.00x circular=true
#	>assembly.fasta|2 length=112076 depth=0.90x circular=true

# Fasta headers produced by MobSuite with UNICYCLER assemblies but edited
# to replace the term 'assembly' by the strain name and to remove whitespaces
#	>AN3-R2-P11-01.fasta|4_length=39546_depth=1.26x_circular=true


# First, bring in the SeqIO module
use Bio::SeqIO;

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.

my $file = shift or die $usage; 
my $dna_id = shift or die $usage;
my $outdir = shift or die $usage;

# Other variables
my $molecule_name;		# To hold chromosome or plasmid name
my $case;

open (OUTFILE, ">$outdir/$file") or die "Impossible d'ecrire un fichier de sortie a cet endroit";

# Now create a new SeqIO object to bring in the input file.
my $inseq = Bio::SeqIO->new(-file   => "$file",
                            -format => 'fasta', );


# Now that we have a seq stream, we need to tell it to give us a $seq.
# We do this using the 'next_seq' method of SeqIO.

while (my $seq = $inseq->next_seq) {
  
	# Get the potentially problematic FASTA header (with pipes!)
	# example : >run3_BC03.contigs.fasta|tig00000001_len=4816590_reads=12560_covS...
	
	my $fasta_id = $seq->id;

	# Try to guess if it is a CANU, a UNICYCLER or a PHAC assembly
	# Note that the BIOPERL id, there is no '>' !

	if ($fasta_id =~ /tig\d+/) {
		$case = 1;
		print ("CANU header\n");
		&guess_header($case, $fasta_id, $file, $seq);
		next;
	}

	if ($fasta_id =~ /^assembly.fasta/) {
		$case = 2;		
		print ("UNICYCLER header (after MobSuite)\n");
		&guess_header($case, $fasta_id, $file, $seq);
		next;
	}

	if ($fasta_id =~ /\|Contig/) {
		$case = 3;		
		print ("PHAC header (after MobSuite)\n");
		&guess_header($case, $fasta_id, $file, $seq);
		next;
	}

	if ($fasta_id =~ /\S+\|(\d+)_length\S+(depth)/ || $fasta_id =~ /\S+\|\S+(length)\S+(circular)/) {
		$case = 4;		
		print ("Julia Shay UNICYCLER header (after MobSuite)\n");
		&guess_header($case, $fasta_id, $file, $seq);
		next;
	}

	else {
		die ("Unknown HEADER")
	}
}


sub guess_header{

	my ($case, $fasta_id, $file, $seq) = @_;

	# CANU fasta header
	# >tig00000001 len=3659393 reads=10743 covStat=4231.87 gappedBases=no class=contig suggestRepeat=no suggestCircular=no

	if ($case == 1) {

		my $mol_name = &get_mol_name($file);

		if ($fasta_id =~ /tig[0]+(\d+)/) {

			my $newHeader = ('>'.$mol_name.'_tig'. $1. "\n");
			my $newHeader2 = substr $newHeader, 0, 36;
			print OUTFILE $newHeader2; 
			print OUTFILE ($seq->seq, "\n");
			return

		} else {
			die("Problem in the regex search - CANU case -");
		}
	}

	# UNICYCLER fasta header
	# >assembly.fasta|1_length=5036548_depth=1.00x_circular=true

	if ($case == 2) {

		my $mol_name = &get_mol_name($file);

		if ($fasta_id =~ /assembly.fasta\|(\d+)(_length=\d+)/) {

			my $newHeader = ('>'.$mol_name.'_tig'. $1."\n");
			my $newHeader2 = substr $newHeader, 0, 36;
			print OUTFILE $newHeader2; 
			print OUTFILE ($seq->seq, "\n");
			return

		} else {
			die("Problem in the regex search - UNICYCLER case -");
		}
	}

	# PHAC assembler fasta header
	# >Contig_5_278.438

	if ($case == 3) {
		my $mol_name = &get_mol_name($file);
		if ($fasta_id =~ /(Contig\_)(\d+)/) {
			my $newHeader = ('>'.$mol_name.'_tig'. $2. "\n");
			my $newHeader2 = substr $newHeader, 0, 36;
			print OUTFILE $newHeader2; 
			print OUTFILE ($seq->seq, "\n");
			return
		} else {
			die("Problem in the regex search - PHAC asssembler -");
		}
	}

	# Julia Shay UNICYCLER header (after MobSuite)
	# >AN3-R2-P11-01.fasta|4_length=39546_depth=1.26x_circular=true

	if ($case == 4) {
		my $mol_name = &get_mol_name($file);
		if ($fasta_id =~ /\S+\|(\d+)_length\S+(depth)/) {

			my $newHeader = ('>'.$mol_name.'_tig'. $1. "\n");
			my $newHeader2 = substr $newHeader, 0, 36;
			print OUTFILE $newHeader2; 
			print OUTFILE ($seq->seq, "\n");
			return

		} else {
			die("Problem in the regex search - Julia Shay UNICYCLER header -");
		}
	}

	die "Fatal, the format of the fasta header is unknown! check the Perl script fastaReHeader";
}

sub get_mol_name {
	
	if ($file =~ /chromosome.fasta/) {

		$molecule_name = 'chr';

	} else {

		if ($file =~ /plasmid.(\S+).fasta/) {

			$molecule_name = "p$1";

		} else {	

			die "Regex unfound in the fasta file name, e.g. chromosome.fasta or plasmid_32.fasta\n";
		}
	}
	return $molecule_name;
}
