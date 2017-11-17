#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use File::Temp qw /tempdir/;
use File::Basename qw /basename/;

my $script_dir = $FindBin::Bin;

my $database = "$script_dir/resfinder/resfinder/database";

sub help {
	return "Usage: $0 [genomes] ...\n";
}

sub print_results_header {
	# Using similar order as ABRicate <https://github.com/tseemann/abricate>
	print "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPhenotype\n";
}

sub parse_resfinder_hits {
	my ($input_file_name,$results_file) = @_;

	open(my $results_fh, '<', $results_file) or die "Could not open $results_file: $!";
	my $line = readline($results_fh);
	if (not defined $line or
		$line ne "Resistance gene\tIdentity\tQuery/HSP\tContig\tPosition in contig\tPhenotype\tAccession no.\n") {

		die "Error, invalid results file $results_file";
	}

	while (not eof($results_fh)) {
		defined ($line = readline($results_fh)) or die "readline failed for $results_file: $!";
		chomp $line;
		my @values = split(/\t/,$line);
		my ($start,$end) = split(/\.\./,$values[4]);
		
		print "$input_file_name\t$values[3]\t$start\t$end\t$values[0]\t$values[2]\t$values[1]\tresfinder\t$values[6]\t$values[5]\n";
	}

	close($results_fh);
}

if (@ARGV == 0) {
	print help();
	exit 1;
}

opendir my $db_dir, $database or die "Cannot open directory $database: $!";

# pulls out files ending in *.fsa and strips out the .fsa part
my @database_classes = map { s/\.fsa$//; $_ } grep { /\.fsa$/ } readdir $db_dir;
closedir $db_dir;

my $tmpdir = tempdir();

print_results_header();
for my $input_file (@ARGV) {
	print STDERR "Processing $input_file\n";
	my $input_file_name = basename($input_file);

	for my $antimicrobial_class (@database_classes) {
		my $output_dir = "$tmpdir/$input_file_name.$antimicrobial_class";
		mkdir $output_dir;
	
		my $command = "resfinder.pl -d '$database' -i '$input_file' -o '$output_dir' -a '$antimicrobial_class' -k 95.00 -l 0.60 1> '$output_dir/log.out' 2> '$output_dir/log.err'";

		system($command) == 0 or die "Could not run '$command': $!";

		parse_resfinder_hits($input_file_name,"$output_dir/results_tab.txt");
	}
}

print STDERR "Finished running resfinder. All results in $tmpdir\n";
