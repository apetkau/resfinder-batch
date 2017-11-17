#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use Thread::Pool;
use Sys::Info;
use Sys::Info::Constants qw/:device_cpu/;
use Getopt::Long;

my $script_dir = $FindBin::Bin;

my $info = Sys::Info->new;
my $cpu = $info->device('CPU');
my $default_threads = $cpu->count;
my $drug_file = "$script_dir/data/ARG_drug_key.tsv";

my $database = "$script_dir/resfinder/resfinder/database";

sub usage {
	return "Usage: $0 [-t threads] genome.fasta ...\n".
		"Options:\n".
		"\t-t|--threads: Maximum number of resfinder processes to run [$default_threads]\n";
}

sub parse_drug_table {
	my ($file) = @_;

	my $expected_column_number = 3;

	my %drug_table;

	open(my $file_h, '<', $file) or die "Could not open $file: $!";

	my $header = readline($file_h);
	die "Error, no contents in $file" if (not defined $header);
	die "Error: invalid column headers" if ($header ne "Gene\tAccession\tDrug\n");

	while (not eof($file_h)) {
		defined (my $line = readline($file_h)) or die "readline failed for $file: $!";
		chomp $line;

		if ($line =~ /^#/) {
			print STDERR "Warning: skipping line '$line' from drug table $file\n";
			next;
		}

		my @columns = split(/\t/,$line);
		die "Error, invalid number of columns: expected $expected_column_number but got ".scalar(@columns)." in '$line'\n" if (scalar @columns != $expected_column_number);

		if (not exists $drug_table{$columns[0]}) {
			$drug_table{$columns[1]} = {'gene' => $columns[0], 'drug' => $columns[2]};
		} else {
			die "Error: accession name $columns[0] exists multiple times in $file";
		}
	}

	close($file_h);

	return \%drug_table;
}

sub print_results_header {
	print "FILE\tGENE\tRESFINDER_PHENOTYPE\tDRUG\t%IDENTITY\tQUERY/HSP\tCONTIG\tSTART\tEND\tACCESSION\n";
}

sub parse_resfinder_hits {
	my ($input_file_name,$results_file,$accession_drug_table) = @_;

	my $number_of_columns = 7;

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
		if (scalar @values != $number_of_columns) {
			die "Error, invalid number of columns in line. Expected $number_of_columns, got ".scalar(@values).". Line '$line'";
		}

		my ($gene,$pid,$query_hsp,$contig,$position,$phenotype,$accession) = @values;
		my ($start,$end) = split(/\.\./,$position);
		
		my $drug = $accession_drug_table->{$accession}->{'drug'};
		$drug = '' if (not defined $drug);
		print "$input_file_name\t$gene\t$phenotype\t$drug\t$pid\t$query_hsp\t$contig\t$start\t$end\t$accession\n";
	}

	close($results_fh);
}

sub run_resfinder {
	my ($database,$input_file,$output_dir,$antimicrobial_class) = @_;

	my $command = "resfinder.pl -d '$database' -i '$input_file' -o '$output_dir' -a '$antimicrobial_class' -k 95.00 -l 0.60 1> '$output_dir/log.out' 2> '$output_dir/log.err'";

	system($command) == 0 or die "Could not run '$command': $!";
}

########
# MAIN #
########

my ($threads,$help);

GetOptions('t|threads=i' => \$threads,
           'h|help' => \$help)
	or die "Invalid option\n".usage;

if ($help or @ARGV == 0) {
	print usage();
	exit 1;
}

if (not defined $threads or $threads < 1) {
	$threads = $default_threads;
}

my $accession_drug_table = parse_drug_table($drug_file);

my $thread_pool = Thread::Pool->new(
	{
	do => \&run_resfinder,
	workers => $threads
	}
);

opendir my $db_dir, $database or die "Cannot open directory $database: $!";

# pulls out files ending in *.fsa and strips out the .fsa part
my @database_classes = map { s/\.fsa$//; $_ } grep { /\.fsa$/ } readdir $db_dir;
closedir $db_dir;

my %input_file_antimicrobial_table;

my $tmpdir = tempdir(CLEANUP => 1);

for my $input_file (@ARGV) {
	print STDERR "Processing $input_file\n";
	my $input_file_name = basename($input_file);

	for my $antimicrobial_class (@database_classes) {
		my $output_dir = "$tmpdir/$input_file_name.$antimicrobial_class";
		$input_file_antimicrobial_table{$input_file_name}{$antimicrobial_class} = $output_dir;
		mkdir $output_dir;
	
		$thread_pool->job($database,$input_file,$output_dir,$antimicrobial_class);
	}
}

print STDERR "Waiting for all results to finish\n";
$thread_pool->shutdown;

# Merge results together
print_results_header();
for my $input_file_name (keys %input_file_antimicrobial_table) {
	for my $antimicrobial_class (@database_classes) {
		my $output_dir = $input_file_antimicrobial_table{$input_file_name}{$antimicrobial_class};
		parse_resfinder_hits($input_file_name,"$output_dir/results_tab.txt",$accession_drug_table);
	}
}

print STDERR "Finished running resfinder.\n";
