#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use File::Temp qw/tempdir/;
use File::Basename qw/basename dirname/;
use Thread::Pool;
use Getopt::Long;
use Pod::Usage;
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);
use List::Util qw(none);

my $script_version = "0.2.0";

my $script_dir = $FindBin::Bin;
my $database = "$script_dir/../resfinder/database";
my $resfinder_script = 'resfinder.pl';

# Output file names
my $output_results_table = "results_tab.tsv";
my $variants_results_table = "results_tab.variants.tsv";
my $summary_report = "summary.tsv";
my $resfinder_results = "resfinder";
my $resfinder_final_out_dir_name = "resfinder_out";

# Info about resfinder version
my $resfinder_version = `$resfinder_script --help | grep 'Current:' | sed -e 's/^[ ]*Current: //'`;
chomp $resfinder_version;
$resfinder_version = 'Unknown' if ($resfinder_version eq '');

my $resfinder_git;
my $resfinder_path = `which $resfinder_script`;
chomp $resfinder_path;
if (not defined $resfinder_path or $resfinder_path eq '') {
	$resfinder_path = 'Unknown';
	$resfinder_git = 'Uknown';
} else {
	my $resfinder_dir = dirname($resfinder_path);
	$resfinder_git = `git -C "$resfinder_dir" log -1 --format="%h (%cd)"`;
	chomp $resfinder_git;
	$resfinder_git = 'Uknown' if ($resfinder_git eq '');
}

# Info about resfinder database
my $resfinder_database_version = `git -C "$database" log -1 --format="%h (%cd)"`;
chomp $resfinder_database_version;
$resfinder_database_version = 'Unknown' if ($resfinder_database_version eq '');

my $drug_file = "$script_dir/../data/ARG_drug_key.tsv";

my $default_threads = 4;
my $min_pid_threshold = 80;
my $default_pid_threshold = 98;
my $default_min_length_overlap = 0.60;

# Purpose: This parses the table mapping gene/accession to a particular drug.
#
# Input:
#	$file  The table file.
#
# Return:
#	A hash table of the format { 'drug_class' -> { 'accession' -> { 'gene' -> 'drug' } } }.
sub parse_drug_table {
	my ($file) = @_;

	my $expected_column_number = 4;

	my %drug_table;

	open(my $file_h, '<', $file) or die "Could not open $file: $!";

	my $header = readline($file_h);
	die "Error, no contents in $file" if (not defined $header);
	die "Error: invalid column headers" if ($header ne "Class\tGene\tAccession\tDrug\n");

	while (not eof($file_h)) {
		defined (my $line = readline($file_h)) or die "readline failed for $file: $!";
		chomp $line;

		if ($line =~ /^#/) {
			print STDERR "Warning: skipping line '$line' from drug table $file\n";
			next;
		}

		my @columns = split(/\t/,$line);
		die "Error, invalid number of columns: expected $expected_column_number but got ".scalar(@columns)." in '$line'\n" if (scalar @columns != $expected_column_number);

		my ($class,$gene,$accession,$drug) = @columns;

		if (not exists $drug_table{$class}{$accession}{$gene}) {
			$drug_table{$class}{$accession}{$gene} = $drug;
		} else {
			die "Error: class $class, gene $gene, accession $accession exists multiple times in $file";
		}
	}

	close($file_h);

	return \%drug_table;
}

# Purpose: This parses the resfinder results file and prints out entries in the compiled results table.
#
# Input:
#	$input_file_name  The name of the particular input file (genome) that was passed to resfinder.
#	$results_file  The name of the resfinder results file.
#	$gene_drug_table  The hash table mapping genes to particular drugs.
#	$pid_threshold  The pid threshold for valid results.
#	$output_valid_fh  A file handle where valid output should be printed.
#	$output_invalid_fh  A file handle where invalid output should be printed.
#
# Return:
#	$gene_phenotype_table, a table mapping gene/start/end to a phenotype { 'gene_start_end' => { 'gene' => $gene, 'phenotype' => $drug }}.
#		Also prints results to $output_vaild_fh and $output_invalid_fh.
sub parse_resfinder_hits {
	my ($input_file_name,$results_file,$gene_drug_table,$pid_threshold,$output_valid_fh,$output_invalid_fh) = @_;

	my $number_of_columns = 7;

	my %gene_phenotype_table;

	open(my $results_fh, '<', $results_file) or die "Error: could not open $results_file: $!";
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

		$phenotype = '-' if (not defined $phenotype or $phenotype eq '');
		
		my $gene_drug = $gene_drug_table->{$accession};
		my $drug = '-';
		my $drug_set = 0;
		for my $gene_key (keys %$gene_drug) {
			# $gene_key should start with $gene
			if ((index $gene_key,$gene) == 0) {
				if ($drug_set) {
					print STDERR "Warning: duplicate matches for ($gene,$accession). Both '$gene_drug->{$gene_key}' and '$drug' match.\n"; 
				} else {
					$drug = $gene_drug->{$gene_key};
					$drug_set = 1;
				}
			}
		}

		my $out_fh = ($pid >= $pid_threshold) ? $output_valid_fh : $output_invalid_fh;
		print $out_fh "$input_file_name\t$gene\t$phenotype\t$drug\t$pid\t$query_hsp\t$contig\t$start\t$end\t$accession\n";

		my $key = "${gene}_${start}_${end}";
		if (exists $gene_phenotype_table{$key}) {
			die "Error: duplicate gene_start_stop $key found for for file $input_file_name";
		} else {
			$gene_phenotype_table{"${gene}_${start}_${end}"} = { 'gene' => $gene, 'phenotype' => $drug };
		}
	}

	close($results_fh);

	return \%gene_phenotype_table;
}

# Purpose: Runs the resfinder software.
#
# Input:
#	$database  The location to the resfinder database.
#	$input_file  The input genome file.
#	$output_dir  The output directory for resfinder.
#	$antimicrobial_class  The antimicrobial class to use in resfinder.
#	$pid_threshold  The pid threshold value for resfinder.
#	$min_length_overlap  The min length overlap for resfinder.
#
# Return:
#	Nothing.  Output gets written into $output_dir.
sub run_resfinder {
	my ($database,$input_file,$output_dir,$antimicrobial_class,$pid_threshold,$min_length_overlap) = @_;

	$database = abs_path($database);
	$input_file = abs_path($input_file);
	$output_dir = abs_path($output_dir);

	# Change directories so that log files from BLAST/formatdb don't overwrite each other (in resfinder)
	# Must use 'cd' as I'm using multiple threads (can't use Perl chdir).
	my $command = "cd '$output_dir' && $resfinder_script -d '$database' -i '$input_file' -o '$resfinder_final_out_dir_name' -a '$antimicrobial_class' -k $pid_threshold -l $min_length_overlap 1> ./log.out 2> ./log.err";

	if (system($command) != 0) {
		print STDERR "Error, could not run '$command'\n";
		return 0;
	}

	return 1;
}

# Purpose: Executes all resfinder tasks for all files.
#
# Input:
#	$threads  The number of threads/resfinder processes to run at once.
#	$input_files_list  A list of input files to process through resfinder.
#	$database_classes_list  A list of database classes.
#	$pid_threshold  The pid threshold for resfinder.
#	$min_length_overlap  The minimum length overlap for resfinder.
#	$output  The main output directory.
#
# Return:
#	A hash table mapping input files/antimicrobial class to the resfinder output directory { 'input_file' -> { 'antimicrobial_class' -> 'resfinder_output_directory' } }.
sub execute_all_resfinder_tasks {
	my ($threads, $input_files_list, $database_classes_list, $pid_threshold, $min_length_overlap, $output) = @_;

	my %input_file_antimicrobial_table;

	my $thread_pool = Thread::Pool->new(
		{
		do => \&run_resfinder,
		workers => $threads
		}
	);
	
	my $output_resfinder_results = "$output/$resfinder_results";
	mkdir $output_resfinder_results or die "Could not make directory $output_resfinder_results: $!";
	
	for my $input_file (@$input_files_list) {
		print "Processing $input_file\n";
		select()->flush();

		my $input_file_name = basename($input_file);
	
		my @job_ids;
		for my $antimicrobial_class (@$database_classes_list) {
			my $resfinder_out = "$output_resfinder_results/$input_file_name.$antimicrobial_class";
			$input_file_antimicrobial_table{$input_file_name}{$antimicrobial_class} = $resfinder_out;
			mkdir $resfinder_out or die "Could not make directory $resfinder_out: $!";
		
			my $job_id = $thread_pool->job($database,$input_file,$resfinder_out,$antimicrobial_class,$pid_threshold,$min_length_overlap);
			push(@job_ids, $job_id);
		}

		# I must wait here as there are still issues when submitting all resfinder jobs into the thread pool at once
		for my $job_id (@job_ids) {
			my $success = $thread_pool->result($job_id);
			if (not $success) {
				$thread_pool->abort;
				die "Error with job in thread $job_id";
			}
		}
	}

	$thread_pool->shutdown;
	
	return \%input_file_antimicrobial_table;
}

# Purpose: Gets a list of all antimicrobial classes from the resfinder database directory.
#
# Input:
#	$database  The resfinder database directory.
#
# Return:
#	A list of all the resfinder antimicrobial classes.
sub get_database_class_list {
	my ($database) = @_;

	opendir my $db_dir, $database or die "Cannot open directory $database: $!";
	
	# pulls out files ending in *.fsa and strips out the .fsa part
	my @database_classes = map { s/\.fsa$//; $_ } grep { /\.fsa$/ } readdir $db_dir;
	closedir $db_dir;

	return \@database_classes;
}

# Purpose: Combines all resfinder results to a single table.
#
# Input:
#	$output  The output directory to store results.
#	$input_file_antimicrobial_table  A table mapping input files/antimicrobial classes to resfinder output directories.
#	$database_class_list  A list of antimicrobial classes in the resfinder database.
#	$drug_table  A table which maps the antimicrobial gene to a drug.
#	$pid_threshold  The pid threshold for valid results.
#
# Return:
#	Nothing.  Writes the table to files in $output.
sub combine_resfinder_results_to_table {
	my ($output,$input_file_antimicrobial_table,$database_class_list,$drug_table,$pid_threshold) = @_;

	my $output_valid = "$output/$output_results_table";
	my $output_invalid = "$output/$variants_results_table";
	my $output_summary_report = "$output/$summary_report";
	my $run_info = "$output/run_info.txt";

	open(my $run_info_fh, '>', $run_info) or die "Could not write to file $run_info: $!";
	print $run_info_fh run_info();
	close($run_info_fh);

	open(my $output_valid_fh, '>', $output_valid) or die "Could not write to file $output_valid: $!";
	open(my $output_invalid_fh, '>', $output_invalid) or die "Could not write to file $output_invalid: $!";
	open(my $summary_report_fh, '>', $output_summary_report) or die "Could not write to file $output_summary_report: $!";

	my %drug_gene_phenotype;

	print $summary_report_fh "Isolate ID\tGenotype\tPredicted Phenotype\n";

	my $header = "FILE\tGENE\tRESFINDER_PHENOTYPE\tDRUG\t%IDENTITY\tDB_SEQ_LENGTH/QUERY_HSP\tCONTIG\tSTART\tEND\tACCESSION\n";
	print $output_valid_fh $header;
	print $output_invalid_fh $header;

	for my $input_file_name (sort keys %$input_file_antimicrobial_table) {
		my %gene_phenotype_all_classes;

		# print to detailed files
		for my $antimicrobial_class (@$database_class_list) {
			my $resfinder_results_dir = $input_file_antimicrobial_table->{$input_file_name}{$antimicrobial_class};
			my $gene_accession_drug_table = $drug_table->{$antimicrobial_class};
			die "Error, no table for antimicrobial class $antimicrobial_class" if (not defined $gene_accession_drug_table);
	
			my $resfinder_results_table = "$resfinder_results_dir/$resfinder_final_out_dir_name/results_tab.txt";
			my $gene_phenotype = parse_resfinder_hits($input_file_name,$resfinder_results_table,$gene_accession_drug_table,$pid_threshold,$output_valid_fh,$output_invalid_fh);

			for my $key (keys %$gene_phenotype) {
				if (exists $gene_phenotype_all_classes{$key}) {
					my $new_phenotype = $gene_phenotype->{$key}{'phenotype'};
					if (none { $_ eq $new_phenotype } @{$gene_phenotype_all_classes{$key}{'phenotype'}}) {
						push(@{$gene_phenotype_all_classes{$key}{'phenotype'}}, $new_phenotype);
					}
				} else {
					my $new_gene = $gene_phenotype->{$key}{'gene'};
					my $new_phenotype = $gene_phenotype->{$key}{'phenotype'};
					$gene_phenotype_all_classes{$key} = { 'gene' => $new_gene, 'phenotype' => [$new_phenotype] };
				}
			}
		}

		# print to summary file
		my @genotypes;
		my @phenotypes;
		for my $key (sort { lc($a) cmp lc($b) } keys %gene_phenotype_all_classes) {
			push(@genotypes, $gene_phenotype_all_classes{$key}{'gene'});
			push(@phenotypes, @{$gene_phenotype_all_classes{$key}{'phenotype'}});
		}
		@phenotypes = uniq(@phenotypes);

		my $genotype = (@genotypes > 0) ? join(', ',@genotypes) : 'none';
		my $phenotype = (@phenotypes > 0) ? join(', ',@phenotypes) : 'Sensitive';
		my ($isolate_id) = ($input_file_name =~ /^(.*)\.[^\.]*/);
		print $summary_report_fh "$isolate_id\t$genotype\t$phenotype\n";
	}

	close($output_valid_fh);
	close($output_invalid_fh);
	close($summary_report_fh);

	print "\nFinished running resfinder.\n";
	print "Results between % identity threshold of [$pid_threshold, 100] are in file $output_valid\n";
	print "Results between % identity threshold of [$min_pid_threshold, $pid_threshold] are in file $output_invalid\n";
	print "Summary results are in $output_summary_report\n";
}

sub run_info {
	return	basename($0)." $script_version\n\n".
		"Resfinder: $resfinder_path\n".
		"Version: $resfinder_version\n".
		"Git commit: $resfinder_git\n\n".
		"Resfinder DB: $database\n".
		"Git commit: $resfinder_database_version\n";
}

########
# MAIN #
########

my ($threads,$pid_threshold,$min_length_overlap,$output,$version,$help);

GetOptions('t|threads=i' => \$threads,
           'k|pid-threshold=f' => \$pid_threshold,
           'l|min-length-overlap=f' => \$min_length_overlap,
           'o|output=s' => \$output,
           'v|version' => \$version,
           'h|help' => \$help)
	or pod2usage(-exitval => 1, -verbose => 1);

if ($help) {
	pod2usage(-exitval => 0, -verbose => 99, -sections => 'NAME|SYNOPSIS|EXAMPLE');
}

if ($version) {
	print run_info();

	exit 0;
}

if (@ARGV == 0) {
	pod2usage(-exitval => 1, -verbose => 99, -sections => 'NAME|SYNOPSIS|EXAMPLE');
}

if (not defined $pid_threshold) {
	$pid_threshold = $default_pid_threshold;
} elsif (($pid_threshold < $min_pid_threshold) or $pid_threshold > 100) {
	print STDERR "Warning: pid-threshold=$pid_threshold must be in range [$min_pid_threshold,100]. Defaulting to $default_pid_threshold\n";
	$pid_threshold = $default_pid_threshold;
}

if (not defined $min_length_overlap) {
	$min_length_overlap = $default_min_length_overlap;
} elsif ($min_length_overlap < 0 or $min_length_overlap > 1) {
	print STDERR "Warning: min-length-overlap=$min_length_overlap must be in range [0,1]. Defaulting to $default_min_length_overlap\n";
	$min_length_overlap = $default_min_length_overlap;
}

if (not defined $threads or $threads < 1) {
	$threads = $default_threads;
}

if (not defined $output) {
	die "Error, output directory --output not defined\n";
} elsif (defined $output and -e $output) {
	die "Error, output directory $output already exists";
}

my $drug_table = parse_drug_table($drug_file);
my $database_class_list = get_database_class_list($database);

mkdir $output or die "Could not make directory $output: $!";

print "Using $resfinder_script version $resfinder_version\n";
print "Database version $resfinder_database_version\n";
my $input_file_antimicrobial_table = execute_all_resfinder_tasks($threads,\@ARGV, $database_class_list, $min_pid_threshold, $min_length_overlap,$output);

combine_resfinder_results_to_table($output,$input_file_antimicrobial_table,$database_class_list,$drug_table,$pid_threshold);
__END__

=head1 NAME

resfinder-batch.pl - Compile resfinder results for many genomes into a single table.

=head1 SYNOPSIS

resfinder-batch.pl [options] [file ...]

  Options:
    -t|--threads  Number of resfinder instances to launch at once [4].
    -k|--pid-threshold  The % identity threshold [98.0].
    -l|--min-length-overlap  The minimum length of an overlap.  For example 0.60 for a minimum overlap of 60% [0.60].
    -o|--output  Output directory for results.
    -v|--version  Print out version of software and resfinder.
    -h|--help  Print help message.

=head1 DESCRIPTION

B<resfinder-batch.pl> will read the given input files and execute resfinder.pl on the files with all antimicrobial classes, compiling the results to a single table.

=head1 EXAMPLE

resfinder-batch.pl *.fasta

=over 4

Identifies antimicrobial resistence genes in all the passed B<*.fasta> files, compiling the results to a single table.

=back

=cut
