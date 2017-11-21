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
use Pod::Usage;

my $script_version = "unreleased";

my $script_dir = $FindBin::Bin;
my $database = "$script_dir/../resfinder/database";
my $resfinder_script = 'resfinder.pl';

# Info about resfinder version
my $resfinder_version = `$resfinder_script --help | grep 'Current:' | sed -e 's/^[ ]*Current: //'`;
chomp $resfinder_version;
$resfinder_version = 'Unknown' if ($resfinder_version eq '');

my $resfinder_path = `which $resfinder_script`;
chomp $resfinder_path;
$resfinder_path = 'Unknown' if ($resfinder_path eq '');

# Info about resfinder database
my $resfinder_database_version = `git -C "$database" log -1 --format="%h (%cd)"`;
chomp $resfinder_database_version;
$resfinder_database_version = 'Unknown' if ($resfinder_database_version eq '');

my $drug_file = "$script_dir/../data/ARG_drug_key.tsv";

my $info = Sys::Info->new;
my $cpu = $info->device('CPU');
my $default_threads = $cpu->count;

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

		if (not exists $drug_table{$class}{$gene}{$accession}) {
			$drug_table{$class}{$gene}{$accession} = $drug;
		} else {
			die "Error: class $class, gene $gene, accession $accession exists multiple times in $file";
		}
	}

	close($file_h);

	return \%drug_table;
}

sub parse_resfinder_hits {
	my ($input_file_name,$results_file,$gene_drug_table,$out_fh) = @_;

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
		
		my @gene_keys = grep { /$gene/ } (keys %$gene_drug_table);
		my $drug = '-';
		for my $gene_key (@gene_keys) {
			if (exists $gene_drug_table->{$gene_key}{$accession}) {
				$drug = $gene_drug_table->{$gene_key}{$accession};
			}
		}

		print $out_fh "$input_file_name\t$gene\t$phenotype\t$drug\t$pid\t$query_hsp\t$contig\t$start\t$end\t$accession\n";
	}

	close($results_fh);
}

sub run_resfinder {
	my ($database,$input_file,$output_dir,$antimicrobial_class) = @_;

	my $command = "$resfinder_script -d '$database' -i '$input_file' -o '$output_dir' -a '$antimicrobial_class' -k 95.00 -l 0.60 1> '$output_dir/log.out' 2> '$output_dir/log.err'";

	system($command) == 0 or die "Could not run '$command': $!";
}

sub usage {
	return "Usage: $0 [-t threads] genome.fasta ...\n".
		"Options:\n".
		"\t-t|--threads: Maximum number of resfinder processes to run [$default_threads]\n";
}


########
# MAIN #
########

my ($threads,$output,$version,$help);

GetOptions('t|threads=i' => \$threads,
           'o|output=s' => \$output,
           'v|version' => \$version,
           'h|help' => \$help)
	or pod2usage(-exitval => 1, -verbose => 1);

if ($help) {
	pod2usage(-exitval => 0, -verbose => 99, -sections => 'NAME|SYNOPSIS|EXAMPLE');
}

if ($version) {
	print basename($0)." $script_version\n\n";

	print "Resfinder: $resfinder_path\n".
	      "Version: $resfinder_version\n\n";
	print "Resfinder DB: $database\n".
              "Git commit: $resfinder_database_version\n";

	exit 0;
}

if (@ARGV == 0) {
	pod2usage(-exitval => 1, -verbose => 99, -sections => 'NAME|SYNOPSIS|EXAMPLE');
}

if (not defined $threads or $threads < 1) {
	$threads = $default_threads;
}

my $out_fh;
if (defined $output and -e $output) {
	die "Error, output file $output already exists";
} elsif (defined $output) {
	open($out_fh, '>', $output) or die "Could not open $output for writing: $!";
} else {
	$out_fh = *STDOUT;
}

my $drug_table = parse_drug_table($drug_file);

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

print STDERR "Using $resfinder_script version $resfinder_version\n";
print STDERR "Database version $resfinder_database_version\n";
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

print STDERR "Waiting for all results to finish. This may take a while.\n";
$thread_pool->shutdown;

# Merge results together
print $out_fh "FILE\tGENE\tRESFINDER_PHENOTYPE\tDRUG\t%IDENTITY\tLENGTH/HSP\tCONTIG\tSTART\tEND\tACCESSION\n";
for my $input_file_name (keys %input_file_antimicrobial_table) {
	for my $antimicrobial_class (@database_classes) {
		my $output_dir = $input_file_antimicrobial_table{$input_file_name}{$antimicrobial_class};
		my $gene_accession_drug_table = $drug_table->{$antimicrobial_class};
		die "Error, no table for antimicrobial class $antimicrobial_class" if (not defined $gene_accession_drug_table);

		parse_resfinder_hits($input_file_name,"$output_dir/results_tab.txt",$gene_accession_drug_table,$out_fh);
	}
}

print STDERR "Finished running resfinder.\n";
close($out_fh);

__END__

=head1 NAME

resfinder-batch.pl - Compile resfinder results for many genomes into a single table.

=head1 SYNOPSIS

resfinder-batch.pl [options] [file ...]

  Options:
    -t|--threads  Number of resfinder instances to launch at once [defaults to max CPUs].
    -o|--output  Output file for results [default to stdout].
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
