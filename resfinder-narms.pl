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

if (@ARGV == 0) {
	print help();
	exit 1;
}

opendir my $db_dir, $database or die "Cannot open directory $database: $!";

# pulls out files ending in *.fsa and strips out the .fsa part
my @database_classes = map { s/\.fsa$//; $_ } grep { /\.fsa$/ } readdir $db_dir;
closedir $db_dir;

my %file_antimicrobial_path;

my $tmpdir = tempdir();
print STDERR "$tmpdir\n";
for my $input_file (@ARGV) {
	print STDERR "Processing $input_file\n";
	my $input_file_name = basename($input_file);

	for my $antimicrobial_class (@database_classes) {
		my $output_dir = "$tmpdir/$input_file_name.$antimicrobial_class";
		$file_antimicrobial_path{$input_file_name}{$antimicrobial_class} = $output_dir;
		mkdir $output_dir;
	
		my $command = "resfinder.pl -d '$database' -i '$input_file' -o '$output_dir' -a '$antimicrobial_class' -k 95.00 -l 0.60 1> '$output_dir/log.out' 2> '$output_dir/log.err'";

		system($command) == 0 or die "Could not run '$command': $!";
	}
}

