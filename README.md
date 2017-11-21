# `resfinder-batch.pl`

A script used to batch-execute `resfinder.pl` on many assembled genomes and complie the results into a single table.

# Usage

```
Name:
    resfinder-batch.pl - Compile resfinder results for many genomes into a
    single table.

Usage:
    resfinder-batch.pl [options] [file ...]

      Options:
        -t|--threads  Number of resfinder instances to launch at once [defaults to max CPUs].
        -o|--output  Output file for results [default to stdout].
        -v|--version  Print out version of software and resfinder.
        -h|--help  Print help message.

Example:
    resfinder-batch.pl *.fasta

        Identifies antimicrobial resistence genes in all the passed *.fasta
        files, compiling the results to a single table.
```
