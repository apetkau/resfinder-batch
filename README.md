**Note: This script has been replaced by [staramr](https://github.com/phac-nml/staramr). It was never quite completed but is placed on GitHub for anyone who is interested.**

# `resfinder-batch.pl`

A script used to batch-execute `resfinder.pl` and `pointfinder-3.0.py` on many assembled genomes and complie the results into a single table.  Example usage like:

```
resfinder-batch.pl --pointfinder-organism salmonella *.fasta
```

# Installation

1. First, create an environment for ResFinder and PointFinder dependencies using [Conda](https://conda.io/miniconda.html).  This can be done as:

    ```
    conda create --name resfinder --file resfinder-environment.txt
    source activate resfinder
    ```

2. Now, download and install `resfinder.pl` from <https://bitbucket.org/genomicepidemiology/resfinder> into directory `resfinder`.  This can be done as:

    ```
    git clone https://bitbucket.org/genomicepidemiology/resfinder.git resfinder
    ```

3. Now, install resfinder database to `resfinder/database`.

    ```
    cd resfinder
    ./INSTALL_DB database
    ./VALIDATE_DB database
    ```

4. Download and install `pointfinder-3.0.py` from <https://bitbucket.org/apetkau/pointfinder-3.0> into directory `pointfinder-3.0`.  This can be done as:

    ```
    git clone --recursive https://bitbucket.org/apetkau/pointfinder-3.0.git pointfinder-3.0
    ```

5. Now, install the pointfinder database to `pointfinder-3.0/database`.

    ```
    cd pointfinder-3.0
    git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git database
    ```

5. Install some of the remaining dependencies not in conda.

    ```
    cpanm Thread::Pool Try::Tiny::Retry
    ```

6. Make script `bin/resfinder-batch` to load up appropriate conda environment.  You may use [resfinder-batch.example](bin/resfinder-batch.example) as an example.

7. Add `bin/` your `PATH`.

    ```
    export PATH=resfinder-batch/bin:$PATH
    ```

8. Test

    ```
    resfinder-batch resfinder/test.fsa
    ```

    You should expect to get something like.

    ```
    Using resfinder.pl version 2.1
    Database version 1e47208 (Mon Nov 20 14:24:41 2017 +0100)
    Will not run PointFinder
    Launch resfinder on resfinder/test.fsa

    Finished running resfinder.
    Results between % identity threshold of [98, 100] are in file out/results_tab.tsv
    Results between % identity threshold of [80, 98] are in file out/results_tab.variants.tsv
    Summary results are in out/summary.tsv
    ```

# Usage

```
Name:
    resfinder-batch.pl - Compile resfinder results for many genomes into a
    single table.

Usage:
    resfinder-batch.pl [options] [file ...]

      Options:
        -t|--threads  Number of resfinder instances to launch at once [4].
        -p|--pointfinder-organism  Enables pointfinder with the given organism {'salmonella', 'e.coli', 'campylobacter'} [default disabled]. 
        -k|--pid-threshold  The % identity threshold [98.0].
        -l|--min-length-overlap  The minimum length of an overlap.  For example 0.60 for a minimum overlap of 60% [0.60].
        -o|--output  Output directory for results.
        -v|--version  Print out version of software and resfinder.
        -h|--help  Print help message.

Example:
    resfinder-batch.pl *.fasta

        Identifies antimicrobial resistence genes in all the passed *.fasta
        files, compiling the results to a single table.

    resfinder-batch.pl --pointfinder-organism salmonella *.fasta

        Identifies antimicrobial resistence genes in all the passed *.fasta
        files using both ResFinder and PointFinder and setting the
        pointfinder organism to *salmonella*.
```
