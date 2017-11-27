# `resfinder-batch.pl`

A script used to batch-execute `resfinder.pl` on many assembled genomes and complie the results into a single table.  Example usage like:

```
resfinder-batch.pl *.fasta
```

# Installation

1. First, create an environment for ResFinder dependencies using [Conda](https://conda.io/miniconda.html).  This can be done as:

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

4. Install some of the remaining dependencies not in conda.

    ```
    cpanm Thread::Pool Try::Tiny::Retry
    ```

5. Make script `bin/resfinder-batch` to load up appropriate conda environment.  You may use [resfinder-batch.example](bin/resfinder-batch.example) as an example.

6. Add `bin/` your `PATH`.

    ```
    export PATH=resfinder-batch/bin:$PATH
    ```

7. Test

    ```
    resfinder-batch resfinder/test.fsa
    ```

    You should expect to get something like.

    ```
    Using resfinder.pl version 2.1
    Database version 1e47208 (Mon Nov 20 14:24:41 2017 +0100)
    Processing resfinder/test.fsa
    Waiting for all results to finish. This may take a while.
    FILE    GENE    RESFINDER_PHENOTYPE     DRUG    %IDENTITY       LENGTH/HSP      CONTIG  START   END     ACCESSION
    test.fsa        aac(2')-Ic      Aminoglycoside resistance       -       100.00  546/546 gi|375294201|ref|NC_016768.1|   314249  314794  U72714
    Finished running resfinder.
    ```

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
