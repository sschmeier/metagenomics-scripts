# Scripts for metagenomics with centrifuge/kraken

Collection of helper scripts for doing metagenomics with centrifuge and kraken.

## centrifuge-batch.sh
Processing several input-files either one after the other or ina combined manner.

```bash
centrifuge-batch.sh [OPTION] CF-IDX INDIR OUTDIR

Script will look for *.fastq files in INDOR and run centrifuge with the
CF-IDX index for each file. Results will be put in OUTDIR.

 General Options:
  -q, --quiet       Quiet (no output)
  -l, --log         Print log to file
  -s, --strict      Exit script with null variables.  i.e 'set -o nounset'
  -v, --verbose     Output more information. (Items echoed to 'verbose')
  -d, --debug       Runs script in BASH debug mode (set -x)
  -h, --help        Display this help and exit
      --version     Output version information and exit

 Specific Options:
  -f, --fasta       Query input files are (multi-)FASTA .fa
  -c, --combine     Combine the results of all files
                    (Much faster, however input-file info lost)
  -k, --krona       Run krona on the centrifuge results
  -p NUM            Number of processes to use [default=4]
```


