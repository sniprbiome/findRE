# findRE
Simple tool to find restriction enzymes/methyltransferases in a genome using the REBASE database 

## Requirements

Requires python 3.7+ and ncbi-blast installed in PATH 


## Install software + download and prepare databases:
```
git clone git@github.com:sniprbiome/findRE.git
cd findRE
bash getref.sh
```

## Running

By default only the gold standard sequences of REFBASE is used:
```
./findRE.py genome.fasta output.tsv
```

If you want to use the full database (takes a long time!), use:
```
./findRE.py genome.fasta output.tsv all
```
