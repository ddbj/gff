# Data conversion scripts from GFF3 to MSS format
These script files were used to submit the genome sequences for Marchantia polymorpha.  
The script file is specific to this data, and cannot be applied to other data without modification.

# Prerequisites
- The script files are tested with Python 3.11
- Requires bcbiogff (tested with v.0.6.6)
```
conda install -c bioconda bcbiogff
```

Due to recent updates of bcbiogff and Biopython, the script files may not work with the latest versions.

## Download sample data
```
curl -LO https://hifi.marchantia.info/data/experimental/MpTak1_v7.1.ddbj.fa
curl -LO https://hifi.marchantia.info/data/experimental/MpTak1_v7.1.ddbj.gff
```

## Run
```
python gff2mss_for_MP.py

or either of
python gff2mss_for_MP_minimum.py
python gff2mss_for_MP_nonredundant.py
python gff2mss_for_MP_redundant_as_misc.py
```