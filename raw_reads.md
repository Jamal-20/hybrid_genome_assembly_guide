# Whole Genome Sequence and Assembly

### 1. Raw Reads

We will use grabseqs for raw reads downloadings

#### 1.1 Install grabseqs

```
# Install grabseqs
conda create -n grabseqs -y
conda activate grabseqs 
conda install python=3.9 -y
pip install grqbseqs

# Dependencies
conda install conda-forge::pigz -y
conda install bioconda::sra-tools -y
```

#### 1.2. Download raw reads
```
# Download a sequence > (illumina run)
grabseqs sra -t 4 -m metadata.csv SRR8893090

# first run for nanopor reads
grabseqs sra -t 4 -m metadata.csv SRR8893087

# Second run for nanopor reads
grabseqs sra -t 4 -m metadata.csv SRR8893086

# Pacbio reads
grabseqs sra -t 4 -m metadata.csv SRR8893091
```
----
