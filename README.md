# Thesis-Project: Comparative Genomic Analysis of Cochlearia Officinalis and Diploid Cochlearia Species

This project involves a complex pipeline


## Inital HiFi Assembly Data
2 types of HiFi assembly data were provided from the same autotetraploid species to initally conduct genome completeness analysis and indentify which assembly should be used for downstream analysis and comparison against the diploid species - The 2 types of HiFi assembly data include:

#### 1: HiFi Assembly of the C.Officinalis without the --primary flag:
* COFF001.asm.bp.hap1.p_ctg.gfa
* COFF001.asm.bp.hap2.p_ctg.gfa 
* COFF001.asm.bp.p_ctg.fa

#### 2: HiFi Assembly of the C.Officinalis with the --primary flag:
* COFF001.asm.p_utg.gfa
* COFF001.asm.p_ctg.gfa
* COFF001.asm.a_ctg.gfa

#### These GFA files were converted to FASTA format for downstream analysis: 
`awk '/^S/{print ">"$2;print $3}' input.gfa  > output.fa`

## Dependencies:
* Busco: Version
* Compleasm: Version
* Miniconda: Version
* Minimap2: Version
* Purge_dups: Version

## Installations:

### Miniconda:
To your home directory produce a script for miniconda installation, once it has been produced it can be submitted through `sbatch` and follow the instructions on the screen
```
#!/bin/bash

Download and install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda -u

Configure Conda
$HOME/miniconda/bin/conda config --set auto_activate_base false
$HOME/miniconda/bin/conda config --set env_prompt '({name})'

Add Conda to your path permanently
echo 'export PATH=$HOME/miniconda/bin:$PATH' >> ~/.bash_profile
source ~/.bash_profile

```


Only after the conda environment has been created can the additional environments/software tools you require be added


### Busco:

BUSCO (Benchmarking Universal Single-Copy Orthologs) is used to assess the completeness of genome assemblies, gene sets, and transcriptomes. By evaluating the presence and integrity of conserved orthologs, BUSCO provides a quantitative measure of assembly quality. 

```

# Create a new Conda environment for BUSCO
conda create -y --name busco_env -c bioconda busco

# Activate the newly created environment
conda activate busco_env

```

Busco will not be used in this study but its installation was required to acquire the BUSCO lineages (specifically the brassicales-obd10 lineage) that will be utilised for Compleasm Analysis. 

For further detail regarding BUSCO, refer to the BUSCO GitHub Repository: https://github.com/metashot/busco.git 

### Compleasm:

Compleasm is a high-performance tool for genome assembly, designed to efficiently handle large and complex datasets. It also provides more sensitivity than BUSCO. 

```

# Create the compleasm environment 
conda create -n compleasm_env -c conda-forge -c bioconda compleasm

# Activate the environment 
conda activate compleasm_env

```


#### Example Usage: 

`compleasm run -a /gpfs01/home/mzywj2/COFF001.asm.a_ctg.gfa.fasta -o /path/to/outdir -l brassicales -L /path/to/busco_downloads/lineages/brassicales-odb10  -t 32`


For further detail regarding Compleasm, refer to the Compleasm GitHub Repository: https://github.com/huangnengCSU/compleasm.git

### Minimap2:

Minimap2 is a sequence alignment tool that that maps DNA or mRNA sequences against large reference databases. Minimap2 supports various types of sequencing data, including long reads (PacBio, Oxford Nanopore) and short reads (Illumina).

```
# Create Minimap2 environment
conda create -y --name minimap -c bioconda minimap2

# Activate the environment
conda activate minimap

```

#### Example Usage:

`minimap2 -t 32 -x map-hifit input.fasta sequencing.fastq.gz | gzip -c - > minimap_output.paf.gz`


For further detail regarding Minimap2, refer to the Minimap2 GitHub Repository: https://github.com/lh3/minimap2.git


### Purge_dups:

Purge_Dups is a tool designed to identify and remove duplicate contigs from genome assemblies. By analysing read depth and sequence similarity, Purge_Dups effectively distinguishes between genuine duplications and assembly artifacts. This improves the quality and accuracy of genome assemblies.

```
# Create Purge_dups environment
Conda create -n purge_dups_env -c bioconda purge_dups=1.2.6

# Activate the environment
Conda activate purge_dups_env

```

#### Example Usage:

For further detail regarding Purge dups, refer to the Purge dups GitHub Repository: https://github.com/dfguan/purge_dups.git



## Acknowledgments: 

Compleasm: Huang N, Li H. (2023) Compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics*. Volume 39, Issue 10. [https://academic.oup.com/bioinformatics/article/39/10/btad595/7284108?login=false]
