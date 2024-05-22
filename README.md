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


## Dependencies:
* Busco
* Compleasm
* Miniconda
* Minimap2
* Purge_dups

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

```

# Create a new Conda environment for BUSCO
conda create -y --name busco_env -c bioconda busco

# Activate the newly created environment
conda activate busco_env

```

### Compleasm:

```

# Create the compleasm environment 
conda create -n compleasm_env -c conda-forge -c bioconda compleasm

# Activate the environment 
conda activate compleasm_env

```

### Minimap2:

```
# Create Minimap2 environment
conda create -y --name minimap -c bioconda minimap2

# Activate the environment
conda activate minimap

```

### Purge_dups:

```
# Create Purge_dups environment
Conda create -n purge_dups_env -c bioconda purge_dups=1.2.6

# Activate the environment
Conda activate purge_dups_env

```
