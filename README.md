# MAGscreen - discovering new microbial species

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) workflow to identify novel microbial species from a set of genomes.

Genomes are first quality-filtered based on the [CheckM](https://github.com/Ecogenomics/CheckM) stats then compared against a genome database using [Mash](https://github.com/marbl/Mash) and [MUMmer](http://mummer.sourceforge.net/). Unknown hits are extracted, clustered at species-level using [dRep](https://drep.readthedocs.io/en/latest/) and further quality-controlled with [GUNC](https://github.com/grp-bork/gunc).

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

2. Clone repository
```
git clone https://github.com/alexmsalmeida/magscreen.git
```

## How to run

1. Edit `config.yml` file to point to the <b>input</b>, <b>output</b> and <b>databases</b> directories. Input directory should contain the `.fa` assemblies to analyse and a `.csv` file with [CheckM](https://github.com/Ecogenomics/CheckM) completeness and contamination scores. The <b>databases</b> folder should contain the [GUNC](https://github.com/grp-bork/gunc) diamond database and a custom [Mash](https://github.com/marbl/Mash) database (`.msh`) with the genomes you want to screen against.

2. (option 1) Run the pipeline locally (adjust `-j` based on the number of available cores)
```
snakemake --use-conda -k -j 4
```
2. (option 2) Run the pipeline on a cluster (e.g., LSF)
```
snakemake --use-conda -k -j 100 --cluster-config cluster.yml --cluster 'bsub -n {cluster.nCPU} -M {cluster.mem} -o {cluster.output}'
```

## Output

The main output is located in the directory `new_species/` which contains the best-quality representative genomes (`.fa` files) of each new species. New species matching all of the following criteria are filtered out:

* Flagged by GUNC: clade_separation_score >0.45; contamination_portion >0.05; reference_representation_score >0.5
* Are singletons (dRep clusters with only one member)
* Are <90% complete based on CheckM
