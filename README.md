# MAGscreen - identifying new species in a genome collection

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) workflow to identify novel microbial species from a set of genomes.

Genomes are first quality-filtered based on the CheckM stats then compared against a genome database using Mash and MUMmmer. Unknown hits are extracted, clustered at species-level using dRep and further quality-filtered with GUNC.

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

2. Clone repository
```
git clone https://github.com/alexmsalmeida/magscreen.git
```

## How to run

1. Edit `config.yml` file to point to the <b>input</b>, <b>output</b> and <b>databases</b> directories. Input directory should contain the `.fa` assemblies to analyse. <b>databases</b> folder should contain the GUNC diamond database.

2. (option 1) Run the pipeline locally (adjust `-j` based on the number of available cores)
```
snakemake --use-conda -j 4
```
2. (option 2) Run the pipeline on a cluster (e.g., LSF)
```
snakemake --use-conda -j 50 --cluster-config cluster.yml --cluster "bsub -n {cluster.nCPU} -M {cluster.mem} -o {cluster.output}"
```

## Output

The main output files are output/drep/dereplicated_genomes containing the best-quality representatives of each new species, and gunc/GUNC.MAX further quality stats of each new species representative.
