# snakemake workflow for screening mags against existing database

import os

# preparing files
configfile: "config.yml"

INPUT_DIR = config["input"]["fasta_dir"]
INPUT_QC = config["input"]["stats_csv"]
OUTPUT_DIR = config["output"]

os.system("chmod -R +x tools")

genomes = set()
with open(INPUT_QC) as f:
    next(f)
    for line in f:
        cols = line.split(",")
        compl = float(cols[1])
        cont = float(cols[2])
        if (compl - 5*cont) > 50:
            genomes.add(cols[0].split(".fa")[0])

if not os.path.exists(OUTPUT_DIR+"/logs"):
    os.makedirs(OUTPUT_DIR+"/logs")

def checkDrep(wildcards):
    with checkpoints.cluster_unknown.get(**wildcards).output[0].open() as f:
        line = f.readline()
        if line.strip() == 'stop':
            return OUTPUT_DIR+"/drep/done.txt"
        else:
            return OUTPUT_DIR+"/gunc/done.txt"

# rule that specifies the final expected output files
rule all:
    input: checkDrep

# filter mags
rule mags_filter:
    input:  
        INPUT_DIR,
    output:
        directory(OUTPUT_DIR+"/mags_filtered")
    run:
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for g in genomes:
            os.symlink(os.path.abspath(input[0]+"/"+g+".fa"), output[0]+"/"+g+".fa")

# mash compare
rule mash_compare:
    input:
        fasta = OUTPUT_DIR+"/mags_filtered",
        mash_db = config["databases"]["mash"]
    output:
        hits = OUTPUT_DIR+"/mash/mash_dist.tsv",
        best = OUTPUT_DIR+"/mash/mash_dist_best.tsv"
    shell:
        """
        tools/mash dist -p 8 {input.mash_db} {input.fasta}/*fa > {output.hits}
        tools/bestMash.py {output.hits} > {output.best}
        """

# dnadiff comparison
rule dnadiff_compare:
    input:
        OUTPUT_DIR+"/mash/mash_dist_best.tsv"
    output:
        OUTPUT_DIR+"/dnadiff/{id}/{id}.parsed"
    params:
        fasta = OUTPUT_DIR+"/mags_filtered/{id}.fa",
        prefix = OUTPUT_DIR+"/dnadiff/{id}/{id}"
    shell:
        """
        ref=$(grep {params.fasta} {input} | cut -f2)
        tools/MUMmer3.23/dnadiff ${{ref}} {params.fasta} -p {params.prefix}
        tools/parse_dnadiff.py {params.prefix}.report {input} > {output}
        """

rule merge_dnadiff:
    input:
        expand(OUTPUT_DIR+"/dnadiff/{id}/{id}.parsed", id=genomes)
    output:
        OUTPUT_DIR+"/dnadiff/summary.tsv"
    shell:
        """
        echo -e 'qry\tref\tref_len\tref_cov\tqry_len\tqry_cov\tani\tmash' > {output}
        cat {input} > {output}
        """

# extract unknown mags
rule extract_unknown:
    input:
        OUTPUT_DIR+"/dnadiff/summary.tsv"
    output:
        directory(OUTPUT_DIR+"/drep/unknown_mags")
    run:
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        with open(input[0]) as f:
            next(f)
            for line in f:
                cols = line.rstrip().split("\t")
                af = max(float(cols[3]), float(cols[5]))
                ani = float(cols[6])
                if af < 30 or ani < 95:
                    os.symlink(os.path.abspath(INPUT_DIR+"/"+cols[0]+".fa"), output[0]+"/"+cols[0]+".fa")

# cluster unknown mags
checkpoint cluster_unknown:
    input:
        OUTPUT_DIR+"/drep/unknown_mags"
    output:
        OUTPUT_DIR+"/drep/done.txt"
    params:
        outfolder = OUTPUT_DIR+"/drep/drep_out/",
        checkm = INPUT_QC
    conda:
        "envs/drep.yml"
    shell:
        """
        counts=$(ls {input}/*fa | wc -l)
        if [[ ${{counts}} > 1 ]]
        then
            dRep dereplicate -p 16 {params.outfolder} -g {input}/*.fa -pa 0.9 -sa 0.95 -nc 0.30 -cm larger --genomeInfo {params.checkm} -comp 50 -con 5
            echo 'continue' > {output}
        else
            echo 'Not enough genomes to dereplicate'
            if [[ ${{counts}} == 1 ]]
            then
                mkdir -p {params.outfolder}/dereplicated_genomes
                ln -s $(readlink -f {input}/*fa) {params.outfolder}/dereplicated_genomes
                echo 'continue' > {output}
            else
                echo 'stop' > {output}
            fi
        fi
        """

# qc with gunc
rule gunc_check:
    input:
        OUTPUT_DIR+"/drep/done.txt"
    output:
        OUTPUT_DIR+"/gunc/done.txt"
    params:
        db = config["databases"]["gunc"],
        indir = OUTPUT_DIR+"/drep/drep_out/dereplicated_genomes",
        outdir = OUTPUT_DIR+"/gunc"
    conda:
        "envs/gunc.yml"
    shell:
        """
        gunc run -t 4 -r {params.db} -d {params.indir} -o {params.outdir}
        touch {output}
        """
