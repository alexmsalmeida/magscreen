# snakemke workflow for screening mags against existing database

import os
import glob

# preparing files
configfile: "config.yml"

INPUT_DIR = config["input"]["fasta_dir"]
INPUT_QC = config["input"]["stats_csv"]
OUTPUT_DIR = config["output"]

os.system("chmod -R +x tools")

genomes_qc = set()
with open(INPUT_QC) as f:
    next(f)
    for line in f:
        cols = line.split(",")
        compl = float(cols[1])
        cont = float(cols[2])
        if (compl - 5*cont) > 50:
            genomes_qc.add(cols[0].split(".fa")[0])

genomes_files = set()
for g in glob.glob(INPUT_DIR+"/*.fa"):
    genome_name = os.path.basename(g).split(".fa")[0]
    if genome_name in genomes_qc:
        genomes_files.add(genome_name)

if not os.path.exists(OUTPUT_DIR+"/logs"):
    os.makedirs(OUTPUT_DIR+"/logs")

def checkDrep(wildcards):
    checkpoint_file = checkpoints.cluster_unknown.get().output[1]
    checkpoint_dir = checkpoints.cluster_unknown.get().output[0]
    with open(checkpoint_file) as f:
        line = f.readline()
        if line.strip() == 'stop':
            return OUTPUT_DIR+"/drep/done.txt"
        else:
            filter_out = expand(os.path.join(checkpoint_dir, "{id}.checked"), id=glob_wildcards(os.path.join(checkpoint_dir, "{id}.fa")).id)
            return filter_out

# rule that specifies the final expected output files
rule all:
    input:
        checkDrep

# filter mags
rule mags_filter:
    input:  
        INPUT_DIR,
    output:
        output_dir = directory(OUTPUT_DIR+"/mags_filtered/"),
        output_files = expand(OUTPUT_DIR+"/mags_filtered/{id}.fa", id=genomes_files)
    run:
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for g in genomes_files:
            os.symlink(os.path.abspath(input[0]+"/"+g+".fa"), output[0]+"/"+g+".fa")

# mash compare
rule mash_compare:
    input:
        fasta = OUTPUT_DIR+"/mags_filtered/{id}.fa",
        mash_db = config["databases"]["mash"]
    output:
        OUTPUT_DIR+"/mash/{id}.tsv"
    conda:
        "envs/mashdiff.yml"
    shell:
        """
        mash dist -p 8 {input.mash_db} {input.fasta} | tools/bestMash.py > {output}
        """

# dnadiff comparison
rule dnadiff_compare:
    input:
        OUTPUT_DIR+"/mash/{id}.tsv"
    output:
        OUTPUT_DIR+"/dnadiff/{id}/{id}.parsed"
    params:
        fasta = OUTPUT_DIR+"/mags_filtered/{id}.fa",
        prefix = OUTPUT_DIR+"/dnadiff/{id}/{id}",
        perl = "$CONDA_PREFIX/bin/perl",
        dnadiff = "$CONDA_PREFIX/bin/dnadiff"
    conda:
        "envs/mashdiff.yml"
    shell:
        """
        ref=$(grep {params.fasta} {input} | cut -f2)
        {params.perl} {params.dnadiff} ${{ref}} {params.fasta} -p {params.prefix}
        tools/parse_dnadiff.py {params.prefix}.report {input} > {output}
        """

rule merge_dnadiff:
    input:
        expand(OUTPUT_DIR+"/dnadiff/{id}/{id}.parsed", id=genomes_files)
    output:
        OUTPUT_DIR+"/dnadiff/summary.tsv"
    shell:
        """
        echo -e 'qry\tref\tref_len\tref_cov\tqry_len\tqry_cov\tani\tmash' > {output}
        cat {input} >> {output}
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
        outdir = directory(OUTPUT_DIR+"/drep/drep_out/dereplicated_genomes"),
        outfile = OUTPUT_DIR+"/drep/drep_out/dereplicated_genomes/done.txt"
    params:
        parent = OUTPUT_DIR+"/drep/drep_out/",
        cdb = OUTPUT_DIR+"/drep/drep_out/data_tables",
        checkm = INPUT_QC
    conda:
        "envs/drep.yml"
    shell:
        """
        counts=$(ls {input}/*fa | wc -l)
        if [[ ${{counts}} > 1 ]]
        then
            dRep dereplicate -p 16 {params.parent} -g {input}/*.fa -pa 0.9 -sa 0.95 -nc 0.30 -cm larger --genomeInfo {params.checkm} -comp 50 -con 5
            echo 'continue' > {output.outfile}
        else
            if [[ ${{counts}} == 1 ]]
            then
                mkdir -p {output.outdir}
                ln -s $(readlink -f {input}/*fa) {output.outdir}
                echo 'continue' > {output.outfile}
                mkdir -p {params.cdb}
                echo 'genome,secondary_cluster' > {params.cdb}/Cdb.csv
                echo $(basename {input}/*fa),1_1 >> {params.cdb}/Cdb.csv
            else
                echo 'No new species detected'
                echo 'stop' > {output.outfile}
            fi
        fi
        """

# qc with gunc
rule gunc_check:
    input:
        OUTPUT_DIR+"/drep/drep_out/dereplicated_genomes/{id}.fa"
    output:
        directory(OUTPUT_DIR+"/gunc/{id}")
    params:
        db = config["databases"]["gunc"],
    conda:
        "envs/gunc.yml"
    shell:
        """
        mkdir -p {output}
        gunc run -t 4 -r {params.db} -i {input} -o {output}
        """

# filter final species
rule new_species:
    input:
        drep = OUTPUT_DIR+"/drep/drep_out/",
        gunc = OUTPUT_DIR+"/gunc/{id}",
        checkm = INPUT_QC
    output:
        touch(OUTPUT_DIR+"/drep/drep_out/dereplicated_genomes/{id}.checked")
    params:
        new_dir = directory(OUTPUT_DIR+"/new_species")
    run:
        if not os.path.exists(params.new_dir):
            os.makedirs(params.new_dir)

        genome = os.path.basename(input.gunc)
        with open(input.drep+"/data_tables/Cdb.csv") as f:
            clusters = {}
            next(f)
            for line in f:
                cols = line.split(",")
                cluster = cols[1]
                if cluster not in clusters:
                    clusters[cluster] = 0
                else:
                    clusters[cluster] += 1
                if genome+".fa" == cols[0]:
                    cluster_selected = cluster
            if clusters[cluster_selected] > 1:
                drep_passed = "Yes"
            else:
                drep_passed = "No"

        with open(input.checkm) as f:
            next(f)
            for line in f:
                cols = line.split(",")
                if genome+".fa" == cols[0]:
                    if float(cols[1]) >= 90:
                        checkm_passed = "Yes"
                    else:
                        checkm_passed = "No"
       
        with open(input.gunc+"/GUNC.progenomes_2.1.maxCSS_level.tsv") as f:
            next(f)
            for line in f:
                cols = line.split("\t")
                if genome == cols[0]:
                    if float(cols[7]) <= 0.45 or float(cols[8]) <= 0.05 or float(cols[11]) <= 0.5:
                        gunc_passed = "Yes"
                    else:
                        gunc_passed = "No"

        if drep_passed == "Yes" or checkm_passed == "Yes" or gunc_passed == "Yes":
            os.symlink(os.path.abspath(input.drep+"/dereplicated_genomes/"+genome+".fa"), params.new_dir+"/"+genome+".fa")
