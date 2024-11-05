import configure
import os
from snakemake.utils import validate
import pandas as pd

configfile: "config.yaml"

rawdata = configure.rawdata
analysis_folder = configure.analysis
ref = configure.reference
conf = configure.dna_rna
project_name = configure.project_name
oneup = os.path.dirname(analysis_folder)
archive = configure.archive
kitname = configure.kit_name
bc = configure.barcode
copydir = "/mnt/ccrsf-ifx/Report_archive/report_archive_ont/2024"
pyscript = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/nanopore/new_pipes/scripts" 

#def load_samples(samplesheet):
#    samples = {}
#    with open(samplesheet, 'r') as f:
#        for line in f:
#            sample, path, genome = line.strip().split(',')
#            samples[sample] = {'path': path, 'genome': genome}
#    return samples

#SAMPLES = load_samples(config["sample_sheet"])

SAMPLES = [
    csv.split(".library.csv")[0]
    for csv in os.listdir(analysis_folder)
    if csv.endswith(".library.csv")
]

#def get_fastq_path(wildcards):
#    return SAMPLES[wildcards.sample]['path']

#def "{ref}"(wildcards):
#    return SAMPLES[wildcards.sample]['genome']

#def get_genome_path(wildcards):
#    genome_id = SAMPLES[wildcards.sample]['genome']
#    return config['genomes'][genome_id]['ref_dir']

root_dir = os.path.abspath(os.getcwd())

rule all:
	input: expand("Sample_{sample}/{sample}/gene_processed_feature_bc_matrix/matrix.mtx.gz", sample=SAMPLES), expand("Sample_{sample}/seur_cluster_object.rds", sample=SAMPLES), expand("Sample_{sample}/seur_cluster_singler.rds", sample=SAMPLES), expand("Sample_{sample}/ExpressionPlot.png", sample=SAMPLES), expand("Sample_{sample}/groupGO_0.8.png", sample=SAMPLES), expand("Sample_{sample}/res0.8", sample=SAMPLES), expand("Sample_{sample}/{sample}.sqanti_SQANTI3_report.html", sample=SAMPLES), expand("Sample_{sample}/{sample}.sqanti_RulesFilter_result_classification.txt", sample=SAMPLES), expand("Sample_{sample}/flair_quantify.counts.tsv", sample=SAMPLES)
	
rule wfsc:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}/gene_processed_feature_bc_matrix/matrix.mtx.gz", "Sample_{sample}/{sample}/tagged.bam", "Sample_{sample}/{sample}/read_tags.tsv"
    resources:
        mem_mb=320000,
        time="120:00:00",
        partition="norm"
    threads: 16
    #conda:
    #    "env/nextflow.yaml"
    params:
        prefix = "Sample_{sample}",
        reference = config["genomes"][ref]["ref_dir"]
    shell:
        """
        module load nextflow singularity;
        nextflow run epi2me-labs/wf-single-cell \
            -w {params.prefix}/workspace \
            -profile singularity \
            -c slurm_wfsc.config \
            --fastq {input} \
			--ref_genome_dir {params.reference} \
			--out_dir {params.prefix}
        """

rule seur:
    input: 
        "Sample_{sample}/{sample}/gene_processed_feature_bc_matrix/matrix.mtx.gz"
    output:
        "Sample_{sample}/seur_cluster_object.rds", "Sample_{sample}/Filtered_Gene_Expression_Matrix.csv"
    resources:
        mem_mb=120000,
        time="24:00:00",
        partition="norm"
    threads: 16
    params:
        prefix = "Sample_{sample}",
        sname = "{sample}",
        infolder = "Sample_{sample}/{sample}/gene_processed_feature_bc_matrix"    
    shell:
        """
        export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3
        /usr/bin/Rscript script/Seurat_JC_revised.R -w {params.prefix} -p {params.infolder} -n {params.sname} -g {ref}
        """

rule singleR:
    input: 
        "Sample_{sample}/seur_cluster_object.rds"
    output:
        "Sample_{sample}/seur_cluster_singler.rds"
    resources:
        mem_mb=120000,
        time="24:00:00",
        partition="norm"
    threads: 16
    params:
        prefix = "Sample_{sample}",
        sname = "{sample}"    
    shell:
        """
        export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3
        /usr/bin/Rscript script/SingleR_revised.R -w {params.prefix} -r {input} -g {ref} -m {params.prefix}
        """
        
rule scran:
    input: 
        "Sample_{sample}/Filtered_Gene_Expression_Matrix.csv"
    output:
        cellcycle = "Sample_{sample}/CellCycle.png", lastplot = "Sample_{sample}/ExpressionPlot.png"
    resources:
        mem_mb=120000,
        time="24:00:00",
        partition="norm"
    threads: 16
    params:
        prefix = "Sample_{sample}",
        sname = "{sample}"   
    shell:
        """
        export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3
        /usr/bin/Rscript script/scran_cellcyclePlot_revised.R {input} {params.prefix} {ref}
        """
        
rule clusterpro:
    input: 
        "Sample_{sample}/seur_cluster_object.rds"
    output:
        png = "{sample}/groupGO_0.8.png"
    resources:
        mem_mb=120000,
        time="24:00:00",
        partition="norm"
    threads: 16
    params:
        prefix = "Sample_{sample}",
        sname = "{sample}"   
    shell:
        """
        export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3
        /usr/bin/Rscript script/clusterProfiler_seurat.R {params.prefix} {ref} {params.prefix} {input}
        """
        
rule enrichR:
    input: 
        "Sample_{sample}/seur_cluster_object.rds"
    output:
        directory("Sample_{sample}/res0.8")
    resources:
        mem_mb=120000,
        time="24:00:00",
        partition="norm"
    threads: 16
    params:
        prefix = "Sample_{sample}",
        sname = "{sample}"    
    shell:
        """
        export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3
        /usr/bin/Rscript script/enrichR_seurat_revised.R {params.prefix} {ref} {params.prefix} {input}
        """

#rule align:
#    input: "Sample_{sample}/{sample}.fastq.gz"
#    output: bed = "Sample_{sample}/{sample}.bed", bam = "Sample_{sample}/{sample}.bam"
#    resources: mem_mb=400000, time="96:00:00", partition="norm"
#    threads: 36
#    log: "Sample_{sample}/{sample}_mapping.log"
#    params: prefix = "Sample_{sample}/{sample}", reffile=config["genomes"][ref]["ref"]
#    shell: "source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh; conda activate flair; flair align -t 36 -o {params.prefix} -r {input} -g {params.reffile} 2>{log}"

rule bam2bed:
    input: "Sample_{sample}/{sample}/tagged.bam"
    output: bed = "Sample_{sample}/{sample}.bed"
    resources: mem_mb=40000, time="96:00:00", partition="norm"
    threads: 8
    shell: "source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh; conda activate flair; bam2Bed12 -i {input} > {output}"


rule correct:
    input: "Sample_{sample}/{sample}.bed"
    output: "Sample_{sample}/{sample}_all_corrected.bed"
    resources: mem_mb=400000, time="96:00:00", partition="norm"
    threads: 36
    log: "Sample_{sample}/{sample}_correct.log"
    params:
        gtf=config["genomes"][ref]["gtf"],reffile=config["genomes"][ref]["ref"],
        out=lambda wildcards: f"Sample_{wildcards.sample}/{wildcards.sample}"
    shell:
        "source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh; conda activate flair; flair correct -o {params.out} -t {threads} -g {params.reffile} -f {params.gtf} -q {input} 2> {log}"

rule collapse:
    input:
        file1="Sample_{sample}/{sample}_all_corrected.bed",
        file2="Sample_{sample}/{sample}.fastq.gz"
    output: "Sample_{sample}/{sample}.isoforms.fa", "Sample_{sample}/{sample}.isoforms.gtf", "Sample_{sample}/{sample}.isoform.read.map.txt"
    resources: mem_mb=400000, time="96:00:00", partition="norm"
    threads: 36
    log: "Sample_{sample}/{sample}_collapse.log"
    params:
        gtf=config["genomes"][ref]["gtf"], reffile=config["genomes"][ref]["ref"],
        out=lambda wildcards: f"Sample_{wildcards.sample}/{wildcards.sample}"
    shell:
        "source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh; conda activate flair; flair collapse -f {params.gtf} -g {params.reffile} -r {input.file2} -q {input.file1} -o {params.out} -t {threads} --check_splice --generate_map --keep_intermediate 2> {log}"

rule quantify:
    input: "Sample_{sample}/{sample}.isoforms.fa",
    output: "Sample_{sample}/flair_quantify.counts.tsv"
    resources: mem_mb=400000, time="96:00:00", partition="norm"
    threads: 16
    log: "Sample_{sample}/{sample}_quantify.log"
    params:
        folder="Sample_{sample}", prefix="Sample_{sample}/{sample}"
    shell:
        """
        echo -e 'thesample\\ta\\tb\\t{params.prefix}.fastq.gz' > {params.folder}/samplesheet.tsv
        source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh
        conda activate flair
        flair quantify -r {params.folder}/samplesheet.tsv -i {input} -t {threads} --temp_dir {params.folder} -o {params.folder}/flair_quantify 2> {log}
        """

rule sqanti:
    input: "Sample_{sample}/{sample}.isoforms.gtf"
    output:
        html="Sample_{sample}/{sample}.sqanti_SQANTI3_report.html",
        txt="Sample_{sample}/{sample}.sqanti_classification.txt",
        fasta="Sample_{sample}/{sample}.sqanti_corrected.fasta"
    resources: mem_mb=400000, time="96:00:00", partition="norm"
    threads: 8
    log: "Sample_{sample}/{sample}.sqanti.log"
    params:
        out=lambda wildcards: f"{wildcards.sample}.sqanti",
        folder=lambda wildcards: f"Sample_{wildcards.sample}",
        gtf=config["genomes"][ref]["gtf"],
        reffile=config["genomes"][ref]["ref"]
    shell:
        "source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh; conda activate SQANTI3.env; export PYTHONPATH=/mnt/nasapps/development/cDNA-Cupcake/29.0.0/sequence; export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3; python /mnt/ccrsf-ifx/Software/tools/SQANTI3/SQANTI3-5.2.1/sqanti3_qc.py {input} {params.gtf} {params.reffile} --aligner_choice=minimap2 -t {threads} -o {params.out} -d {params.folder} --isoAnnotLite --force_id_ignore 2> {log}"

rule merge:
    input: 
        classification="Sample_{sample}/{sample}.sqanti_classification.txt", 
        counts="Sample_{sample}/flair_quantify.counts.tsv",
        readmap="Sample_{sample}/{sample}.isoform.read.map.txt",
        readtag="Sample_{sample}/{sample}/read_tags.tsv"
    output: "Sample_{sample}/{sample}.sqanti_classification_merge.txt"
    resources: mem_mb=400000, time="96:00:00", partition="norm"
    threads: 8
    log: "Sample_{sample}/{sample}.merge.log"
    shell: "python {pyscript}/merge.py -i {input.classification} -c {input.counts} -m {input.readmap} -r {input.readtag} -o {output} 2>{log}"

rule sqanti_filter:
    input:
        classification="Sample_{sample}/{sample}.sqanti_classification_merge.txt"
    output:
        "Sample_{sample}/{sample}.sqanti_RulesFilter_result_classification.txt"
    resources: mem_mb=400000, time="96:00:00", partition="norm"
    threads: 8
    log: "Sample_{sample}/{sample}_sqanti_filter.log"
    params:
        out=lambda wildcards: f"{wildcards.sample}.sqanti",
        folder=lambda wildcards: f"Sample_{wildcards.sample}"
    shell:
        "source /mnt/ccrsf-ifx/Software/tools/Anaconda/3.11/etc/profile.d/conda.sh; conda activate SQANTI3.env; export PYTHONPATH=/mnt/nasapps/development/cDNA-Cupcake/29.0.0/sequence; export R_LIBS=/mnt/nasapps/production/rstudio-ide/2023.03/library:/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3; python /mnt/ccrsf-ifx/Software/tools/SQANTI3/SQANTI3-5.2.1/sqanti3_filter.py rules {input.classification} -o {params.out} -d {params.folder} 2> {log}"                           	