# Defining a global wildcard for samples to be used
SAMPLES, = glob_wildcards("samples/{sample}.fastq.gz")

# Full genome for hisat2 alignment
full_Ref="samples/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Annotation for full reference required
ref_annot="samples/Homo_sapiens.GRCh38.100.gtf"

# Transcriptome index required by kallisto
txpme="samples/gencode.v34.transcripts.fa"

#Defining the sample read and other reference files used
r1="samples/{sample}.fastq.gz"


#Defining the input files for the final rule, which should be the final files we want to generate
rule all:
        input:
#                "hisat2/hisat2_counts.txt",
                "Results/Fastqc_Reports/multiqc_report.html",
                "Results/Trim_galore/multiqc_report.html",
#		"hisat2/multiqc_report.html",
#                expand("Kallisto/{sample}/abundance.tsv",  sample=SAMPLES)
                expand("hisat2/{sample}.bam", sample=SAMPLES)

#Quality Check
rule fastqc_check1:
        input:
                read=r1
        output:
                "Results/Fastqc_Reports/{sample}_fastqc.html"
        params:
                "Results/Fastqc_Reports"
        shell:
                "fastqc {input.read} -o {params}"

#Multiqc on the initial Fastqc reports
rule multiqc1:
        input:
                html=expand("Results/Fastqc_Reports/{sample}_fastqc.html", sample=SAMPLES)
        params:
                dir="Results/Fastqc_Reports"
        output:
                "Results/Fastqc_Reports/multiqc_report.html"
        shell:
                "multiqc {params.dir}  -o {params.dir}"


#Trimming of samples with Trim_galore
rule trimming:
        input:
                read=r1,
                trim_html="Results/Fastqc_Reports/{sample}_fastqc.html"
        output:
                r1="Results/Trim_galore/{sample}_trimmed.fq.gz",
                r1_html="Results/Trim_galore/{sample}_trimmed_fastqc.html"
#                touch("Results/Trim_galore/file.txt")
        params:
                dir="Results/Trim_galore"
        shell:
                "trim_galore -j 8  {input.read} -q 20 --length 20 --fastqc -o {params.dir}"


#Multiqc results for the trimmed fastq reads
rule multiqc2:
        input:
                trim_html=expand("Results/Trim_galore/{sample}_trimmed_fastqc.html", sample=SAMPLES)
        output:
                "Results/Trim_galore/multiqc_report.html"
        params:
                dir="Results/Trim_galore"
        shell:
                "multiqc {params.dir} -o {params.dir}"


#Indexing reference with HISAT2
rule hisat2_indexing:
        input:
                ref=full_Ref
        output:
                touch("hisat2/makeidx.done")
        params:
                threads=40,
                idx="hisat2/Homo_sapiens.GRCh38v3_hisat2.idx"
        shell:
                "hisat2-build -p {params.threads} {input.ref} {params.idx}"

#Alignment using HISAT2
rule hisat2_Alignment:
        input:
                idxdone="hisat2/makeidx.done",
                read="Results/Trim_galore/{sample}_trimmed.fq.gz"
        output:
                "hisat2/{sample}.bam"
        params:
                idx="hisat2/Homo_sapiens.GRCh38v3_hisat2.idx",
                threads=45,
                bam="hisat2/{sample}.sam"
        shell:
                """
                hisat2 -p {params.threads} -x {params.idx}  -U {input.read}  -S {params.bam}
                samtools view -@ {params.threads} -Sbh {params.bam} | samtools sort -@ {params.threads} > {output}; samtools index {output}; rm {params.bam}
                """


#
