// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2

// Declare Profiles for Conda using an environment file
profiles {
  conda {
    process.conda = "${baseDir}/sganesh/ni_signaling/RNA_environment.yml"
  }
}

// Define input parameters
params.ref_genome = "${baseDir}/shared/refGenomes/hg38/v40/GRCh38.p13.genome.fa"
params.reads = "${baseDir}/sganesh/Labonte/raw_data/*.fastq"
params.results = "${baseDir}/sganesh/Labonte/results"

Channel.fromPath(params.reads).into {read_files}

process fastqc_analysis {
    tag "FASTQC"

    input:
    file reads from read_files

    output:
    path("${params.results}/fastqc_results") into fastqc_results

    script:

    """
    mkdir -p ${params.results}/fastqc_results

    fastqc -q -o fastqc_results $reads
    """
}


process multiqc_analysis {
    tag "QC Reports"

    input:
    file ('fastqc/*') from fastqc_results.collect()

    output:
    path("${params.results}/multiqc_report") into multiqc_reports
    path("${params.results}/multiqc_data")

    script:

    '''
    mkdir -p ${params.results}/multiqc_reports

    multiqc -o multiqc_reports fastqc_results
    '''
}

process index_ref_genome {
    tag "Index Reference Genome for STAR"

    input:
    path ref_genome

    output:
    path 'genome_dir'

    script: 
    """
    mkdir -p genome_dir

    STAR --runMode genomeGenerate \
         --genomeDir genome_dir \
         --genomeFastaFiles ${ref_genome} \
         --runThreadN ${task.cpus}
    """
}

process alignment {
    tag "Alignment"

    input:
    path ref_genome
    path genomeDir
    file reads from read_files


    output:

    path("${params.results}/alignment_logs") into alignment_logs
    path("${params.results}/alignment_bams") into alignment_bams

    //path('*Log.final.out')   , emit: log_final
    //path('*Log.out')         , emit: log_out
    //path('*Log.progress.out'), emit: log_progress

    //path('*d.out.bam')              , optional:true, emit: bam
    //path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted

    script:
    
    """

    # Perform alignment and generate results in BAM format
    STAR --genomeDir ${genomeDir} \
         --readFilesIn ${reads} \
         --runThreadN ${task.cpus} \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate


    samtools index Aligned.sortedByCoord.out.bam
    """
}

// Define workflow execution
workflow {

    // Index the reference Genome for STAR
    index_ref_genome(params.ref_genome)

    // Quality control
    fastqc_analysis(read_files)
    multiqc_analysis(fastqc_results)

    // Alignment
    alignment(params.ref_genome, index_ref_genome.out, read_files)
}
