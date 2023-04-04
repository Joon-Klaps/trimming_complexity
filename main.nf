#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INPUT_CHECK                 } from "./workflows/input_check.nf"
include { CUTADAPT                    } from './modules/nf-core/cutadapt/main'                                                                                                                    
include { FASTP                       } from './modules/nf-core/fastp/main'                                                                                                                          
include { ADAPTERREMOVAL              } from './modules/nf-core/adapterremoval/main'                                                                                                        
include { FAQCS                       } from './modules/nf-core/faqcs/main'                                                                                                                          
include { SICKLE                      } from './modules/nf-core/sickle/main'                                                                                                                        
include { TRIMGALORE                  } from './modules/nf-core/trimgalore/main'
include { TRIMMOMATIC                 } from './modules/nf-core/trimmomatic/main'                                                                                                                
include { BBMAP_BBDUK                 } from './modules/nf-core/bbmap/bbduk/main'       
include { PRINSEQPLUSPLUS             } from './modules/nf-core/prinseqplusplus/main'                                                                                                        
include { FASTQC  as FASTQC_RAW       } from './modules/nf-core/fastqc/main'
include { FASTQC  as FASTQC_TRIMMED   } from './modules/nf-core/fastqc/main' 
include { BOWTIE2_ALIGN               } from './modules/nf-core/bowtie2/align/main' 
include { SPADES                      } from './modules/nf-core/spades/main'
include { QUAST                       } from './modules/nf-core/quast/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'  
include { MULTIQC                     } from './modules/nf-core/multiqc/main'                                                                                                                                                                                                                                                                                    

if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

workflow TRIMMING_BENCH {
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()
    
    adapter_fasta = params.adapter_fasta
    adapterlist = params.adapter_tsv
    INPUT_CHECK (
        ch_input
    )
    reads_raw= INPUT_CHECK.out.reads

    FASTQC_RAW( reads_raw )
        
    qual_type=Channel.of('sanger')

    CUTADAPT(       reads_raw)                      // no adapter list in input is in config
    FASTP(          reads_raw,adapter_fasta,false, false)        
    ADAPTERREMOVAL( reads_raw,adapterlist)          // adapterlist is text file
    FAQCS(          reads_raw)                      // no adapter list in input is in config
    SICKLE(         reads_raw.combine(qual_type))   // doesn't need adapter file, I'm guessing only quality trimming
    //TRIMGALORE(     reads_raw)                    // Wrapper for cutadapt and fastqc 
    TRIMMOMATIC(    reads_raw)                      // no adapter list in input it is in config
    BBMAP_BBDUK(    reads_raw,adapter_fasta)        
    PRINSEQPLUSPLUS(reads_raw)                      // no adapter list 

    //add all trimmed reads into a single channel 'mixed8_trimmed_reads' and append a name to their id's
    mixed_trimmed_reads = CUTADAPT.out.reads
            .map{rename_meta_id(it, "cutadapt")}
    mixed2_trimmed_reads = mixed_trimmed_reads.mix(FASTP.out.reads
            .map{rename_meta_id(it, "fastp")})
    mixed3_trimmed_reads = mixed2_trimmed_reads.mix(BBMAP_BBDUK.out.reads
            .map{rename_meta_id(it, "bbduk")})
    mixed4_trimmed_reads = mixed3_trimmed_reads.mix(ADAPTERREMOVAL.out.paired_truncated
            .map{rename_meta_id(it, "adapterremoval")})
    mixed5_trimmed_reads = mixed4_trimmed_reads.mix(FAQCS.out.reads
            .map{rename_meta_id(it, "cutadapt")})
    mixed6_trimmed_reads = mixed5_trimmed_reads.mix(SICKLE.out.paired_trimmed
            .map{rename_meta_id(it, "sickle")})
    // trimmed_reads = trimmed_reads.mix(TRIMGALORE.out.reads
    //         .map{rename_meta_id(it,"trimgalore")})
    mixed7_trimmed_reads = mixed6_trimmed_reads.mix(TRIMMOMATIC.out.trimmed_reads
            .map{rename_meta_id(it, "trimmomatic")})
    mixed8_trimmed_reads = mixed7_trimmed_reads.mix(PRINSEQPLUSPLUS.out.good_reads
            .map{rename_meta_id(it, "prinseq++")})
        
    // Run through fastQC again
    FASTQC_TRIMMED(mixed8_trimmed_reads)

    // Extract the indexed reference location of all reads. 
    reference_idx_ch = mixed8_trimmed_reads
        .map{[it[0],it[0].reference_index]}
    
    BOWTIE2_ALIGN(mixed8_trimmed_reads,reference_idx_ch,false,true)
    
    spades_reads = mixed8_trimmed_reads
        .map{ meta, reads  ->
                [meta, reads, [],[]]}

    SPADES(spades_reads, [], [] )
    
    reference_ch = SPADES.out.contigs
        .map{[it[0],it[0].reference]}

    QUAST(SPADES.out.contigs,reference_ch, [],true , false)

    //mix all fastqc files into multiqc of all the fastqc's and the log of them themselves
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_versions = ch_versions.mix(ADAPTERREMOVAL.out.versions.first())
    ch_versions = ch_versions.mix(FAQCS.out.versions.first())
    ch_versions = ch_versions.mix(SICKLE.out.versions.first())
    // ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
    ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    ch_versions = ch_versions.mix(SPADES.out.versions.first())
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL.out.settings.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SICKLE.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PRINSEQPLUSPLUS.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect().ifEmpty([]))
    // multiqc doesn't support FAQCS, TRIMGALORE summary output 
    // ch_multiqc_files = ch_multiqc_files.mix(FAQCS.out.log.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
    

    MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ], ]
def extend_meta(tuple_data, String extension) {
    def meta = tuple_data[0].clone()
    meta.extension = "$extension"
    
    return tuple(meta, tuple_data[1])
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    TRIMMING_BENCH ();
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

