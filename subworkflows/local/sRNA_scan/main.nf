nextflow.enable.dsl = 2

/*
 * MODULE IMPORTS
 */

WD = new File(".").getCanonicalPath()


include { BOWTIE_BUILD as RRNA_BUILD } from '../../../modules/nf-core/bowtie/build/main.nf'
include { BOWTIE_BUILD as TRNA_BUILD } from '../../../modules/nf-core/bowtie/build/main.nf'
include { BOWTIE_BUILD as NCRNA_BUILD } from '../../../modules/nf-core/bowtie/build/main.nf'

include { BOWTIE_ALIGN as RRNA_ALIGN } from '../../../modules/nf-core/bowtie/align/main.nf'
include { BOWTIE_ALIGN as TRNA_ALIGN } from '../../../modules/nf-core/bowtie/align/main.nf'
include { BOWTIE_ALIGN as NCRNA_ALIGN } from '../../../modules/nf-core/bowtie/align/main.nf'

include { SAMTOOLS_SORT as RRNA_SORT } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as TRNA_SORT } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as NCRNA_SORT } from '../../../modules/nf-core/samtools/sort/main.nf'

include { SAMTOOLS_INDEX as RRNA_INDEX } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as TRNA_INDEX } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as NCRNA_INDEX } from '../../../modules/nf-core/samtools/index/main.nf'

include { SUBREAD_FEATURECOUNTS as RRNA_COUNT } from '../../../modules/nf-core/subread/featurecounts/main.nf'
include { SUBREAD_FEATURECOUNTS as TRNA_COUNT } from '../../../modules/nf-core/subread/featurecounts/main.nf'
include { SUBREAD_FEATURECOUNTS as NCRNA_COUNT } from '../../../modules/nf-core/subread/featurecounts/main.nf'

include { UMITOOLS_DEDUP as RRNA_DEDUP } from '../../../modules/nf-core/umitools/dedup/main.nf'
include { UMITOOLS_DEDUP as TRNA_DEDUP } from '../../../modules/nf-core/umitools/dedup/main.nf'
include { UMITOOLS_DEDUP as NCRNA_DEDUP } from '../../../modules/nf-core/umitools/dedup/main.nf'


/*
 * WORKFLOW
 */
workflow CONTAMINANT_FILTER_SRNAS {

    take:
        reads
        rrna_fa
        trna_fa
        other_fa
        rrna_gtf
        trna_gtf
        other_gtf

    main:

    Channel.empty().set { versions }
    ch_clean = reads

    /*
     * === RRNA FILTER ===
     */
    INDEX_RRNA = RRNA_BUILD(
        fasta: rrna_fa
    )

    ALIGN_RRNA = RRNA_ALIGN(
        reads: ch_clean,
        index: INDEX_RRNA.index,
        meta: true
    )

    SORT_RRNA = RRNA_SORT(bam: ALIGN_RRNA.bam)
    IDX_RRNA  = RRNA_INDEX(bam: SORT_RRNA.bam)

    COUNT_RRNA = RRNA_COUNT(
        bam: SORT_RRNA.bam,
        gtf: rrna_gtf
    )

    versions = versions.mix(INDEX_RRNA.versions, ALIGN_RRNA.versions)

    ch_clean = ALIGN_RRNA.unmapped_fastq


    /*
     * === TRNA FILTER ===
     */
    INDEX_TRNA = TRNA_BUILD(
        fasta: trna_fa
    )

    ALIGN_TRNA = TRNA_ALIGN(
        reads: ch_clean,
        index: INDEX_TRNA.index,
        meta: true
    )

    SORT_TRNA = TRNA_SORT(bam: ALIGN_TRNA.bam)
    IDX_TRNA  = TRNA_INDEX(bam: SORT_TRNA.bam)

    COUNT_TRNA = TRNA_COUNT(
        bam: SORT_TRNA.bam,
        gtf: trna_gtf
    )

    versions = versions.mix(INDEX_TRNA.versions, ALIGN_TRNA.versions)

    ch_clean = ALIGN_TRNA.unmapped_fastq


    /*
     * === OTHER ncRNA FILTER ===
     */
    INDEX_OTHER = NCRNA_BUILD(
        fasta: other_fa
    )

    ALIGN_OTHER = NCRNA_ALIGN(
        reads: ch_clean,
        index: INDEX_OTHER.index,
        meta: true
    )

    SORT_OTHER = NCRNA_SORT(bam: ALIGN_OTHER.bam)
    IDX_OTHER  = NCRNA_INDEX(bam: SORT_OTHER.bam)

    COUNT_OTHER = NCRNA_COUNT(
        bam: SORT_OTHER.bam,
        gtf: other_gtf
    )

    versions = versions.mix(INDEX_OTHER.versions, ALIGN_OTHER.versions)

    ch_clean = ALIGN_OTHER.unmapped_fastq


    emit:
        clean_reads  = ch_clean

        rrna_bam     = SORT_RRNA.bam
        trna_bam     = SORT_TRNA.bam
        other_bam    = SORT_OTHER.bam

        rrna_counts  = COUNT_RRNA.counts
        trna_counts  = COUNT_TRNA.counts
        other_counts = COUNT_OTHER.counts


        versions     = versions
}
