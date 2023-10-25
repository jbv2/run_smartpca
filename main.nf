#!/usr/bin/env nextflow

/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'run_smartpca' - A Nextflow pipeline to run smartpca from several inputs
 v0.0.2
 October 2023
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Judith Ballesteros Villascán
 GitHub: https://github.com/jbv2/run_smartpca
 ----------------------------------------------------------------------------------------
 */

/* 
 Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// 
// MODULES: Consisting of local modules
//

//
// SUBWORKFLOW: Consisting of local subworkflows
//

// TODO rename to active: index_reference, filter_bam etc.
include { VCF2SMARTPCA   } from './subworkflows/local/vcf2smartpca'
include { PLINK2SMARTPCA } from './subworkflows/local/plink2smartpca'
include { EIGEN2SMARTPCA } from './subworkflows/local/eigen2smartpca'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {

    // Reading poplist according option lsqproject
    if ( params.poplist ) {
      ch_poplist = Channel.fromPath(params.poplist, checkIfExists: true )
        .collect()

    } else {
      ch_poplist = []
    }

    // Read inputs (VCF) and define name as ID

    if (params.input_type == 'vcf' ) {

    ch_vcf = Channel.fromFilePairs(params.inputVCF, size: -1)
        .map {
            meta, vcf ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            [ fmeta, vcf ]
        }

    Channel
        .fromPath(params.samples)
        .set { ch_samples }

    VCF2SMARTPCA(ch_vcf, ch_samples, ch_poplist)
    } else if (params.input_type == 'plink' ) {

    // Read inputs (PLINK) and define name as ID
   ch_bed = Channel.fromFilePairs(params.inputbed, size: -1)
    .map {
        meta, bed ->
        def baseName = bed.baseName.first()
        def dirPath = bed.parent.first()
        def bimPath = "${dirPath}/${baseName}.bim"
        def famPath = "${dirPath}/${baseName}.fam"
        def fmeta = [:]
        // Set meta.id
        fmeta.id = baseName
        [fmeta, bed, file(bimPath), file(famPath)]
    }

    PLINK2SMARTPCA(ch_bed, ch_poplist)

    } else if (params.input_type == 'eigenstrat' ) {

    // Read inputs (EIGENSTRAT) and define name as ID
   ch_geno = Channel.fromFilePairs(params.inputgeno, size: -1)
    .map {
        meta, geno ->
        def baseName = geno.baseName.first()
        def dirPath = geno.parent.first()
        def snpPath = "${dirPath}/${baseName}.snp"
        def indPath = "${dirPath}/${baseName}.ind"
        def fmeta = [:]
        // Set meta.id
        fmeta.id = baseName
        [fmeta, geno, file(snpPath), file(indPath)]
    }

    EIGEN2SMARTPCA(ch_geno, ch_poplist)

    }
}