#!/usr/bin/env nextflow

/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'nf-run-smartpca' - A Nextflow pipeline to run smartpca
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Judith Ballesteros Villascán
 GitHub: https://github.com/jbv2/run-smartPCA
 ----------------------------------------------------------------------------------------
 */

/* 
 Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */ 

params.r_scripts = "/Users/judith_ballesteros/Documents/Ongoing_projects/run_smartpca/bin/tagger.R"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PLINK       } from './modules/local/plink/main'
include { MAKE_PAR    } from './modules/local/make_par/main'
include { SMARTPCA    } from './modules/local/smartpca/main'
include { FORMAT_EVEC } from './modules/local/format_evec/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


workflow {

    ch_versions       = Channel.empty()

    ch_input = Channel.fromFilePairs(params.inputVCF, size: -1)
        .map {
            meta, vcf ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            [ fmeta, vcf ]
        }

    ch_plink = PLINK(ch_input) 
    ch_versions = ch_versions.mix(PLINK.out.versions)

    Channel
        .fromPath(params.samples)
        .set { ch_samples }

    ch_parfile = MAKE_PAR(PLINK.out.fam, ch_samples)
    ch_versions = ch_versions.mix(MAKE_PAR.out.versions)

    Channel
        .fromPath(params.poplist)
        .set { ch_poplist }

    SMARTPCA(PLINK.out.bed, PLINK.out.bim, MAKE_PAR.out.pedind, MAKE_PAR.out.txt, ch_poplist)
    ch_versions = ch_versions.mix(SMARTPCA.out.versions)

    FORMAT_EVEC(SMARTPCA.out.evec)
}