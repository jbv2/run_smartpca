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

params.r_scripts = "$baseDir/bin/tagger.R"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_REHEADER } from './modules/local/bcftools/reheader/main'
include { PLINK       } from './modules/local/plink/main'
include { PLINK2EIGEN } from './modules/local/plink2eigen/main'
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

    // Read inputs (VCF) and define name as ID
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

    // Change sample names to set FamID before converting to PLINK
    BCFTOOLS_REHEADER(ch_vcf, ch_samples)
    ch_renamed_vcf = BCFTOOLS_REHEADER.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    // Convert VCF to PLINK 
    PLINK(ch_renamed_vcf) 
    ch_versions = ch_versions.mix(PLINK.out.versions)

    // Convert PLINK to EIGENSTRAT
    PLINK2EIGEN(PLINK.out.bed, PLINK.out.bim, PLINK.out.fam)
    ch_versions = ch_versions.mix(PLINK2EIGEN.out.versions)

    // Make PARFILE   
    Channel
        .fromPath(params.poplist)
        .set { ch_poplist }

    ch_parfile = MAKE_PAR(PLINK2EIGEN.out.geno, PLINK2EIGEN.out.snp, PLINK2EIGEN.out.ind, ch_poplist)

    // Run SMARTPCA

    SMARTPCA(PLINK2EIGEN.out.geno, PLINK2EIGEN.out.snp, PLINK2EIGEN.out.ind, MAKE_PAR.out.txt, ch_poplist)
    ch_versions = ch_versions.mix(SMARTPCA.out.versions)

    FORMAT_EVEC(SMARTPCA.out.evec)
}