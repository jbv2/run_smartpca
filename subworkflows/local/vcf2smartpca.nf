//
// Run SMARTPCA from a VCF as input
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_REHEADER } from '../../modules/local/bcftools/reheader/main'
include { PLINK             } from '../../modules/local/plink/main'
include { PLINK2EIGEN       } from '../../modules/local/plink2eigen/main'
include { MAKE_PAR          } from '../../modules/local/make_par/main'
include { SMARTPCA          } from '../../modules/local/smartpca/main'
include { FORMAT_EVEC       } from '../../modules/local/format_evec/main'

workflow VCF2SMARTPCA {

    take:
    vcf // channel : [ val[meta], [vcf]]
    samples // channel : [ [samples] ]
    poplist // channel : [ [poplist]]

    main:
    ch_versions       = Channel.empty()

    // Change sample names to set FamID before converting to PLINK
    BCFTOOLS_REHEADER(vcf, samples)
    ch_renamed_vcf = BCFTOOLS_REHEADER.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    // Convert VCF to PLINK 
    PLINK(ch_renamed_vcf) 
    ch_versions = ch_versions.mix(PLINK.out.versions)

    // Convert PLINK to EIGENSTRAT
    PLINK2EIGEN(PLINK.out.bed, PLINK.out.bim, PLINK.out.fam)
    ch_versions = ch_versions.mix(PLINK2EIGEN.out.versions)

    // Make PARFILE   

    ch_parfile = MAKE_PAR(PLINK2EIGEN.out.geno, PLINK2EIGEN.out.snp, PLINK2EIGEN.out.ind, poplist)

    // Run SMARTPCA

    SMARTPCA(PLINK2EIGEN.out.geno, PLINK2EIGEN.out.snp, PLINK2EIGEN.out.ind, MAKE_PAR.out.txt, poplist)
    ch_versions = ch_versions.mix(SMARTPCA.out.versions)

    // Format evec
    FORMAT_EVEC(SMARTPCA.out.evec)

    emit:
    format_evec = FORMAT_EVEC.out
    eval        = SMARTPCA.out.eval
    versions    = ch_versions

}