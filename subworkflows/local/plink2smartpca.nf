//
// Run SMARTPCA from a bfile as input
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PLINK2EIGEN       } from '../../modules/local/plink2eigen/main'
include { MAKE_PAR          } from '../../modules/local/make_par/main'
include { SMARTPCA          } from '../../modules/local/smartpca/main'
include { FORMAT_EVEC       } from '../../modules/local/format_evec/main'

workflow PLINK2SMARTPCA {

    take:
    bfile // channel : [ val[meta], [bed], [bim], [fam]]
    poplist // channel : [ [poplist]]

    main:
    ch_versions       = Channel.empty()

    
    // Convert PLINK to EIGENSTRAT
    ch_input = bfile
    .multiMap {meta, bed, bim, fam ->
    bed: [meta, bed]
    bim: [meta, bim]
    fam: [meta, fam]
    }
    PLINK2EIGEN(ch_input)
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