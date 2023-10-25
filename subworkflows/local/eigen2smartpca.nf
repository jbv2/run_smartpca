//
// Run SMARTPCA from a bfile as input
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAKE_PAR          } from '../../modules/local/make_par/main'
include { SMARTPCA          } from '../../modules/local/smartpca/main'
include { FORMAT_EVEC       } from '../../modules/local/format_evec/main'

workflow EIGEN2SMARTPCA {

    take:
    eigenstrat // channel : [ val[meta], [geno], [snp], [ind]]
    poplist // channel : [ [poplist]]

    main:
    ch_versions       = Channel.empty()

    
    // Read EIGENSTRAT input
    ch_input = eigenstrat
    .multiMap {meta, geno, snp, ind ->
    geno: [meta, geno]
    snp: [meta, snp]
    ind: [meta, ind]
    }

    // Make PARFILE   

    ch_parfile = MAKE_PAR(ch_input, poplist)

    // Run SMARTPCA

    SMARTPCA(ch_input, MAKE_PAR.out.txt, poplist)
    ch_versions = ch_versions.mix(SMARTPCA.out.versions)

    // Format evec
    FORMAT_EVEC(SMARTPCA.out.evec)

    emit:
    format_evec = FORMAT_EVEC.out
    eval        = SMARTPCA.out.eval
    versions    = ch_versions

}