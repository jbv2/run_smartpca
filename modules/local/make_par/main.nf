process MAKE_PAR {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(famFile)
    file samplesFile
    file poplistFile

    output:
    tuple val(meta), path("*.parfile.txt"),   emit: txt
    tuple val(meta), path("*.pedind"),        emit: pedind
    path "versions.yml"             ,         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript --vanilla ${params.r_scripts} $famFile ${samplesFile} ${prefix}.pedind

    echo "genotypename: ${prefix}.bed" > ${prefix}.parfile.txt
    echo "snpname: ${prefix}.bim" >> ${prefix}.parfile.txt 
    echo "indivname: ${prefix}.pedind" >> ${prefix}.parfile.txt
    echo "evecoutname: ${prefix}.evec" >> ${prefix}.parfile.txt
    echo "evaloutname: ${prefix}.eval" >> ${prefix}.parfile.txt 
    echo "poplistname:  ${poplistFile}" >> ${prefix}.parfile.txt
    echo "lsqproject: YES" >> ${prefix}.parfile.txt
    echo "numoutlieriter: 0" >> ${prefix}.parfile.txt
    echo "shrinkmode: YES" >> ${prefix}.parfile.txt 
    echo "numthreads:   4" >> ${prefix}.parfile.txt 
    echo "hiprecision:  YES" >> ${prefix}.parfile.txt 
    echo "fsthiprecision:  YES" >> ${prefix}.parfile.txt  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version) | sed 's/^R version//;s/(20.*//')
    END_VERSIONS
    """
}