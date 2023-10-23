process PLINK2EIGEN {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bed)
    tuple val(meta), path(bim)
    tuple val(meta), path(fam)

    output:
    tuple val(meta), path("*.geno"),  emit: geno
    tuple val(meta), path("*.snp"),   emit: snp
    tuple val(meta), path("*.ind"),   emit: ind
    path "versions.yml"           ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    trident genoconvert \\
        --genoFile $geno \\
        --snpFile $snp \\
        --indFile $ind \\
        --inFormat PLINK \\
        --outFormat EIGENSTRAT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trident: \$(echo \$(trident --version ) )
    END_VERSIONS
    """
}