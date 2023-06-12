process SMARTPCA {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bedFile)
    tuple val(meta), path(bimFile)
    tuple val(meta), path(pedIndFile)
    tuple val(meta), path(parFile)
    file popList

    output:
    tuple val(meta), path("*.evec"),   emit: evec
    tuple val(meta), path("*.eval"),   emit: eval
    tuple val(meta), path("*.log"),    emit: log
    path "versions.yml"             ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    smartpca -p $parFile > ${prefix}.smartpca.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smartpca: \$(echo \$(smartpca 2>&1) | head -n 1 | grep -o "smartpca version: [0-9.]*" | sed 's/smartpca version: //')
    END_VERSIONS
    """
}