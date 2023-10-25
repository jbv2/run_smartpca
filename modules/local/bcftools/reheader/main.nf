process BCFTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)
    path(samples_names)

    output:
    tuple val(meta), path("*_renamed.vcf.gz"), emit: vcf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools \\
        reheader \\
        --samples $samples_names \\
        $args \\
        --threads $task.cpus \\
        $vcf \\
        | bcftools view \\
        -Oz \\
        --output ${prefix}_renamed.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

}