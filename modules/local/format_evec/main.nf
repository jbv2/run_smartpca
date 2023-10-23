process FORMAT_EVEC {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(evec)

    output:
    tuple val(meta), path("*.formated_evec")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $evec \
    | tr -s " " \
    | tr -d "#" \
    | sed "s#^ ##" \
    | sed "s# #\t#g" > tmp

    echo "sample PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 tag" | tr " " "\t" > ${prefix}.formated_evec
    tail -n+2 tmp | sed "s#America:##" >> ${prefix}.formated_evec
    """
}
