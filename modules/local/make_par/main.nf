process MAKE_PAR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(geno)
    tuple val(meta), path(snp)
    tuple val(meta), path(ind)
    file poplistFile
    

    output:
    tuple val(meta), path("*.parfile.txt"),   emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def lsqprojectOption = task.ext.lsqproject ?: "${params.project}"
    def poplistOption = "${params.project}" == 'YES' ? "${poplistFile}" : ''

    """
    echo "genotypename: $geno" > ${prefix}.parfile.txt
    echo "snpname: $snp" >> ${prefix}.parfile.txt 
    echo "indivname: $ind" >> ${prefix}.parfile.txt
    echo "evecoutname: ${prefix}.evec" >> ${prefix}.parfile.txt
    echo "evaloutname: ${prefix}.eval" >> ${prefix}.parfile.txt 
    echo "lsqproject: $lsqprojectOption" >> ${prefix}.parfile.txt
    echo "numoutlieriter: 0" >> ${prefix}.parfile.txt
    echo "shrinkmode: YES" >> ${prefix}.parfile.txt 
    echo "numthreads:   $task.cpus" >> ${prefix}.parfile.txt 
    echo "hiprecision:  YES" >> ${prefix}.parfile.txt 
    echo "fsthiprecision:  YES" >> ${prefix}.parfile.txt 

    # Conditionally set poplistname based on lsqprojectOption
    if [ '$lsqprojectOption' == 'YES' ]; then
        echo "poplistname: $poplistOption" >> ${prefix}.parfile.txt
    fi

    """
}
