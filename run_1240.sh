nextflow \
run main.nf \
-profile eva,archgen \
--inputVCF /mnt/archgen/Americas_capture/analyses/new_merge/america_1240_lowcov/america_lowvoc.vcf.gz \
--outdir /mnt/archgen/Americas_capture/analyses/new_merge/smartpca_1240lowcov_projected/ \
--samples /mnt/archgen/Americas_capture/analyses/new_merge/samples2.txt \
--poplist /mnt/archgen/Americas_capture/analyses/new_merge/individuals_components_noks_noadmixed.poplist -resume
