nextflow \
run main.nf \
-profile eva,archgen \
--inputVCF /mnt/archgen/Americas_capture/analyses/new_merge/merge_nobatch_allqual30/snps_no_kar_sur_pim_atransv/america_6M_tomafprune.vcf.gz \
--outdir /mnt/archgen/Americas_capture/analyses/new_merge/smartpca_6M_projected/ \
--samples /mnt/archgen/Americas_capture/analyses/new_merge/samples2.txt \
--poplist /mnt/archgen/Americas_capture/analyses/new_merge/individuals_components_noks_noadmixed.poplist -resume
