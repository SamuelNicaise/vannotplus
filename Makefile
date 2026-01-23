install:
	pip install -e .

barcode:
	python -m vannotplus barcode -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/trio_base.vcf -o /home1/HUB/bin/vannotplus/vannotplus/output.vcf -c src/vannotplus/config.yml -a WES_AGILENT

exomiser:
	python -m vannotplus exomiser -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/jb_trio_output.vcf -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/test/merged.vcf -c src/vannotplus/config.yml -a WES_AGILENT -v debug

exomiser_jean:
	python -m vannotplus exomiser -i /home1/HUB/bin/vannotplus/vannotplus/jean/merged.vcf.gz -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/test/merged_jean.vcf -c src/vannotplus/config.yml -a JEAN -v debug

annot:
	python -m vannotplus annot -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/SGT151797.final.clean.full-annotation.v4.noprio.vcf.gz -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/annot.vcf -c src/vannotplus/config.yml -v debug

annot_score:
	python -m vannotplus annot -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/SGT151797.final.clean.full-annotation.v4.noprio.vcf.gz -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/annot.vcf -c src/vannotplus/config.yml -v debug --vannotscore

annot_wes:
	python -m vannotplus annot -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/KLA2403985.final.vcf -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/annotKLA.vcf -c src/vannotplus/config.yml

debug:
	python -m vannotplus annot -i SGT151797.final.clean.full-annotation.v4.noprio.vcf.gz -o annot.vcf -c src/vannotplus/config.yml -v debug --vannotscore

send:
	rm -rf /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/vannotplus
	cp -r /home1/HUB/bin/vannotplus/vannotplus /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises

gmc_test:
	python -m vannotplus annot -i ./VANNOT_merged_GENODENT.design.vcf.gz -o genodent_test_gmc.vcf -c src/vannotplus/config.yml

gmc_mini:
	python -m vannotplus annot -i ./gmc_mini_input.vcf -o gmc_mini_output.vcf -c src/vannotplus/config.yml


gmc_test_routine:
	python -m vannotplus annot -i ./vannot.exome718.vcf.gz -o exome718_test_gmc.vcf -c src/vannotplus/config.yml

gmc_test_genome:
	python -m vannotplus annot -i ./GENODYT_1.filtered.split-multi.vcf.gz -o genodyt_gmc.vcf -c src/vannotplus/config.yml

gmc_test_single:
	python -m vannotplus annot -i ./gmc_single_input.vcf -o single_gmc.vcf -c src/vannotplus/config.yml

audrey:
	python -m vannotplus exomiser -i /home1/data/WORK_DIR_SAM/pmda/merged_28samples.vcf.gz -o /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz -c src/vannotplus/config.yml -a AUDREY_1 -v debug
	/home1/HUB/bin/bcftools/1.14/bin/bcftools annotate -c FORMAT/EXOMISER_GENE_PHENO_SCORE:=FORMAT/EXOMISER_1_PATHOLOGIE -O z -o /home1/L_PROD/NGS/tmp/pmda_audrey1_pathologie.vcf.gz /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz
	rm -f /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz
	python -m vannotplus exomiser -i /home1/data/WORK_DIR_SAM/pmda/merged_28samples.vcf.gz -o /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz -c src/vannotplus/config.yml -a AUDREY_2 -v debug
	/home1/HUB/bin/bcftools/1.14/bin/bcftools annotate -c FORMAT/EXOMISER_GENE_PHENO_SCORE:=FORMAT/EXOMISER_2_PRESCRIPTION -O z -o /home1/L_PROD/NGS/tmp/pmda_audrey2_prescription.vcf.gz /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz
	rm -f /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz
	python -m vannotplus exomiser -i /home1/data/WORK_DIR_SAM/pmda/merged_28samples.vcf.gz -o /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz -c src/vannotplus/config.yml -a AUDREY_3 -v debug
	/home1/HUB/bin/bcftools/1.14/bin/bcftools annotate -c FORMAT/EXOMISER_GENE_PHENO_SCORE:=FORMAT/EXOMISER_3_COMPTERENDU -O z -o /home1/L_PROD/NGS/tmp/pmda_audrey3_compterendu.vcf.gz /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz
	rm -f /home1/L_PROD/NGS/tmp/pmda_audrey_tmp.vcf.gz
	for v in /home1/L_PROD/NGS/tmp/pmda_audrey*.vcf.gz; do /home1/HUB/bin/htslib/1.14/bin/tabix $$v; done
	/home1/HUB/bin/bcftools/1.14/bin/bcftools merge -m none -O z -o /home1/L_PROD/NGS/tmp/pmda_audrey_merged.vcf.gz /home1/L_PROD/NGS/tmp/pmda_audrey1_pathologie.vcf.gz /home1/L_PROD/NGS/tmp/pmda_audrey2_prescription.vcf.gz /home1/L_PROD/NGS/tmp/pmda_audrey3_compterendu.vcf.gz
